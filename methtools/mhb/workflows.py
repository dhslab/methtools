# methtools/mhb/workflows.py
"""
This module defines high-level workflows for Methylation Haplotype Block (MHB)
analysis, orchestrating calls to the core API and processing steps.
It is intended to be called by the CLI command definitions.
"""
import logging
from pathlib import Path
from typing import List, Optional, Union
import shutil

import pandas as pd
# import numpy as np # No longer directly used in the top-level functions after refactor
from pybedtools import BedTool

from ..utils.system import get_temporary_directory, check_tool_version # Added check_tool_version
from . import api as mhb_api
from .parallel_processor import find_mhb_for_all_regions_parallel

logger = logging.getLogger(__name__)


def _validate_common_inputs(
    regions_bed_file_path: Path,
    cpg_annotations_bed_path: Path,
    output_dir_path: Path
) -> bool:
    """Validates common input files and directories."""
    common_files_to_check = {
        "Regions BED file": regions_bed_file_path,
        "CpG annotations BED file": cpg_annotations_bed_path,
    }
    for name, file_path in common_files_to_check.items():
        if not file_path.is_file():
            logger.error(f"Input Error: {name} not found at '{file_path}'. Workflow cannot proceed.")
            return False
    try:
        output_dir_path.mkdir(parents=True, exist_ok=True)
        logger.info(f"Ensured output directory exists: '{output_dir_path.resolve()}'")
    except OSError as e:
        logger.error(f"Directory Error: Could not create or access output directory '{output_dir_path}': {e}. Workflow cannot proceed.")
        return False
    return True


def _validate_and_prepare_output_formats(output_formats_req: Optional[List[str]]) -> Optional[List[str]]:
    """Validates and prepares the list of output formats."""
    if output_formats_req is None:
        output_formats_req = ["tsv", "bed"]
    
    valid_formats = []
    for fmt in output_formats_req:
        fmt_lower = fmt.strip().lower()
        if fmt_lower in ["tsv", "bed"]:
            valid_formats.append(fmt_lower)
        else:
            logger.warning(f"Configuration Warning: Unsupported output format '{fmt}' specified. It will be ignored.")
    if not valid_formats:
        logger.error("Configuration Error: No valid output formats specified (supported: tsv, bed). Aborting workflow.")
        return None
    logger.info(f"Validated output formats to be generated: {valid_formats}")
    return valid_formats


def _mhb_core_processing_and_saving(
    actual_methylation_calls_files: List[Path],
    regions_bed_file_path: Path,
    cpg_annotations_bed_path: Path,
    output_dir_path: Path,
    num_threads: int,
    output_prefix: str,
    output_formats_validated: List[str],
    mld_pairing_threshold: int,
    block_r2_soft_threshold: float,
    block_r2_hard_threshold: float,
    block_min_cpgs_in_block: int,
    block_tolerance_gaps: int
) -> bool:
    """
    Core MHB processing steps from intersecting regions to saving results.
    This function corresponds to steps 4-8 of the original combined workflow.
    """
    try:
        # --- Step 4: Intersect Analysis Regions with All CpG Positions ---
        logger.info("Step 4: Intersecting analysis regions with CpG annotations to identify relevant CpG sites...")
        regions_df_for_mhb_analysis = pd.read_csv(regions_bed_file_path, sep="\t", header=None, usecols=[0, 1, 2])
        regions_df_for_mhb_analysis.columns = ["chrom", "start", "end"]
        
        regions_bedtool = BedTool.from_dataframe(regions_df_for_mhb_analysis)
        cpg_annotations_bedtool = BedTool(str(cpg_annotations_bed_path))

        intersected_cpgs_df = regions_bedtool.intersect(
            cpg_annotations_bedtool, wao=True
        ).to_dataframe(
            disable_auto_names=True, header=None,
            names=["analysis_region_chrom", "analysis_region_start", "analysis_region_end", 
                   "cpg_chrom", "cpg_start", "cpg_end", "overlap_bp"]
        )
        
        relevant_cpgs_df = intersected_cpgs_df[intersected_cpgs_df["overlap_bp"] > 0].copy()
        relevant_cpgs_df.rename(columns={"cpg_chrom": "chrom", "cpg_start": "cpg_pos"}, inplace=True)
        all_workflow_cpg_positions_df = relevant_cpgs_df[["chrom", "cpg_pos"]].drop_duplicates().reset_index(drop=True)

        if all_workflow_cpg_positions_df.empty:
            logger.error("No CpG sites found within the specified analysis regions after intersection with annotations. Aborting core processing.")
            return False
        logger.info(f"Identified {len(all_workflow_cpg_positions_df)} unique CpG sites within the analysis regions.")

        # --- Step 5: Generate Consolidated Methylation Calls Matrix ---
        logger.info("Step 5: Generating consolidated methylation matrix from call files...")
        calls_df = mhb_api.generate_methylation_matrix(actual_methylation_calls_files)
        if calls_df.empty:
            logger.error("Consolidated methylation matrix is empty after processing call files. Aborting core processing.")
            return False
        logger.info(f"Consolidated methylation matrix created with {calls_df.shape[0]} methylation calls.")

        # --- Step 6: Collapse CpG Positions ---
        logger.info("Step 6: Collapsing CpG positions to a common reference strand...")
        if "ref_position" not in calls_df.columns or "ref_strand" not in calls_df.columns:
            logger.error("Required columns 'ref_position' or 'ref_strand' not found in methylation matrix for collapsing. Aborting core processing.")
            return False
            
        calls_df["ref_position_collapsed"] = calls_df.apply(
            lambda row: mhb_api.collapse_cpg_positions(row["ref_position"], row["ref_strand"]),
            axis=1
        )
        logger.info("CpG positions successfully collapsed.")

        # --- Step 7: Call MHBs in Parallel ---
        logger.info("Step 7: Calling Methylation Haplotype Blocks in parallel across regions...")
        all_blocks_df = find_mhb_for_all_regions_parallel(
            bed_df=regions_df_for_mhb_analysis, 
            calls_df=calls_df,
            all_regions_cpg_positions_df=all_workflow_cpg_positions_df,
            mld_pairing_threshold=mld_pairing_threshold,
            block_r2_soft_threshold=block_r2_soft_threshold,
            block_r2_hard_threshold=block_r2_hard_threshold,
            block_min_cpgs_in_block=block_min_cpgs_in_block,
            block_tolerance_gaps=block_tolerance_gaps,
            num_threads=num_threads
        )

        if all_blocks_df.empty:
            logger.warning("No MHBs were identified in any of the processed regions.")
        else:
            logger.info(f"Successfully identified a total of {all_blocks_df.shape[0]} MHBs across all processed regions.")

        # --- Step 8: Save Results ---
        logger.info("Step 8: Saving results...")
        if not all_blocks_df.empty:
            for fmt in output_formats_validated:
                output_file_name = f"{output_prefix}.mhb.{fmt}"
                output_file_path = output_dir_path / output_file_name
                try:
                    if fmt == "tsv":
                        all_blocks_df.to_csv(output_file_path, sep="\t", index=False)
                        logger.info(f"Saved MHB blocks (TSV) to: '{output_file_path}'")
                    elif fmt == "bed":
                        bed_output_df = all_blocks_df[[
                            "chrm", "first_cpg_pos", "last_cpg_pos", 
                            "block_id", 
                            "num_cpgs"  
                        ]].copy()
                        bed_output_df.rename(columns={"block_id": "name", "num_cpgs": "score"}, inplace=True)
                        bed_output_df.to_csv(output_file_path, sep="\t", index=False, header=False)
                        logger.info(f"Saved MHB blocks (BED) to: '{output_file_path}'")
                except Exception as e:
                    logger.error(f"Failed to save results in '{fmt}' format to '{output_file_path}': {e}", exc_info=True)
        else:
            logger.info("No MHB blocks to save as none were identified.")
        
        return True

    except FileNotFoundError as e:
        logger.error(f"Core Processing Error: Missing file - {e}", exc_info=True)
        return False
    except pd.errors.EmptyDataError as e:
        logger.error(f"Core Processing Error: Empty input file encountered - {e}", exc_info=True)
        return False
    except KeyError as e:
        logger.error(f"Core Processing Error: Missing expected column in data - {e}", exc_info=True)
        return False
    except Exception as e:
        logger.error(f"Core Processing Error: An unexpected error occurred - {e}", exc_info=True)
        return False


def mhb_workflow_from_bams(
    fasta_file_path: Path,
    bam_file_paths: List[Path],
    regions_bed_file_path: Path,
    cpg_annotations_bed_path: Path,
    output_dir_path: Path,
    custom_temp_dir_base: Optional[Path] = None,
    num_threads: int = 1,
    output_prefix: str = "mhb_analysis",
    output_formats: Optional[List[str]] = None,
    mld_pairing_threshold: int = 10,
    block_r2_soft_threshold: float = 0.5,
    block_r2_hard_threshold: float = 0.4,
    block_min_cpgs_in_block: int = 4,
    block_tolerance_gaps: int = 1,
    keep_temp_files: bool = False,
    save_extracted_calls_to: Optional[Path] = None
) -> bool:
    """
    Orchestrates the MHB calling workflow starting from BAM files.
    Includes methylation extraction and prerequisite tool checks.
    """
    logger.info("Starting MHB Workflow (from BAMs)...")
    logger.info(f"Output directory set to: '{output_dir_path.resolve()}'")
    logger.info(f"Number of threads for parallel tasks: {num_threads}")

    output_formats_validated = _validate_and_prepare_output_formats(output_formats)
    if not output_formats_validated:
        return False

    if not fasta_file_path.is_file():
        logger.error(f"Input Error: FASTA reference not found at '{fasta_file_path}'. Workflow cannot proceed.")
        return False
    if not bam_file_paths:
        logger.error("Input Error: No BAM files provided. Workflow cannot proceed.")
        return False
    for bam_path in bam_file_paths:
        if not bam_path.is_file():
            logger.error(f"Input Error: Input BAM file not found: '{bam_path}'. Workflow cannot proceed.")
            return False
            
    if not _validate_common_inputs(regions_bed_file_path, cpg_annotations_bed_path, output_dir_path):
        return False

    # --- Prerequisite Tool Checks ---
    logger.info("Checking for prerequisite command-line tools...")
    try:
        # Define minimum versions (these are examples, adjust as necessary)
        samtools_min_version = "1.9"
        modkit_min_version = "0.1.10" # Version that includes 'extract calls'
        bedtools_min_version = "2.31.1"  # Version required for accurate BED file operations

        logger.info(f"Checking for samtools (minimum version {samtools_min_version})...")
        check_tool_version("samtools", samtools_min_version)

        logger.info(f"Checking for modkit (minimum version {modkit_min_version})...")
        check_tool_version("modkit", modkit_min_version)

        logger.info(f"Checking for bedtools (minimum version {bedtools_min_version})...")
        check_tool_version("bedtools", bedtools_min_version)
        logger.info("Prerequisite tools check passed.")

    except RuntimeError as e:
        # check_tool_version logs details of the failure.
        logger.error(f"Tool version check failed: {e}. Workflow cannot proceed.")
        return False

    main_temp_dir_path: Optional[Path] = None
    actual_methylation_calls_files: List[Path] = []
    user_specified_calls_output_dir_for_cleanup_check: Optional[Path] = None

    try:
        # --- Setup Temporary Directory ---
        logger.info("Setting up temporary directory for BAM workflow...")
        temp_dir_base_for_creation = custom_temp_dir_base if custom_temp_dir_base else output_dir_path
        main_temp_dir_path = get_temporary_directory(
            user_specified_dir=temp_dir_base_for_creation, prefix="mhb_from_bam_run_"
        )
        if not main_temp_dir_path or not main_temp_dir_path.exists():
            logger.error("Critical Error: Failed to create a main temporary directory for the BAM workflow. Aborting.")
            main_temp_dir_path = None
            return False
        logger.info(f"Main temporary directory for this run: '{main_temp_dir_path.resolve()}'")

        # --- Extract Methylation Calls ---
        logger.info("Extracting methylation calls from BAM files...")
        methylation_extraction_output_dir: Path
        if save_extracted_calls_to:
            methylation_extraction_output_dir = Path(save_extracted_calls_to).resolve()
            try:
                methylation_extraction_output_dir.mkdir(parents=True, exist_ok=True)
                logger.info(f"Extracted methylation call files will be saved to specified directory: '{methylation_extraction_output_dir}'")
                user_specified_calls_output_dir_for_cleanup_check = methylation_extraction_output_dir
            except OSError as e:
                logger.error(f"Directory Error: Could not create or access directory for saving extracted calls '{methylation_extraction_output_dir}': {e}. Aborting.")
                return False
        else:
            methylation_extraction_output_dir = main_temp_dir_path / "modkit_extracted_calls"
            methylation_extraction_output_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Extracted methylation call files will be saved to temporary directory: '{methylation_extraction_output_dir}'")

        actual_methylation_calls_files = mhb_api.extract_methylation_calls(
            bam_files=bam_file_paths,
            bed_file=regions_bed_file_path,
            output_dir=methylation_extraction_output_dir,
            fasta_ref=fasta_file_path,
            threads=num_threads
        )
        if not actual_methylation_calls_files:
            logger.error("Methylation extraction (modkit) step failed to produce any output files. Aborting workflow.")
            return False
        logger.info(f"Successfully generated {len(actual_methylation_calls_files)} methylation call files in '{methylation_extraction_output_dir}'.")

        # --- Call Core Processing ---
        if not _mhb_core_processing_and_saving(
            actual_methylation_calls_files=actual_methylation_calls_files,
            regions_bed_file_path=regions_bed_file_path,
            cpg_annotations_bed_path=cpg_annotations_bed_path,
            output_dir_path=output_dir_path,
            num_threads=num_threads,
            output_prefix=output_prefix,
            output_formats_validated=output_formats_validated,
            mld_pairing_threshold=mld_pairing_threshold,
            block_r2_soft_threshold=block_r2_soft_threshold,
            block_r2_hard_threshold=block_r2_hard_threshold,
            block_min_cpgs_in_block=block_min_cpgs_in_block,
            block_tolerance_gaps=block_tolerance_gaps
        ):
            logger.error("Core MHB processing and saving failed for BAM workflow.")
            return False

        logger.info("MHB Workflow (from BAMs) completed successfully.")
        return True

    except Exception as e:
        logger.error(f"MHB Workflow (from BAMs) failed with an unexpected error: {e}", exc_info=True)
        return False
    finally:
        # --- Cleanup Operations for BAM workflow ---
        if (user_specified_calls_output_dir_for_cleanup_check and
            user_specified_calls_output_dir_for_cleanup_check.exists() and
            main_temp_dir_path and main_temp_dir_path.exists()):
            logger.info(f"Checking for intermediate BAM/BAI files in '{user_specified_calls_output_dir_for_cleanup_check}' to move to temporary storage.")
            destination_for_moved_intermediates = main_temp_dir_path / "moved_intermediate_bams_from_saved_calls"
            try:
                destination_for_moved_intermediates.mkdir(parents=True, exist_ok=True)
            except OSError as e:
                logger.warning(f"Could not create directory '{destination_for_moved_intermediates}' for stashing intermediate BAMs. "
                               f"They may remain in '{user_specified_calls_output_dir_for_cleanup_check}'. Error: {e}")
                destination_for_moved_intermediates = None

            if destination_for_moved_intermediates and actual_methylation_calls_files:
                moved_count = 0
                for tsv_file_path in actual_methylation_calls_files:
                    if tsv_file_path.parent == user_specified_calls_output_dir_for_cleanup_check:
                        base_name_for_intermediates = tsv_file_path.stem
                        bam_to_move = user_specified_calls_output_dir_for_cleanup_check / f"{base_name_for_intermediates}.bam"
                        bai_to_move = user_specified_calls_output_dir_for_cleanup_check / f"{base_name_for_intermediates}.bam.bai"
                        for intermediate_fpath in [bam_to_move, bai_to_move]:
                            if intermediate_fpath.exists():
                                try:
                                    target_intermediate_path = destination_for_moved_intermediates / intermediate_fpath.name
                                    shutil.move(str(intermediate_fpath), str(target_intermediate_path))
                                    logger.debug(f"Moved '{intermediate_fpath.name}' to '{target_intermediate_path}'.")
                                    moved_count +=1
                                except Exception as e_move:
                                    logger.error(f"Failed to move intermediate file '{intermediate_fpath}' to "
                                                 f"'{destination_for_moved_intermediates}': {e_move}", exc_info=True)
                if moved_count > 0:
                    logger.info(f"Moved {moved_count} intermediate BAM/BAI files from '{user_specified_calls_output_dir_for_cleanup_check}' to '{destination_for_moved_intermediates}'.")
                else:
                    logger.info(f"No intermediate BAM/BAI files found in or moved from '{user_specified_calls_output_dir_for_cleanup_check}'.")

        if main_temp_dir_path and main_temp_dir_path.exists():
            if keep_temp_files:
                logger.info(f"Temporary files for BAM workflow are being kept at: '{main_temp_dir_path}'.")
            else:
                logger.info(f"Attempting to clean up temporary directory for BAM workflow: '{main_temp_dir_path}'.")
                try:
                    shutil.rmtree(main_temp_dir_path)
                    logger.info(f"Successfully removed temporary directory for BAM workflow: '{main_temp_dir_path}'.")
                except Exception as e:
                    logger.error(f"Cleanup Error: Failed to remove temporary directory '{main_temp_dir_path}' for BAM workflow: {e}", exc_info=True)
        elif main_temp_dir_path:
             logger.warning(f"Cleanup Warning: Main temporary directory '{main_temp_dir_path}' for BAM workflow was assigned but does not exist (or was already removed).")


def mhb_workflow_from_precomputed_calls(
    precomputed_calls_files: Union[Path, List[Path]],
    regions_bed_file_path: Path,
    cpg_annotations_bed_path: Path,
    output_dir_path: Path,
    custom_temp_dir_base: Optional[Path] = None,
    num_threads: int = 1,
    output_prefix: str = "mhb_analysis",
    output_formats: Optional[List[str]] = None,
    mld_pairing_threshold: int = 10,
    block_r2_soft_threshold: float = 0.5,
    block_r2_hard_threshold: float = 0.4,
    block_min_cpgs_in_block: int = 4,
    block_tolerance_gaps: int = 1,
    keep_temp_files: bool = False
) -> bool:
    """
    Orchestrates the MHB calling workflow starting from precomputed methylation call TSV files.
    """
    logger.info("Starting MHB Workflow (from precomputed calls)...")
    logger.info(f"Output directory set to: '{output_dir_path.resolve()}'")
    logger.info(f"Number of threads for parallel tasks: {num_threads}")

    output_formats_validated = _validate_and_prepare_output_formats(output_formats)
    if not output_formats_validated:
        return False

    # --- Validate Precomputed Calls Input ---
    actual_methylation_calls_files: List[Path] = []
    if precomputed_calls_files is None:
        logger.error("Input Error: 'precomputed_calls_files' were not provided. Workflow cannot proceed.")
        return False
    
    if isinstance(precomputed_calls_files, Path):
        input_path = Path(precomputed_calls_files).resolve()
        if input_path.is_file():
            if input_path.suffix.lower() == ".tsv":
                actual_methylation_calls_files = [input_path]
            else:
                logger.error(f"Input Error: Provided single precomputed file '{input_path}' is not a .tsv file.")
                return False
        elif input_path.is_dir():
            logger.info(f"Searching for '*.tsv' files in precomputed calls directory: '{input_path}'")
            actual_methylation_calls_files = sorted(list(input_path.glob("*.tsv")))
            if not actual_methylation_calls_files:
                logger.error(f"Input Error: Provided directory for precomputed files '{input_path}' contains no '*.tsv' files.")
                return False
            logger.info(f"Found {len(actual_methylation_calls_files)} .tsv files in directory '{input_path}'.")
        else:
            logger.error(f"Input Error: Provided path for 'precomputed_calls_files' ('{input_path}') is not a valid file or directory.")
            return False
    elif isinstance(precomputed_calls_files, list):
        temp_list = []
        all_files_valid = True
        for item in precomputed_calls_files:
            item_path = Path(item).resolve()
            if item_path.is_file() and item_path.suffix.lower() == ".tsv":
                temp_list.append(item_path)
            else:
                logger.error(f"Input Error: Item '{item_path}' in 'precomputed_calls_files' list is not a valid .tsv file or does not exist.")
                all_files_valid = False
                break
        if not all_files_valid:
            return False
        actual_methylation_calls_files = sorted(temp_list)
    else:
        logger.error("Input Error: 'precomputed_calls_files' must be a Path (to a file or directory) or a list of file paths.")
        return False
    
    if not actual_methylation_calls_files:
        logger.error("Input Error: No valid precomputed methylation call files to process.")
        return False
    logger.info(f"Using {len(actual_methylation_calls_files)} precomputed methylation call file(s).")

    if not _validate_common_inputs(regions_bed_file_path, cpg_annotations_bed_path, output_dir_path):
        return False

    main_temp_dir_path: Optional[Path] = None
    try:
        # --- Setup Temporary Directory ---
        logger.info("Setting up temporary directory for precomputed calls workflow...")
        temp_dir_base_for_creation = custom_temp_dir_base if custom_temp_dir_base else output_dir_path
        main_temp_dir_path = get_temporary_directory(
            user_specified_dir=temp_dir_base_for_creation, prefix="mhb_from_precomp_run_"
        )
        if not main_temp_dir_path or not main_temp_dir_path.exists():
            logger.error("Critical Error: Failed to create a main temporary directory for the precomputed calls workflow. Aborting.")
            main_temp_dir_path = None
            return False
        logger.info(f"Main temporary directory for this run: '{main_temp_dir_path.resolve()}'")
        
        # --- Call Core Processing ---
        if not _mhb_core_processing_and_saving(
            actual_methylation_calls_files=actual_methylation_calls_files,
            regions_bed_file_path=regions_bed_file_path,
            cpg_annotations_bed_path=cpg_annotations_bed_path,
            output_dir_path=output_dir_path,
            num_threads=num_threads,
            output_prefix=output_prefix,
            output_formats_validated=output_formats_validated,
            mld_pairing_threshold=mld_pairing_threshold,
            block_r2_soft_threshold=block_r2_soft_threshold,
            block_r2_hard_threshold=block_r2_hard_threshold,
            block_min_cpgs_in_block=block_min_cpgs_in_block,
            block_tolerance_gaps=block_tolerance_gaps
        ):
            logger.error("Core MHB processing and saving failed for precomputed calls workflow.")
            return False

        logger.info("MHB Workflow (from precomputed calls) completed successfully.")
        return True

    except Exception as e:
        logger.error(f"MHB Workflow (from precomputed calls) failed with an unexpected error: {e}", exc_info=True)
        return False
    finally:
        # --- Cleanup Operations for Precomputed Calls workflow ---
        if main_temp_dir_path and main_temp_dir_path.exists():
            if keep_temp_files:
                logger.info(f"Temporary files for precomputed calls workflow are being kept at: '{main_temp_dir_path}'.")
            else:
                logger.info(f"Attempting to clean up temporary directory for precomputed calls workflow: '{main_temp_dir_path}'.")
                try:
                    shutil.rmtree(main_temp_dir_path)
                    logger.info(f"Successfully removed temporary directory for precomputed calls workflow: '{main_temp_dir_path}'.")
                except Exception as e:
                    logger.error(f"Cleanup Error: Failed to remove temporary directory '{main_temp_dir_path}' for precomputed calls workflow: {e}", exc_info=True)
        elif main_temp_dir_path:
             logger.warning(f"Cleanup Warning: Main temporary directory '{main_temp_dir_path}' for precomputed calls workflow was assigned but does not exist (or was already removed).")