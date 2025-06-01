# methtools/mhb/workflows.py
"""
This module defines high-level workflows for Methylation Haplotype Block (MHB)
analysis, orchestrating calls to the core API and processing steps.
It is intended to be called by the CLI command definitions.
"""
import logging
from pathlib import Path
from typing import List, Optional, Union # Added Union
import shutil

import pandas as pd
import numpy as np
from pybedtools import BedTool

from ..utils.system import get_temporary_directory
from . import api as mhb_api
from .parallel_processor import find_mhb_for_all_regions_parallel

logger = logging.getLogger(__name__)


def mhb_calling_workflow(
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
    skip_methylation_extraction: bool = False,
    precomputed_calls_files: Optional[Union[Path, List[Path]]] = None, # Updated type hint
    keep_temp_files: bool = False
) -> bool:
    """
    Orchestrates the complete Methylation Haplotype Block (MHB) calling workflow.

    Args:
        fasta_file_path: Path to the reference genome FASTA file.
        bam_file_paths: List of paths to input BAM/CRAM files.
        regions_bed_file_path: Path to a BED file defining regions for MHB analysis.
        cpg_annotations_bed_path: Path to a BED file with all CpG positions.
        output_dir_path: Path to the directory where final results will be saved.
        custom_temp_dir_base: Optional path to a base directory for temporary files.
        num_threads: Number of threads for parallel processing steps.
        output_prefix: Prefix for output result files.
        output_formats: List of desired output formats (e.g., ["tsv", "bed"]).
                        Defaults to ["tsv", "bed"] if None.
        mld_pairing_threshold: Minimum read coverage for mLD calculation.
        block_r2_soft_threshold: Soft r-squared threshold for block definition.
        block_r2_hard_threshold: Hard r-squared threshold for block definition.
        block_min_cpgs_in_block: Minimum number of CpGs to constitute a block.
        block_tolerance_gaps: Allowed number of r-squared values between soft and
                              hard thresholds before terminating a block.
        skip_methylation_extraction: If True, assumes methylation call TSV files
                                     are already generated and provided via
                                     `precomputed_calls_files`.
        precomputed_calls_files: Optional. Can be a single Path to a .tsv file,
                                 a Path to a directory containing .tsv files,
                                 or a list of Paths to .tsv files.
                                 Used if `skip_methylation_extraction` is True.
        keep_temp_files: If True, temporary files will not be deleted after the run.

    Returns:
        True if the workflow completed successfully, False otherwise.
    """
    logger.info("Starting MHB Calling Workflow...")
    logger.info(f"Output directory set to: '{output_dir_path.resolve()}'")
    logger.info(f"Number of threads for parallel tasks: {num_threads}")

    if output_formats is None:
        output_formats = ["tsv", "bed"]

    # --- 1. Validate Input Parameters ---
    logger.info("Step 1: Validating input parameters...")
    required_files_for_run = {
        "FASTA reference": fasta_file_path,
        "Regions BED file": regions_bed_file_path,
        "CpG annotations BED file": cpg_annotations_bed_path,
    }
    for name, file_path in required_files_for_run.items():
        if not file_path.is_file():
            logger.error(f"Input Error: {name} not found at '{file_path}'. Workflow cannot proceed.")
            return False
    
    actual_methylation_calls_files: List[Path] = [] # Initialize

    if skip_methylation_extraction:
        if precomputed_calls_files is None:
            logger.error("Input Error: 'skip_methylation_extraction' is True, but no 'precomputed_calls_files' were provided.")
            return False
        
        if isinstance(precomputed_calls_files, Path):
            input_path = Path(precomputed_calls_files).resolve() # Ensure it's a Path and absolute
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
                item_path = Path(item).resolve() # Ensure it's a Path and absolute
                if item_path.is_file() and item_path.suffix.lower() == ".tsv":
                    temp_list.append(item_path)
                else:
                    logger.error(f"Input Error: Item '{item_path}' in 'precomputed_calls_files' list is not a valid .tsv file or does not exist.")
                    all_files_valid = False
                    break
            if not all_files_valid:
                return False
            actual_methylation_calls_files = sorted(temp_list) # Sort for consistent order
        else:
            logger.error("Input Error: 'precomputed_calls_files' must be a Path (to a file or directory) or a list of file paths.")
            return False
        
        if not actual_methylation_calls_files:
            logger.error("Input Error: No valid precomputed methylation call files to process after evaluating 'precomputed_calls_files'.")
            return False
        logger.info(f"Using {len(actual_methylation_calls_files)} precomputed methylation call file(s).")

    else: # We need to extract calls, so BAM files are required
        if not bam_file_paths:
            logger.error("Input Error: No BAM files provided and 'skip_methylation_extraction' is False. Workflow cannot proceed.")
            return False
        for bam_path in bam_file_paths:
            if not bam_path.is_file():
                logger.error(f"Input Error: Input BAM file not found: '{bam_path}'. Workflow cannot proceed.")
                return False
    
    try:
        output_dir_path.mkdir(parents=True, exist_ok=True)
        logger.info(f"Ensured output directory exists: '{output_dir_path.resolve()}'")
    except OSError as e:
        logger.error(f"Directory Error: Could not create or access output directory '{output_dir_path}': {e}. Workflow cannot proceed.")
        return False

    valid_formats = []
    for fmt in output_formats:
        fmt_lower = fmt.strip().lower()
        if fmt_lower in ["tsv", "bed"]:
            valid_formats.append(fmt_lower)
        else:
            logger.warning(f"Configuration Warning: Unsupported output format '{fmt}' specified. It will be ignored.")
    if not valid_formats:
        logger.error("Configuration Error: No valid output formats specified (supported: tsv, bed). Aborting workflow.")
        return False
    output_formats_validated = valid_formats
    logger.info(f"Validated output formats to be generated: {output_formats_validated}")

    # --- 2. Setup Temporary Directory ---
    logger.info("Step 2: Setting up temporary directory...")
    temp_dir_base_for_creation = custom_temp_dir_base if custom_temp_dir_base else output_dir_path / "tmp_methtools_workflow"
    
    main_temp_dir_path = get_temporary_directory(
        user_specified_dir=temp_dir_base_for_creation,
        prefix="mhb_run_"
    )
    if not main_temp_dir_path or not main_temp_dir_path.exists():
        logger.error("Critical Error: Failed to create a main temporary directory for the workflow. Aborting.")
        return False
    logger.info(f"Main temporary directory for this run: '{main_temp_dir_path.resolve()}'")

    try:
        # --- 3. Extract Methylation Calls (if not skipped earlier) ---
        if not skip_methylation_extraction: # This re-checks because actual_methylation_calls_files is set above if skipping
            logger.info("Step 3: Extracting methylation calls from BAM files...")
            modkit_output_dir = main_temp_dir_path / "modkit_extracted_calls"
            modkit_output_dir.mkdir(parents=True, exist_ok=True)

            actual_methylation_calls_files = mhb_api.extract_methylation_calls(
                bam_files=bam_file_paths,
                bed_file=regions_bed_file_path,
                output_dir=modkit_output_dir,
                fasta_ref=fasta_file_path,
                threads=num_threads 
            )
            if not actual_methylation_calls_files:
                logger.error("Methylation extraction (modkit) step failed to produce any output files. Aborting workflow.")
                return False
            logger.info(f"Successfully generated {len(actual_methylation_calls_files)} methylation call files in '{modkit_output_dir}'.")
        # If we skipped extraction, actual_methylation_calls_files is already populated and validated.

        # --- 4. Intersect Analysis Regions with All CpG Positions ---
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
            logger.error("No CpG sites found within the specified analysis regions after intersection with annotations. Aborting workflow.")
            return False
        logger.info(f"Identified {len(all_workflow_cpg_positions_df)} unique CpG sites within the analysis regions.")

        # --- 5. Generate Consolidated Methylation Calls Matrix ---
        logger.info("Step 5: Generating consolidated methylation matrix from call files...")
        calls_df = mhb_api.generate_methylation_matrix(actual_methylation_calls_files)
        if calls_df.empty:
            logger.error("Consolidated methylation matrix is empty after processing call files. Aborting workflow.")
            return False
        logger.info(f"Consolidated methylation matrix created with {calls_df.shape[0]} methylation calls.")

        # --- 6. Collapse CpG Positions ---
        logger.info("Step 6: Collapsing CpG positions to a common reference strand...")
        if "ref_position" not in calls_df.columns or "ref_strand" not in calls_df.columns:
            logger.error("Required columns 'ref_position' or 'ref_strand' not found in methylation matrix for collapsing. Aborting workflow.")
            return False
            
        calls_df["ref_position_collapsed"] = calls_df.apply(
            lambda row: mhb_api.collapse_cpg_positions(row["ref_position"], row["ref_strand"]),
            axis=1
        )
        logger.info("CpG positions successfully collapsed.")

        # --- 7. Call MHBs in Parallel ---
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

        # --- 8. Save Results ---
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
                            "chrm", "first_cpg_pos", "last_cpg_pos", "region_start", "region_end"
                        ]].copy()
                        
                        bed_output_df.to_csv(output_file_path, sep="\t", index=False, header=False)
                        logger.info(f"Saved MHB blocks (BED) to: '{output_file_path}'")
                except Exception as e:
                    logger.error(f"Failed to save results in '{fmt}' format to '{output_file_path}': {e}", exc_info=True)
        else:
            logger.info("No MHB blocks to save as none were identified.")
        
        logger.info("MHB Calling Workflow completed successfully.")
        return True

    except FileNotFoundError as e:
        logger.error(f"Workflow Error: Missing file - {e}", exc_info=True)
        return False
    except pd.errors.EmptyDataError as e:
        logger.error(f"Workflow Error: Empty input file encountered - {e}", exc_info=True)
        return False
    except KeyError as e:
        logger.error(f"Workflow Error: Missing expected column in data - {e}", exc_info=True)
        return False
    except Exception as e:
        logger.error(f"Workflow Error: An unexpected error occurred - {e}", exc_info=True)
        return False
    finally:
        # --- 9. Cleanup Temporary Files ---
        if main_temp_dir_path and main_temp_dir_path.exists():
            if keep_temp_files:
                logger.info(f"Temporary files and directories are being kept at: '{main_temp_dir_path}' as per configuration.")
            else:
                logger.info(f"Attempting to clean up temporary directory: '{main_temp_dir_path}'")
                try:
                    shutil.rmtree(main_temp_dir_path)
                    logger.info(f"Successfully removed temporary directory: '{main_temp_dir_path}'")
                except Exception as e:
                    logger.error(f"Cleanup Error: Failed to remove temporary directory '{main_temp_dir_path}': {e}", exc_info=True)
        elif main_temp_dir_path: 
             logger.warning(f"Cleanup Warning: Main temporary directory '{main_temp_dir_path}' was assigned but does not exist (or was already removed).")

