# methtools/mhb/api.py
"""
This module provides the core API for processing BAM files to extract
methylation data and subsequently call Methylation Haplotype Blocks (MHBs).
"""
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import floor
from pathlib import Path
from typing import List

import pandas as pd
import numpy as np
# from pybedtools import BedTool # This import was in block_code.py but not used by the integrated functions. Add if needed.

# Import shared utilities from other parts of the methtools package
from ..utils.system import run_command # Assuming system.py is in methtools/utils/

# Get a logger specific to this module for contextual logging
logger = logging.getLogger(__name__)

# --- Functions for extracting methylation calls from BAM files ---

def extract_methylation_calls(
    bam_files: List[Path],
    bed_file: Path,
    output_dir: Path,
    fasta_ref: Path,
    threads: int = 8,
    dry_run: bool = False,
) -> List[Path]:
    """
    Public API function to extract read-level methylation calls from BAM files.

    This function orchestrates the process of subsetting BAM files to
    specific regions and running modkit to generate methylation call files (TSVs),
    which can be used for downstream MHB analysis.

    Args:
        bam_files: List of paths to input BAM files.
        bed_file: Path to a BED file defining regions of interest.
        output_dir: Directory to save temporary and final output files.
        fasta_ref: Path to the reference FASTA file.
        threads: Total number of CPU threads to use for the entire run.
        dry_run: If True, log commands that would be run without executing them.

    Returns:
        A list of paths to the successfully generated TSV result files.
    """
    logger.info("Starting methylation extraction process...")

    # Here you can add validation for inputs, e.g., check if files exist
    for file_path_check in [bed_file, fasta_ref, *bam_files]:
        if not file_path_check.exists():
            msg = f"Input file not found: {file_path_check}"
            logger.error(msg)
            raise FileNotFoundError(msg)

    # Delegate the core processing to the internal parallelization function
    output_tsvs = _run_parallel_bam_processing(
        bam_files=bam_files,
        bed_file=bed_file,
        tmp_dir=output_dir,
        fasta_ref=fasta_ref,
        max_total_threads=threads,
        dry_run=dry_run,
    )

    if output_tsvs:
        logger.info("Successfully generated all result files for methylation extraction.")
    else:
        logger.warning(
            "Methylation extraction process completed, but no output files were generated."
        )

    return output_tsvs


def _run_parallel_bam_processing(
    bam_files: List[Path],
    bed_file: Path,
    tmp_dir: Path,
    fasta_ref: Path,
    max_total_threads: int,
    dry_run: bool,
) -> List[Path]:
    """Internal function to manage the parallel processing of multiple BAM files."""
    num_bams = len(bam_files)
    if num_bams == 0:
        logger.warning("No BAM files provided for methylation extraction processing.")
        return []

    max_workers = min(num_bams, max_total_threads)
    threads_per_job = max(1, floor(max_total_threads / max_workers))

    logger.info(f"Starting parallel BAM processing for {num_bams} BAM file(s).")
    logger.info(
        f"Distributing tasks across {max_workers} worker(s) with "
        f"{threads_per_job} thread(s) per worker (total threads: {max_total_threads})."
    )

    tmp_dir.mkdir(parents=True, exist_ok=True)

    successful_results = []
    failed_jobs = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_bam = {
            executor.submit(
                _process_bam_for_extraction, # Renamed for clarity
                bam, bed_file, tmp_dir, fasta_ref, threads_per_job, dry_run
            ): bam for bam in bam_files
        }

        for future in as_completed(future_to_bam):
            bam_path = future_to_bam[future]
            try:
                result_path = future.result()
                successful_results.append(result_path)
                logger.debug(f"Successfully processed {bam_path.name} for methylation extraction.")
            except Exception:
                logger.error(f"Job for {bam_path.name} (methylation extraction) failed during execution.")
                failed_jobs.append(bam_path.name)

    logger.info("-" * 50)
    logger.info("Parallel BAM processing for methylation extraction finished.")
    logger.info(f"  {len(successful_results)} jobs succeeded.")
    if failed_jobs:
        logger.warning(f"  {len(failed_jobs)} extraction jobs failed: {', '.join(failed_jobs)}")
    logger.info("-" * 50)

    return successful_results


def _process_bam_for_extraction( # Renamed for clarity
    bam_file: Path,
    bed_file: Path,
    tmp_dir: Path,
    fasta_ref: Path,
    threads_per_job: int,
    dry_run: bool,
) -> Path:
    """Internal worker function to process a single BAM file for methylation extraction."""
    prefix = bam_file.stem
    tmp_subset_bam = tmp_dir / f"{prefix}.subset.bam"

    logger.info(f"Worker started for BAM (methylation extraction): {bam_file.name}")

    samtools_view_cmd = [
        "samtools", "view", "-b",
        "-@", str(threads_per_job),
        "-L", str(bed_file),
        "-o", str(tmp_subset_bam),
        str(bam_file),
    ]
    samtools_index_cmd = ["samtools", "index", str(tmp_subset_bam)]

    try:
        run_command(samtools_view_cmd, dry_run=dry_run, workdir=tmp_dir)
        run_command(samtools_index_cmd, dry_run=dry_run, workdir=tmp_dir)
    except Exception as e:
        logger.error(
            f"Samtools processing failed for {bam_file.name} (methylation extraction).", exc_info=True
        )
        raise RuntimeError(f"Samtools processing failed for {bam_file.name}") from e

    tsv_output = tmp_subset_bam.with_suffix(".tsv")
    modkit_cmd_list = [
        "modkit", "extract", "calls",
        "--threads", str(threads_per_job),
        "--mapped-only",
        "--reference", str(fasta_ref),
        "--cpg",
        "--include-bed", str(bed_file),
        str(tmp_subset_bam),
        str(tsv_output),
    ]

    try:
        run_command(modkit_cmd_list, dry_run=dry_run, workdir=tmp_dir)
        logger.info(f"Worker finished for BAM (methylation extraction): {bam_file.name}, output: {tsv_output}")
    except Exception as e:
        logger.error(
            f"Modkit extraction failed for {tmp_subset_bam.name}.", exc_info=True
        )
        raise RuntimeError(f"Modkit extraction failed for {tmp_subset_bam.name}") from e

    return tsv_output

# --- Functions for MHB analysis ---

def generate_methylation_matrix(methylation_calls_files: List[Path]) -> pd.DataFrame:
    """Generates a consolidated methylation matrix from multiple TSV files.

    This function reads methylation call files (typically from 'modkit extract'),
    transforms the 'call_code' into a binary methylation status (0 for
    unmethylated, 1 for methylated), adds a sample identifier, and merges
    the data from all provided files into a single pandas DataFrame.

    Args:
        methylation_calls_files: A list of pathlib.Path objects, where each
            path points to a TSV file containing methylation calls.
            Expected columns in TSV: "read_id", "ref_position", "ref_strand",
            "call_code".

    Returns:
        A pandas DataFrame containing the merged methylation data.
        The DataFrame will have the following columns:
        ["read_id", "ref_position", "ref_strand", "sample_name", "methylation_call"]
        where "methylation_call" is 0 (unmethylated) or 1 (methylated).
        Returns an empty DataFrame with these columns if no input files
        are provided or if no data is processed.
    """
    processed_dataframes = []
    expected_input_cols = ["read_id", "ref_position", "ref_strand", "call_code"]
    final_output_cols = ["read_id", "ref_position", "ref_strand", "sample_name", "methylation_call"]

    if not methylation_calls_files:
        logger.warning("No methylation call files provided to generate_methylation_matrix. Returning empty DataFrame.")
        return pd.DataFrame(columns=final_output_cols)

    for file_path in methylation_calls_files:
        if not file_path.is_file():
            logger.warning(f"File not found in generate_methylation_matrix: {file_path}. Skipping.")
            continue

        sample_name = file_path.stem
        
        try:
            current_df = pd.read_csv(
                file_path,
                sep="\t",
                usecols=expected_input_cols,
            )

            if current_df.empty:
                logger.info(f"File {file_path} is empty or contains no relevant data for generate_methylation_matrix. Skipping.")
                continue

            methylation_map = {"-": 0, "m": 1, "h": 1} # Assuming 'h' (CHG) is also treated as methylated
            current_df["methylation_call"] = current_df["call_code"].map(methylation_map)
            
            current_df.dropna(subset=["methylation_call"], inplace=True)
            current_df["methylation_call"] = current_df["methylation_call"].astype(int)

            current_df["sample_name"] = sample_name
            
            processed_dataframes.append(current_df[final_output_cols])

        except pd.errors.EmptyDataError:
            logger.warning(f"File {file_path} is empty (EmptyDataError). Skipping in generate_methylation_matrix.")
            continue
        except KeyError as e:
            logger.error(f"Missing expected column in {file_path}: {e}. Skipping in generate_methylation_matrix.")
            continue
        except Exception as e:
            logger.error(f"Error processing file {file_path} in generate_methylation_matrix: {e}. Skipping.", exc_info=True)
            continue

    if not processed_dataframes:
        logger.info("No data processed from files in generate_methylation_matrix. Returning empty DataFrame.")
        return pd.DataFrame(columns=final_output_cols)

    methylation_matrix_df = pd.concat(processed_dataframes, ignore_index=True)
    
    return methylation_matrix_df[final_output_cols]


def collapse_cpg_positions(ref_position: int, ref_strand: str) -> int:
    """
    Collapses CpG positions to the plus strand C's coordinate.

    For CpG dinucleotides, the methylation call can be on the C of the plus
    strand or the C of the minus strand (which is G on the plus strand).
    This function standardizes these to the coordinate of the C on the plus strand.

    Args:
        ref_position: The 0-based reference position of the C in the CpG.
        ref_strand: The strand of the call ("+" or "-").

    Returns:
        The 0-based coordinate of the C on the plus strand for the CpG pair.
    """
    # If the call is on the plus strand, the position is already the C.
    # If on the minus strand, the call is on C (paired with G on plus).
    # The G on plus is at ref_position, so C on plus is ref_position - 1.
    return ref_position if ref_strand == "+" else ref_position - 1


def calculate_mld_from_pairs(methylation_df: pd.DataFrame, threshold: int = 10) -> pd.DataFrame:
    """
    Calculates methylation linkage disequilibrium (mLD) (r-squared)
    between adjacent CpG sites within reads.

    The input DataFrame should have columns representing CpG positions,
    and rows representing individual reads (or read segments). Values should
    be binary (0 for unmethylated, 1 for methylated).

    Args:
        methylation_df: Pandas DataFrame with CpG methylation data.
                        Rows are reads, columns are CpG positions.
        threshold: The minimum number of reads covering a pair of CpGs
                   required to calculate r-squared. If fewer, r-squared is NaN.

    Returns:
        A pandas DataFrame with a single row, where columns are "pos1_pos2"
        and values are the calculated r-squared for that pair.
    """
    if methylation_df.shape[1] < 2: # Need at least two CpG sites to form a pair
        logger.warning("Less than two CpG sites found in the input for mLD calculation. Returning empty DataFrame.")
        return pd.DataFrame()

    cpg_pairs = [(methylation_df.columns[i], methylation_df.columns[i + 1]) for i in range(methylation_df.shape[1] - 1)]
    r2_dict = {}

    for cpg1_pos, cpg2_pos in cpg_pairs:
        pair_df = methylation_df[[cpg1_pos, cpg2_pos]].dropna()
        pair_name = f"{cpg1_pos}_{cpg2_pos}"

        if pair_df.shape[0] < threshold:
            r2_dict[pair_name] = np.nan # Store as scalar, will become a list in DataFrame constructor if needed
            continue

        # Ensure data is numeric for correlation
        col1_values = pd.to_numeric(pair_df[cpg1_pos], errors="coerce")
        col2_values = pd.to_numeric(pair_df[cpg2_pos], errors="coerce")
        
        # Drop NaNs that might have been introduced by coerce
        valid_pair_df = pd.DataFrame({"c1": col1_values, "c2": col2_values}).dropna()

        if valid_pair_df.shape[0] < threshold:
            r2_dict[pair_name] = np.nan
            continue
            
        if valid_pair_df.shape[0] < 2: # np.corrcoef needs at least 2 observations
             r2_dict[pair_name] = np.nan
             continue

        # Check for zero variance in either column, which makes correlation undefined or 0
        if np.var(valid_pair_df["c1"]) == 0 or np.var(valid_pair_df["c2"]) == 0:
            correlation = 0.0 # Or np.nan, depending on desired behavior for no variance
        else:
            correlation_matrix = np.corrcoef(valid_pair_df["c1"], valid_pair_df["c2"])
            correlation = correlation_matrix[0, 1]
        
        r_squared = correlation**2
        r2_dict[pair_name] = r_squared
    
    # If r2_dict is empty, return an empty DataFrame
    if not r2_dict:
        return pd.DataFrame()

    return pd.DataFrame([r2_dict]) # Create a single-row DataFrame


def call_methylation_blocks(
    r_squared_df: pd.DataFrame,
    chromosome: str = "",
    soft_threshold: float = 0.5,
    hard_threshold: float = 0.4,
    min_cpgs_in_block: int = 4,
    tolerance_gaps: int = 1,
) -> pd.DataFrame:
    """
    Identifies methylation haplotype blocks (MHBs) from r-squared values.

    This function iterates through r-squared values between adjacent CpG sites
    and groups them into blocks based on defined thresholds and tolerance.

    Args:
        r_squared_df: A single-row DataFrame where columns are "pos1_pos2"
                      representing CpG pairs and values are their r-squared.
        chromosome: The chromosome name for the blocks (optional).
        soft_threshold: r-squared value above which a CpG pair strongly
                        supports being in a block.
        hard_threshold: r-squared value below which a CpG pair definitively
                        breaks a block. Pairs between hard and soft thresholds
                        can be tolerated up to `tolerance_gaps`.
        min_cpgs_in_block: Minimum number of CpG sites required to form a valid block.
        tolerance_gaps: Number of consecutive CpG pairs with r-squared between
                        `hard_threshold` and `soft_threshold` that are allowed
                        before a block is terminated.

    Returns:
        A pandas DataFrame where each row is an MHB, with columns:
        "chrm", "first_cpg_pos", "last_cpg_pos", "block_id",
        "all_cpg_pos" (comma-separated), "r2_scores" (comma-separated),
        "block_length", "num_cpgs".
    """
    if r_squared_df.empty or r_squared_df.shape[1] == 0:
        logger.info("Empty r-squared DataFrame provided to call_methylation_blocks. Returning empty DataFrame.")
        return pd.DataFrame(columns=["chrm", "first_cpg_pos", "last_cpg_pos", "block_id", "all_cpg_pos", "r2_scores", "block_length", "num_cpgs"])

    blocks = []
    current_block_cpgs = []
    current_block_r2_scores = []
    tolerated_gaps_count = 0

    # Iterate through the r-squared values (pairs of CpGs)
    # r_squared_df is expected to be a single row DataFrame
    for pair_name, r2_value_series in r_squared_df.items():
        r2_value = r2_value_series.iloc[0] # Get scalar r2 value

        if pd.isna(r2_value): # Skip pairs with NaN r-squared (e.g., due to insufficient data)
            # If we skip a pair, it breaks the current block
            if len(current_block_cpgs) >= min_cpgs_in_block:
                # Finalize and add the block before the gap
                first_cpg = min(current_block_cpgs)
                last_cpg = max(current_block_cpgs)
                blocks.append({
                    "chrm": chromosome,
                    "first_cpg_pos": first_cpg,
                    "last_cpg_pos": last_cpg,
                    "block_id": f"{chromosome}_{first_cpg}_{last_cpg}",
                    "all_cpg_pos": ",".join(map(str, sorted(list(set(current_block_cpgs))))),
                    "r2_scores": ",".join(map(lambda x: f"{x:.3f}", current_block_r2_scores)),
                    "block_length": last_cpg - first_cpg + 1,
                    "num_cpgs": len(set(current_block_cpgs)),
                })
            # Reset for the next potential block
            current_block_cpgs = []
            current_block_r2_scores = []
            tolerated_gaps_count = 0
            continue

        pos1, pos2 = map(int, pair_name.split("_"))

        # Decision logic for extending or terminating a block
        if r2_value >= soft_threshold:
            if not current_block_cpgs: # Start of a new block
                current_block_cpgs.extend([pos1, pos2])
            else:
                current_block_cpgs.append(pos2) # Extend current block
            current_block_r2_scores.append(r2_value)
            tolerated_gaps_count = 0 # Reset tolerance on a strong link
        elif r2_value >= hard_threshold: # Between hard and soft threshold
            if current_block_cpgs and tolerated_gaps_count < tolerance_gaps:
                current_block_cpgs.append(pos2)
                current_block_r2_scores.append(r2_value)
                tolerated_gaps_count += 1
            else: # Block broken by weak link or tolerance exceeded
                if len(current_block_cpgs) >= min_cpgs_in_block:
                    first_cpg = min(current_block_cpgs)
                    last_cpg = max(current_block_cpgs)
                    blocks.append({
                        "chrm": chromosome,
                        "first_cpg_pos": first_cpg,
                        "last_cpg_pos": last_cpg,
                        "block_id": f"{chromosome}_{first_cpg}_{last_cpg}",
                        "all_cpg_pos": ",".join(map(str, sorted(list(set(current_block_cpgs))))),
                        "r2_scores": ",".join(map(lambda x: f"{x:.3f}", current_block_r2_scores)),
                        "block_length": last_cpg - first_cpg + 1,
                        "num_cpgs": len(set(current_block_cpgs)),
                    })
                current_block_cpgs = []
                current_block_r2_scores = []
                tolerated_gaps_count = 0
        else: # r2_value < hard_threshold, block broken
            if len(current_block_cpgs) >= min_cpgs_in_block:
                first_cpg = min(current_block_cpgs)
                last_cpg = max(current_block_cpgs)
                blocks.append({
                    "chrm": chromosome,
                    "first_cpg_pos": first_cpg,
                    "last_cpg_pos": last_cpg,
                    "block_id": f"{chromosome}_{first_cpg}_{last_cpg}",
                    "all_cpg_pos": ",".join(map(str, sorted(list(set(current_block_cpgs))))),
                    "r2_scores": ",".join(map(lambda x: f"{x:.3f}", current_block_r2_scores)),
                    "block_length": last_cpg - first_cpg + 1,
                    "num_cpgs": len(set(current_block_cpgs)),
                })
            current_block_cpgs = []
            current_block_r2_scores = []
            tolerated_gaps_count = 0

    # Add the last block if it's valid
    if len(current_block_cpgs) >= min_cpgs_in_block:
        first_cpg = min(current_block_cpgs)
        last_cpg = max(current_block_cpgs)
        blocks.append({
            "chrm": chromosome,
            "first_cpg_pos": first_cpg,
            "last_cpg_pos": last_cpg,
            "block_id": f"{chromosome}_{first_cpg}_{last_cpg}",
            "all_cpg_pos": ",".join(map(str, sorted(list(set(current_block_cpgs))))),
            "r2_scores": ",".join(map(lambda x: f"{x:.3f}", current_block_r2_scores)),
            "block_length": last_cpg - first_cpg + 1,
            "num_cpgs": len(set(current_block_cpgs)),
        })

    if not blocks:
        return pd.DataFrame(columns=["chrm", "first_cpg_pos", "last_cpg_pos", "block_id", "all_cpg_pos", "r2_scores", "block_length", "num_cpgs"])

    return pd.DataFrame(blocks)

