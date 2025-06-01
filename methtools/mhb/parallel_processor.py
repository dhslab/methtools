import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Any, Optional
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import the sibling api module directly.
# This works when parallel_processor.py is treated as part of a package.
from . import api as mhb_api

# Get a logger for this module.
# This logger will inherit its configuration (handlers, level, formatter)
# from the root logger, which should be configured by the main application
# (e.g., your CLI script calling setup_logging()).
logger = logging.getLogger(__name__)


def find_mhb_in_region(
    calls_df: pd.DataFrame,
    region_chrom: str,
    region_start: int,
    region_end: int,
    all_regions_cpg_positions_df: pd.DataFrame,
    mld_pairing_threshold: int,
    block_r2_soft_threshold: float,
    block_r2_hard_threshold: float,
    block_min_cpgs_in_block: int,
    block_tolerance_gaps: int
) -> pd.DataFrame:
    """
    Processes a single genomic region to find Methylation Haplotype Blocks (MHBs).

    Args:
        calls_df: DataFrame containing methylation calls.
        region_chrom: Chromosome of the region.
        region_start: Start coordinate of the region.
        region_end: End coordinate of the region.
        all_regions_cpg_positions_df: DataFrame with CpG positions for all regions.
        mld_pairing_threshold: Min reads for mLD calculation.
        block_r2_soft_threshold: Soft r-squared threshold for block extension.
        block_r2_hard_threshold: Hard r-squared threshold for block termination.
        block_min_cpgs_in_block: Min CpGs to form a block.
        block_tolerance_gaps: Tolerated r-squared gaps.

    Returns:
        A pandas DataFrame of MHBs for the region, or an empty DataFrame.
    """
    empty_block_df_columns = [
        "chrm", "first_cpg_pos", "last_cpg_pos", "block_id",
        "all_cpg_pos", "r2_scores", "block_length", "num_cpgs",
        "region_chrom", "region_start", "region_end", "region_id"
    ]
    # logger.debug(f"Processing region: {region_chrom}:{region_start}-{region_end}")

    region_specific_cpgs_df = all_regions_cpg_positions_df[
        (all_regions_cpg_positions_df["chrom"] == region_chrom) &
        (all_regions_cpg_positions_df["cpg_pos"] >= region_start) &
        (all_regions_cpg_positions_df["cpg_pos"] <= region_end)
    ][["cpg_pos"]].drop_duplicates()
    
    region_cpg_positions = sorted(region_specific_cpgs_df["cpg_pos"].tolist())

    if not region_cpg_positions:
        logger.info(f"No CpG positions found in region {region_chrom}:{region_start}-{region_end}. Skipping.")
        return pd.DataFrame(columns=empty_block_df_columns)

    region_methylation_data = calls_df[
        calls_df["ref_position_collapsed"].isin(region_cpg_positions)
    ]

    if region_methylation_data.empty:
        logger.info(f"No methylation calls found for CpGs in region {region_chrom}:{region_start}-{region_end}. Skipping.")
        return pd.DataFrame(columns=empty_block_df_columns)

    region_methylation_data = region_methylation_data[
        ["read_id", "ref_position_collapsed", "sample_name", "methylation_call"]
    ].rename(columns={"ref_position_collapsed": "cpg_pos"})

    region_pivot_mtx = region_methylation_data.pivot_table(
        index=["read_id", "sample_name"],
        columns="cpg_pos",
        values="methylation_call"
    )

    all_region_cpg_cols = pd.Index(region_cpg_positions)
    missing_cpg_positions = all_region_cpg_cols.difference(region_pivot_mtx.columns)
    
    if not missing_cpg_positions.empty:
        missing_df = pd.DataFrame(
            index=region_pivot_mtx.index,
            columns=missing_cpg_positions,
            data=np.nan
        )
        region_pivot_mtx = pd.concat([region_pivot_mtx, missing_df], axis=1)
    
    region_pivot_mtx = region_pivot_mtx.reindex(columns=all_region_cpg_cols)

    # Use the module-level imported mhb_api
    r_squared_df = mhb_api.calculate_mld_from_pairs(
        region_pivot_mtx, threshold=mld_pairing_threshold
    )

    if r_squared_df.empty:
        logger.info(f"No r-squared values calculated for region {region_chrom}:{region_start}-{region_end}. Skipping.")
        return pd.DataFrame(columns=empty_block_df_columns)
        
    block_df = mhb_api.call_methylation_blocks(
        r_squared_df,
        chromosome=region_chrom,
        soft_threshold=block_r2_soft_threshold,
        hard_threshold=block_r2_hard_threshold,
        min_cpgs_in_block=block_min_cpgs_in_block,
        tolerance_gaps=block_tolerance_gaps
    )

    if block_df.empty:
        return pd.DataFrame(columns=empty_block_df_columns)

    block_df["region_chrom"] = region_chrom
    block_df["region_start"] = region_start
    block_df["region_end"] = region_end
    block_df["region_id"] = f"{region_chrom}_{region_start}_{region_end}"
    
    return block_df


def find_mhb_for_all_regions_parallel(
    bed_df: pd.DataFrame,
    calls_df: pd.DataFrame,
    all_regions_cpg_positions_df: pd.DataFrame,
    mld_pairing_threshold: int,
    block_r2_soft_threshold: float,
    block_r2_hard_threshold: float,
    block_min_cpgs_in_block: int,
    block_tolerance_gaps: int,
    num_threads: int = 1
) -> pd.DataFrame:
    """
    Orchestrates finding MHBs across multiple regions in parallel.

    Args:
        bed_df: DataFrame defining genomic regions. Expected columns are
                typically the first three: chromosome, start, end.
        calls_df: DataFrame with all methylation calls.
        all_regions_cpg_positions_df: DataFrame with all CpG positions.
        mld_pairing_threshold: Threshold for mLD calculation.
        block_r2_soft_threshold: Soft r2 threshold for blocks.
        block_r2_hard_threshold: Hard r2 threshold for blocks.
        block_min_cpgs_in_block: Minimum CpGs per block.
        block_tolerance_gaps: Tolerance for gaps in blocks.
        num_threads: Number of threads to use for parallel processing.

    Returns:
        A pandas DataFrame concatenating all MHBs found across all regions.
    """
    logger.info(f"Starting MHB detection for {bed_df.shape[0]} regions using {num_threads} thread(s).")
    
    regions_blocks_df_list = []
    
    if num_threads < 1:
        logger.warning("Number of threads specified is less than 1. Defaulting to 1 thread.")
        num_threads = 1

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        future_to_region_info = {} # Store more info for better error reporting if needed
        for index, row in bed_df.iterrows():
            # Assuming bed_df columns are 0: chrom, 1: start, 2: end
            # If they have names, use row["chrom_column_name"], etc.
            chrom = row.iloc[0] 
            start = row.iloc[1]
            end = row.iloc[2]

            future = executor.submit(
                find_mhb_in_region, 
                calls_df.copy(), # Pass a copy to ensure thread safety if modifications occur
                chrom,
                start,
                end,
                all_regions_cpg_positions_df.copy(), # Same as above
                mld_pairing_threshold,
                block_r2_soft_threshold,
                block_r2_hard_threshold,
                block_min_cpgs_in_block,
                block_tolerance_gaps
            )
            future_to_region_info[future] = {"id": f"{chrom}:{start}-{end}", "index": index}

        processed_regions_count = 0
        for future in as_completed(future_to_region_info):
            region_info = future_to_region_info[future]
            region_id_str = region_info["id"]
            try:
                single_region_blocks_df = future.result()
                if not single_region_blocks_df.empty:
                    regions_blocks_df_list.append(single_region_blocks_df)
                processed_regions_count +=1
                # Log progress periodically
                if processed_regions_count % 10 == 0 or processed_regions_count == bed_df.shape[0]:
                     logger.info(f"Processed {processed_regions_count}/{bed_df.shape[0]} regions...")
            except Exception as exc:
                logger.error(f"Region {region_id_str} (index {region_info['index']}) generated an exception: {exc}", exc_info=True)

    final_empty_block_df_columns = [
        "chrm", "first_cpg_pos", "last_cpg_pos", "block_id",
        "all_cpg_pos", "r2_scores", "block_length", "num_cpgs",
        "region_chrom", "region_start", "region_end", "region_id"
    ]

    if regions_blocks_df_list:
        all_blocks_df = pd.concat(regions_blocks_df_list, ignore_index=True)
        logger.info(f"MHB detection complete. Found {len(all_blocks_df)} blocks in {processed_regions_count} successfully processed regions (out of {bed_df.shape[0]} total).")
    else:
        logger.warning("No blocks found in any region after parallel processing.")
        all_blocks_df = pd.DataFrame(columns=final_empty_block_df_columns)
        
    return all_blocks_df
