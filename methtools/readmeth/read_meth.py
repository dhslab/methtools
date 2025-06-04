import sys
import pandas as pd
import numpy as np
import pysam
import logging


def main(bam: list, haplotypes: list, out_file: str, regions_bed: str, sample: str):
    # Initialize logging configuration
    logger = logging.getLogger(__name__)

    if len(bam) == 1:
        logger.info(
            "Single BAM file provided. No haplotype consideration will be taken."
        )
        bams_dict = {sample: {"single": bam[0]}}
        haplotypes = ["single"]
    else:
        logger.info(
            "Two BAM files provided. Assuming ref and alt haplotypes unless overridden."
        )
        bams_dict = {sample: {haplotypes[0]: bam[0], haplotypes[1]: bam[1]}}

    # Add logging statements to provide information about the script's progress
    logger.info("Starting methylation data processing")
    logger.info(f"BAM files: {bam}")
    logger.info(f"Output file: {out_file}")
    logger.info(f"Regions BED file: {regions_bed}")
    logger.info(f"Sample name: {sample}")

    regions_df = pd.read_csv(regions_bed, sep="\t", header=None)

    # Log the number of regions being processed
    logger.info(f"Number of regions to process: {len(regions_df)}")

    all_reads_df_list = []

    columns = [
        "chrom",
        "start",
        "end",
        "read_name",
        "total_cpgs",
        "meth_cpgs",
        "unmeth_cpgs",
        "xm_tag",
        "sample",
        "read_pair",
        "supplementary",
    ]
    if len(bam) > 1:
        columns.insert(8, "haplotype")  # Add haplotype column only for ref/alt scenario

    all_reads_df = []
    for sample in bams_dict:
        logger.info(f"Processing sample: {sample}")
        for idx, row in regions_df.iterrows():
            chrom, start, end = row[0], row[1], row[2]
            # logger.info(f"Processing region: {chrom}:{start}-{end}")
            for haplotype in haplotypes:
                bam_file = bams_dict[sample][haplotype]
                bamfile = pysam.AlignmentFile(bam_file, "rb")
                for read in bamfile.fetch(chrom, start, end):
                    read_name = read.query_name
                    if read.has_tag("XM"):
                        # get read positions overlapping the region
                        try:
                            overlapping_pos = [
                                pair[0]
                                for pair in read.get_aligned_pairs()
                                if pair[1] >= start and pair[1] <= end
                            ]

                            # get min max values for overlapping positions
                            min_pos = np.min(overlapping_pos)
                            max_pos = np.max(overlapping_pos)
                            # subset xm_tag to only include positions within the region
                            xm_tag = read.get_tag("XM")
                            xm_tag = xm_tag[min_pos:max_pos]
                            if len(xm_tag) == 0:
                                continue
                            unmeth_cpgs = xm_tag.count("z")
                            meth_cpgs = xm_tag.count("Z")
                            total_cpgs = unmeth_cpgs + meth_cpgs
                            if read.is_read1:
                                read_nmber = "read1"
                            elif read.is_read2:
                                read_nmber = "read2"
                            elif read.is_read1 == False and read.is_read2 == False:
                                read_nmber = "unpaired"
                            # check if supplementary alignment
                            if read.is_supplementary:
                                read_sup = True
                            else:
                                read_sup = False
                            read_list = [
                                chrom,
                                start,
                                end,
                                read_name,
                                total_cpgs,
                                meth_cpgs,
                                unmeth_cpgs,
                                xm_tag,
                                sample,
                                read_nmber,
                                read_sup,
                            ]
                            if len(bam) > 1:
                                read_list.insert(8, haplotype)
                            all_reads_df.append(read_list)
                        except:
                            continue

    # Log completion of processing
    logger.info("Methylation data processing completed")

    # create a dataframe
    reads_df = pd.DataFrame(all_reads_df, columns=columns)
    reads_df.to_csv(out_file, index=False, sep="\t")


if __name__ == "__main__":
    raise RuntimeError(
        "This module is not meant to be run directly. Use the CLI entry point."
    )
