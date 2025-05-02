import os
import pandas as pd
from pybedtools import BedTool
import pysam
import logging
import sys

pysam.set_verbosity(0)
import warnings

warnings.simplefilter("ignore", UserWarning)


def setup_logging():
    """
    Configures the logging to write INFO and above messages to stderr.
    """
    # Create a logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)  # Set the minimum logging level

    # Define a logging format
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # Create a StreamHandler for stderr
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.INFO)  # Set handler level
    stderr_handler.setFormatter(formatter)  # Set formatter

    # Add the handler to the logger
    logger.addHandler(stderr_handler)


# Initialize logging configuration
setup_logging()


def get_read_strand(read):
    """
    Determine the originating strand for a bisulfite sequencing read.

    For directional bisulfite libraries, the following conventions are applied:
      - Paired-end reads:
          * Read1:
              - Not reversed (read.is_reverse == False): top strand.
              - Reversed (read.is_reverse == True): bottom strand.
          * Read2:
              - Reversed (read.is_reverse == True): top strand.
              - Not reversed (read.is_reverse == False): bottom strand.
      - Single-end reads:
          * Not reversed: top strand.
          * Reversed: bottom strand.

    If the read's flags are ambiguous (e.g., neither read1 nor read2 is set in a paired read),
    the function returns 'unknown'.

    Parameters:
        read (pysam.AlignedSegment): A read from a BAM/SAM file.

    Returns:
        str: "top" if the read is from the top strand,
             "bottom" if from the bottom strand,
             "unknown" if it cannot be determined.
    """
    # Check if the read is part of a paired-end experiment.
    if read.is_paired:
        # In a properly paired read, either is_read1 or is_read2 should be True.
        if read.is_read1:
            # For read1 in a directional library:
            #   - If the read is not reversed, it originates from the top strand.
            #   - If the read is reversed, it originates from the bottom strand.
            return "top" if not read.is_reverse else "bottom"
        elif read.is_read2:
            # For read2 in a directional library:
            #   - If the read is reversed, it originates from the top strand.
            #   - If the read is not reversed, it originates from the bottom strand.
            return "top" if read.is_reverse else "bottom"
        else:
            # In case neither is set (this is unusual), we cannot decide.
            return "unknown"
    else:
        # For single-end reads, assume:
        #   - Not reversed implies top strand.
        #   - Reversed implies bottom strand.
        return "top" if not read.is_reverse else "bottom"


def classify_read_allele(read_base, read_strand, snp_ref, snp_alt):
    """
    Classify the allele supported by a read at a SNP position based on bisulfite
    conversion logic similar to the C++ implementation.

    Parameters:
        read_base (str): The base observed in the read at the SNP position.
        read_strand (str): "top" or "bottom", indicating the strand of the read.
        snp_ref (str): The reference allele (e.g., 'C').
        snp_alt (str): The alternate allele (e.g., 'T').

    Returns:
        str: "ref" if the read supports the reference allele,
             "alt" if it supports the alternate allele,
             "unknown" if the read is ambiguous or uninformative.
    """
    # Ensure that all bases are uppercase.
    read_base = read_base.upper()
    snp_ref = snp_ref.upper()
    snp_alt = snp_alt.upper()

    # === Step 1: Early Discard for Ambiguous Cases ===
    # For C/T SNPs: if the read is from the top strand, the observed T might be due to bisulfite conversion.
    if (
        (snp_ref == "C" and snp_alt == "T") or (snp_ref == "T" and snp_alt == "C")
    ) and read_strand == "top":
        return "unknown"

    # For G/A SNPs: if the read is from the bottom strand, the observed A might be due to bisulfite conversion.
    if (
        (snp_ref == "G" and snp_alt == "A") or (snp_ref == "A" and snp_alt == "G")
    ) and read_strand == "bottom":
        return "unknown"

    # === Step 2: Define Allowed Bases for Each Allele ===
    # For the reference allele:
    if snp_ref == "C" and snp_alt != "T" and read_strand == "top":
        # On the top strand, if the allele is C (and not a C/T pair), allow both C and T.
        allowed_ref = {"C", "T"}
    elif snp_ref == "G" and snp_alt != "A" and read_strand == "bottom":
        # On the bottom strand, if the allele is G (and not a G/A pair), allow both G and A.
        allowed_ref = {"G", "A"}
    else:
        allowed_ref = {snp_ref}

    # For the alternate allele:
    if snp_alt == "C" and snp_ref != "T" and read_strand == "top":
        allowed_alt = {"C", "T"}
    elif snp_alt == "G" and snp_ref != "A" and read_strand == "bottom":
        allowed_alt = {"G", "A"}
    else:
        allowed_alt = {snp_alt}

    # === Step 3: Assign the Read Allele Based on the Observed Base ===
    if read_base in allowed_ref:
        return "ref"
    elif read_base in allowed_alt:
        return "alt"
    else:
        return "unknown"


def main(args):
    bam_path = args.bam_path
    ref_fasta = args.ref_fasta
    vcf_path = args.vcf_path
    regions_path = args.regions_path
    results_path = args.results_path
    bp_distance = args.bp_distance
    out_prefix = args.out_prefix

    # Ensure the results directory exists
    os.makedirs(results_path, exist_ok=True)

    # Prepare SNPs overlapping Regions of Interest
    logging.info("Preparing SNPs overlapping Regions of Interest")
    regions_df = pd.read_csv(regions_path, sep="\t", header=None)[[0, 1, 2]]

    regions_df_bed = BedTool.from_dataframe(regions_df).sort()
    regions_df_merged = regions_df_bed.merge(d=bp_distance).to_dataframe()
    regions_df_merged["block_id"] = (
        regions_df_merged["chrom"]
        + "_"
        + regions_df_merged["start"].astype(str)
        + "_"
        + regions_df_merged["end"].astype(str)
    )
    # overlap with original regions
    mergeds_with_nonmerged_regions = (
        BedTool.from_dataframe(regions_df_merged)
        .intersect(regions_path, wa=True, wb=True)
        .to_dataframe(
            disable_auto_names=True,
            header=None,
            names=[
                "merged_chrom",
                "merged_start",
                "merged_end",
                "merged_block_id",
                "regions_chrom",
                "regions_start",
                "regions_end",
                "regions_block_id",
            ],
        )
    )

    # intersect with vcf
    vcf_with_regions = (
        BedTool(vcf_path)
        .intersect(BedTool.from_dataframe(regions_df_merged), wa=True, wb=True)
        .to_dataframe(
            disable_auto_names=True,
            header=None,
        )
    )
    vcf_with_regions = vcf_with_regions[[0, 1, 3, 4, 8, 9, 10, 11]]
    vcf_with_regions.columns = [
        "snp_chrom",
        "snp_pos",
        "ref",
        "alt",
        "merged_chrom",
        "merged_start",
        "merged_end",
        "merged_block_id",
    ]
    # merge with mergeds_with_nonmerged_regions
    vcf_with_regions = vcf_with_regions.merge(
        mergeds_with_nonmerged_regions,
        left_on=["merged_chrom", "merged_start", "merged_end", "merged_block_id"],
        right_on=["merged_chrom", "merged_start", "merged_end", "merged_block_id"],
        how="left",
    )
    # add regions_id
    vcf_with_regions["regions_id"] = (
        vcf_with_regions["regions_chrom"]
        + "_"
        + vcf_with_regions["regions_start"].astype(str)
        + "_"
        + vcf_with_regions["regions_end"].astype(str)
    )

    # loop over each block and then snps overlapping with it
    # open bam file
    # set default base quality and mapping quality
    base_qual = 0
    map_qual = 0
    read_coverage_threshold = 5
    vaf_threshold = 0.15
    bamfile = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_fasta)

    # classifcation counter
    classification_counter = {
        "ref": 0,
        "alt": 0,
        "unknown": 0,
    }

    # logging.info(f"Classifying reads at SNPs in {vcf_with_regions.shape[0]} regions blocks")
    # log threshold
    logging.info(
        (
            f"Classifying reads at {vcf_with_regions.shape[0]} Common SNPs overlapping {len(regions_df_merged)} regions blocks\n"
            "Filtering thresholds:\n"
            "\n"
            "#######################\n"
            f"Base quality threshold: {base_qual}\n"
            f"Mapping quality threshold:\t{map_qual}\n"
            f"Read coverage threshold:\t{read_coverage_threshold}\n"
            f"VAF threshold:\t>= {0.5-vaf_threshold} and <= {0.5 + vaf_threshold}\n"
            "#######################\n"
        )
    )
    # the idea of merged blocks is to pool reads from all regions that are within N bp from each other in a list and to avoid using the same read twice
    all_counter = 0
    het_counter = 0
    ref_merged_block_reads_dict = {}
    het_snps_dict = {}
    for merged_block in vcf_with_regions["merged_block_id"].unique():
        ref_merged_block_reads_dict[merged_block] = {}
        ref_merged_block_reads_list = []
        alt_merged_block_reads_list = []
        merged_block_df = vcf_with_regions[
            vcf_with_regions["merged_block_id"] == merged_block
        ]
        for regions_idx, regions_row in merged_block_df.iterrows():
            all_counter += 1
            regions_id = regions_row["regions_id"]
            regions_chrom = regions_row["regions_chrom"]
            regions_start = regions_row["regions_start"]
            regions_end = regions_row["regions_end"]
            # iterate over snps in block
            regions_snps = merged_block_df[merged_block_df["regions_id"] == regions_id]
            for snp_idx, snp_row in regions_snps.iterrows():
                snp_chrom = snp_row["snp_chrom"]
                snp_pos = snp_row["snp_pos"]
                snp_ref = snp_row["ref"]
                snp_alt = snp_row["alt"]
                snp_id = f"{snp_chrom}_{snp_pos}_{snp_ref}_{snp_alt}"
                snp_ref_merged_block_reads_list = []
                snp_alt_merged_block_reads_list = []
                # iterate over reads pileup at snp position
                for pileup in bamfile.pileup(snp_chrom, snp_pos - 1, snp_pos):
                    if pileup.pos == snp_pos - 1:
                        for pileupread in pileup.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                # get read id, strand, mate number
                                read_id = pileupread.alignment.qname
                                read_pair = (
                                    "read1"
                                    if pileupread.alignment.is_read1
                                    else (
                                        "read2"
                                        if pileupread.alignment.is_read2
                                        else "unpaired"
                                    )
                                )
                                # Get the base from the read at the SNP position
                                read_base = pileupread.alignment.query_sequence[
                                    pileupread.query_position
                                ]
                                # get the base quality from the read at the SNP position
                                read_base_qual = pileupread.alignment.query_qualities[
                                    pileupread.query_position
                                ]
                                if read_base_qual < base_qual:
                                    print(
                                        f"base qual too low: {read_base_qual}, skipping"
                                    )
                                    continue
                                # get the mapping quality from the read
                                read_map_qual = pileupread.alignment.mapping_quality
                                if read_map_qual < map_qual:
                                    print(
                                        f"map qual too low: {read_map_qual}, skipping"
                                    )
                                    continue
                                # get read strand
                                read_strand = get_read_strand(pileupread.alignment)
                                # make a function to classify read as ref or alt based on strand, base, ref, alt
                                read_allele = classify_read_allele(
                                    read_base, read_strand, snp_ref, snp_alt
                                )
                                # add to classification counter
                                classification_counter[read_allele] += 1
                                # append to list if not already in list
                                if read_allele == "ref":
                                    if read_id not in ref_merged_block_reads_list:
                                        snp_ref_merged_block_reads_list.append(read_id)
                                elif read_allele == "alt":
                                    if read_id not in alt_merged_block_reads_list:
                                        snp_alt_merged_block_reads_list.append(read_id)
                # if ref and alt reads are more than set read threshold, and vaf is within vaf threshold, add to ref_merged_block_reads_list/alt_merged_block_reads_list
                ref_coverage = len(snp_ref_merged_block_reads_list)
                alt_coverage = len(snp_alt_merged_block_reads_list)
                total_coverage = ref_coverage + alt_coverage
                if (
                    ref_coverage >= read_coverage_threshold
                    and alt_coverage >= read_coverage_threshold
                ):
                    snp_vaf = alt_coverage / total_coverage
                    if snp_vaf >= vaf_threshold and snp_vaf <= 1 - vaf_threshold:
                        ref_merged_block_reads_list += snp_ref_merged_block_reads_list
                        alt_merged_block_reads_list += snp_alt_merged_block_reads_list
                        # add to snps_dict (ref, alt, vaf, coverage)
                        het_snps_dict[snp_id] = {
                            "chrom": snp_chrom,
                            "pos": snp_pos,
                            "ref": snp_ref,
                            "alt": snp_alt,
                            "vaf": snp_vaf,
                            "coverage": total_coverage,
                            "regions_id": regions_id,
                        }
                        het_counter += 1
            # log progress
            if all_counter % 500 == 0:
                logging.info(
                    f"Processed {all_counter} SNPs, {het_counter} valid heterozygous ones found"
                )

        # add to dict
        # drop any read that is in both ref and alt
        shared_reads = set(ref_merged_block_reads_list).intersection(
            set(alt_merged_block_reads_list)
        )
        ref_merged_block_reads_list = [
            x for x in ref_merged_block_reads_list if x not in shared_reads
        ]
        alt_merged_block_reads_list = [
            x for x in alt_merged_block_reads_list if x not in shared_reads
        ]
        # add to dict if not empty
        if (
            len(ref_merged_block_reads_list) > 0
            and len(alt_merged_block_reads_list) > 0
        ):
            ref_merged_block_reads_dict[merged_block][
                "ref"
            ] = ref_merged_block_reads_list
            ref_merged_block_reads_dict[merged_block][
                "alt"
            ] = alt_merged_block_reads_list

    # close bam file
    bamfile.close()

    # drop any block that has no reads
    ref_merged_block_reads_dict = {
        k: v for k, v in ref_merged_block_reads_dict.items() if len(v) > 0
    }

    # open bam file again and write reads to files
    logging.info("Writing reads to files")
    bamfile = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_fasta)
    outbam_ref = pysam.AlignmentFile(
        f"{results_path}/{out_prefix}.ref.bam",
        "wb",
        template=bamfile,
    )
    outbam_alt = pysam.AlignmentFile(
        f"{results_path}/{out_prefix}.alt.bam",
        "wb",
        template=bamfile,
    )
    bams_dict = {
        "ref": outbam_ref,
        "alt": outbam_alt,
    }

    # write reads to files
    for merged_block in vcf_with_regions["merged_block_id"].unique():
        merged_block_df = vcf_with_regions[
            vcf_with_regions["merged_block_id"] == merged_block
        ]
        added_reads = {
            "ref": {
                "read1": [],
                "read2": [],
                "unpaired": [],
            },
            "alt": {
                "read1": [],
                "read2": [],
                "unpaired": [],
            },
        }
        for regions_idx, regions_row in merged_block_df.iterrows():
            if merged_block not in ref_merged_block_reads_dict:
                continue
            regions_id = regions_row["regions_id"]
            regions_chrom = regions_row["regions_chrom"]
            regions_start = regions_row["regions_start"]
            regions_end = regions_row["regions_end"]
            # iterate over snps in block
            regions_snps = merged_block_df[merged_block_df["regions_id"] == regions_id]
            for snp_idx, snp_row in regions_snps.iterrows():
                snp_chrom = snp_row["snp_chrom"]
                snp_pos = snp_row["snp_pos"]
                snp_ref = snp_row["ref"]
                snp_alt = snp_row["alt"]
                # iterate over reads pileup at snp position
                for pileup in bamfile.pileup(snp_chrom, snp_pos - 1, snp_pos):
                    if pileup.pos == snp_pos - 1:
                        for pileupread in pileup.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                # get read id, strand, mate number
                                read_id = pileupread.alignment.qname
                                read_pair = (
                                    "read1"
                                    if pileupread.alignment.is_read1
                                    else (
                                        "read2"
                                        if pileupread.alignment.is_read2
                                        else "unpaired"
                                    )
                                )
                                if (
                                    read_id
                                    in ref_merged_block_reads_dict[merged_block]["ref"]
                                ):
                                    read_allele = "ref"
                                elif (
                                    read_id
                                    in ref_merged_block_reads_dict[merged_block]["alt"]
                                ):
                                    read_allele = "alt"
                                else:
                                    continue
                                # write read to file
                                if read_id not in added_reads[read_allele][read_pair]:
                                    added_reads[read_allele][read_pair].append(read_id)
                                    bams_dict[read_allele].write(pileupread.alignment)

    # close bam files
    bamfile.close()
    outbam_ref.close()
    outbam_alt.close()

    # sort and index bam files using samtools and os.system
    logging.info("Sorting and indexing BAM files")
    os.system(
        f"samtools sort {results_path}/{out_prefix}.ref.bam > {results_path}/{out_prefix}.ref.sorted.bam"
    )
    os.system(f"samtools index {results_path}/{out_prefix}.ref.sorted.bam")
    os.system(
        f"samtools sort {results_path}/{out_prefix}.alt.bam > {results_path}/{out_prefix}.alt.sorted.bam"
    )
    os.system(f"samtools index {results_path}/{out_prefix}.alt.sorted.bam")

    # remove unsorted bam files
    os.remove(f"{results_path}/{out_prefix}.ref.bam")
    os.remove(f"{results_path}/{out_prefix}.alt.bam")

    # make df of hets and write to file
    het_snps_df = pd.DataFrame(het_snps_dict).T
    het_snps_df.reset_index(inplace=True, drop=False)
    het_snps_df.columns = [
        "snp_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "vaf",
        "coverage",
        "regions_id",
    ]
    # reorder columns
    het_snps_df = het_snps_df[
        ["chrom", "pos", "snp_id", "ref", "alt", "vaf", "coverage", "regions_id"]
    ]

    # write to file
    logging.info("Writing SNPs to file")
    het_snps_df.to_csv(
        f"{results_path}/{out_prefix}.hets.snps.tsv",
        sep="\t",
        index=False,
        header=True,
    )

    # done
    logging.info("Done")


if __name__ == "__main__":
    raise RuntimeError(
        "This module is not meant to be run directly. Use the CLI entry point."
    )
