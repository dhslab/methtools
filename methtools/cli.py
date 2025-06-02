import typer
from pathlib import Path
from typing import List, Optional, Union
import logging # Added for logging levels

from .snpsplit.snp_split import main as snpsplit_main
from .readmeth.read_meth import main as readmeth_main
from types import SimpleNamespace

# Import workflow functions and logging setup
from .mhb.workflows import mhb_workflow_from_bams, mhb_workflow_from_precomputed_calls
from .utils.logconfig import setup_logging


app = typer.Typer(context_settings={"help_option_names": ["-h", "--help"]})


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Enable verbose logging (DEBUG level). Default is WARNING and above.", show_default=False)
):
    """
    Methtools: A toolkit for methylation analysis.
    
    Default logging shows WARNING and ERROR messages.
    Use -v or --verbose for more detailed log output (DEBUG level).
    """
    if verbose:
        setup_logging(level=logging.DEBUG)
        logging.debug("Verbose logging enabled (DEBUG level).")
    else:
        # If not verbose, set logging to WARNING level
        # This will show WARNING, ERROR, and CRITICAL messages only.
        setup_logging(level=logging.WARNING)

    if ctx.invoked_subcommand is None:
        typer.echo("Welcome to Methtools! Use 'methtools --help' for available commands.")


@app.command()
def snpsplit(
    bam_path: Path = typer.Option(..., help="Path to input BAM/CRAM file", exists=True, file_okay=True, dir_okay=False, readable=True),
    ref_fasta: Path = typer.Option(
        ..., 
        help="Path to reference FASTA file",
        exists=True, file_okay=True, dir_okay=False, readable=True
    ),
    vcf_path: Path = typer.Option(..., help="Path to input VCF file", exists=True, file_okay=True, dir_okay=False, readable=True),
    regions_path: Path = typer.Option(..., help="Path to input Regions BED file", exists=True, file_okay=True, dir_okay=False, readable=True),
    results_path: Path = typer.Option(
        Path("./results_snpsplit"), help="Path to output directory.", writable=True
    ),
    bp_distance: int = typer.Option(
        2000, help="Distance in bp for merging Regions."
    ),
    out_prefix: str = typer.Option(
        "snpsplit_out", help="Prefix for output files."
    ),
):
    """Split BAM reads by SNPs (WGBS)."""
    logging.info("Running SNP split analysis...") 
    results_path.mkdir(parents=True, exist_ok=True) 
    
    args = SimpleNamespace(
        bam_path=str(bam_path),
        ref_fasta=str(ref_fasta),
        vcf_path=str(vcf_path),
        regions_path=str(regions_path),
        results_path=str(results_path),
        bp_distance=bp_distance,
        out_prefix=out_prefix,
    )
    snpsplit_main(args)
    logging.info("SNP split analysis completed.")


@app.command()
def readmeth(
    bam: str = typer.Option(..., help="Comma-separated list of one or two BAM files. E.g., 'file1.bam,file2.bam' or 'file.bam'."),
    haplotypes: str = typer.Option(
        "ref,alt", help="Comma-separated haplotype names for the two BAM files (if two BAMs provided)."
    ),
    out_file: Path = typer.Option(..., help="Output file path for the methylation data.", writable=True),
    regions_bed: Path = typer.Option(..., help="Regions BED file.", exists=True, file_okay=True, dir_okay=False, readable=True),
    sample: str = typer.Option(..., help="Sample name."),
):
    """Extract read level methylation data from BAM files (WGBS)."""
    logging.info("Running methylation data processing...") 
    
    if out_file.parent and not out_file.parent.exists():
        out_file.parent.mkdir(parents=True, exist_ok=True)

    bam_list_str = bam.split(",")
    bam_list_paths = []
    for b_str in bam_list_str:
        b_path = Path(b_str.strip())
        if not b_path.exists():
            logging.error(f"Input BAM file not found: {b_path}") 
            raise typer.Exit(code=1)
        bam_list_paths.append(b_path)
        
    haplotype_list = haplotypes.split(",")
    
    readmeth_main(
        bam=[str(p) for p in bam_list_paths], 
        haplotypes=haplotype_list,
        out_file=str(out_file),
        regions_bed=str(regions_bed),
        sample=sample,
    )
    logging.info("Methylation data processing completed.")


@app.command(context_settings={"help_option_names": ["-h", "--help"]})
def call_mhb(
    # Distinguishing Input Options - changed to accept comma-separated strings
    bam_str: Optional[str] = typer.Option(None, "--bam", help="Comma-separated list of input BAM files (e.g., aln1.bam,aln2.bam). Triggers BAM processing workflow.", show_default=False),
    methylation_calls: Optional[Path] = typer.Option(None, "--methylation-calls", help="Path to precomputed methylation call file(s) (single TSV, directory of TSVs, or file listing TSVs). Triggers precomputed calls workflow.", exists=True, readable=True, show_default=False),

    # Options Specific to BAM Workflow
    fasta: Optional[Path] = typer.Option(None, "--fasta", help="Path to the reference FASTA file (required if --bam is used).", exists=True, file_okay=True, dir_okay=False, readable=True, show_default=False),
    methylation_calls_dir: Optional[Path] = typer.Option(None, "--methylation-calls-dir", help="Directory to save extracted methylation call TSVs when using --bam (optional with --bam). If provided, directory will be created if it doesn't exist.", file_okay=False, dir_okay=True, writable=True, show_default=False),

    # Common Required Options
    regions: Path = typer.Option(..., help="BED file defining regions for MHB analysis.", exists=True, file_okay=True, dir_okay=False, readable=True),
    cpgs: Path = typer.Option(..., help="BED file with all CpG positions.", exists=True, file_okay=True, dir_okay=False, readable=True),
    output_dir: Path = typer.Option(..., help="Directory where final MHB results will be saved. Will be created if it doesn't exist.", file_okay=False, dir_okay=True, writable=True),

    # Common Optional Parameters
    temp_dir: Optional[Path] = typer.Option(None, "--temp-dir", help="Custom base directory for temporary files. If not set, a temporary sub-directory within --output-dir will be used.", file_okay=False, dir_okay=True, writable=True, show_default=False),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads for parallel processing."),
    output_prefix: str = typer.Option("mhb_analysis", "--output-prefix", help="Prefix for output result files."),
    output_formats_str: Optional[str] = typer.Option(None, "--output-formats", help="Comma-separated list of desired output formats (e.g., tsv,bed). Defaults to 'tsv,bed' if not specified.", show_default=False), 
    mld_min_reads: int = typer.Option(10, "--mld-min-reads", help="Minimum read coverage for mLD calculation."),
    block_r2_soft: float = typer.Option(0.5, "--block-r2-soft", help="Soft r-squared threshold for block definition."),
    block_r2_hard: float = typer.Option(0.4, "--block-r2-hard", help="Hard r-squared threshold for block definition."),
    block_min_cpgs: int = typer.Option(4, "--block-min-cpgs", help="Minimum number of CpGs to constitute a block."),
    block_gaps_tolerance: int = typer.Option(1, "--block-gaps-tolerance", help="Allowed number of r-squared values between soft and hard thresholds before terminating a block."),
    keep_temp_files: bool = typer.Option(False, "--keep-temp-files", help="Keep temporary files after the run.")
):
    """
    Call Methylation Haplotype Blocks (MHBs) from BAM files (ONT).

    This command can start from BAM files (requiring alignment data and a reference FASTA)
    or from precomputed methylation call files (typically TSV format).
    Specify your input method using either --bam (comma-separated) or --methylation-calls.
    """
    logging.info("Initiating MHB calling workflow...")

    # --- Parse comma-separated inputs ---
    bam_list: Optional[List[Path]] = None
    if bam_str:
        bam_list = [Path(p.strip()) for p in bam_str.split(',')]

    parsed_output_formats: Optional[List[str]] = None
    if output_formats_str:
        parsed_output_formats = [fmt.strip().lower() for fmt in output_formats_str.split(',')]
    # --- Input Validation ---
    if bam_list and methylation_calls:
        logging.error("Error: --bam and --methylation-calls are mutually exclusive. Please provide only one input type.")
        raise typer.Exit(code=1)

    if not bam_list and not methylation_calls:
        logging.error("Error: You must specify an input method. Use either --bam (with --fasta) or --methylation-calls.")
        raise typer.Exit(code=1)

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logging.error(f"Could not create or access output directory: {output_dir}. Error: {e}")
        raise typer.Exit(code=1)
    
    if methylation_calls_dir: 
        try:
            methylation_calls_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logging.error(f"Could not create or access directory for saving methylation calls: {methylation_calls_dir}. Error: {e}")
            raise typer.Exit(code=1)

    success = False
    if bam_list:
        logging.info("Processing with BAM file input.")
        if not fasta:
            logging.error("Error: --fasta is required when using --bam.")
            raise typer.Exit(code=1)
        
        for bam_file_path in bam_list: # bam_list is now a List[Path]
            if not bam_file_path.is_file() or not bam_file_path.stat().st_size > 0 : 
                 logging.error(f"Input BAM file not found or is empty: {bam_file_path}")
                 raise typer.Exit(code=1)

        success = mhb_workflow_from_bams(
            fasta_file_path=fasta,
            bam_file_paths=bam_list, # Pass the parsed list of Paths
            regions_bed_file_path=regions,
            cpg_annotations_bed_path=cpgs,
            output_dir_path=output_dir,
            custom_temp_dir_base=temp_dir,
            num_threads=threads,
            output_prefix=output_prefix,
            output_formats=parsed_output_formats, # Pass the parsed list
            mld_pairing_threshold=mld_min_reads,
            block_r2_soft_threshold=block_r2_soft,
            block_r2_hard_threshold=block_r2_hard,
            block_min_cpgs_in_block=block_min_cpgs,
            block_tolerance_gaps=block_gaps_tolerance,
            keep_temp_files=keep_temp_files,
            save_extracted_calls_to=methylation_calls_dir
        )
    elif methylation_calls:
        logging.info("Processing with precomputed methylation calls input.")
        if fasta:
            logging.warning("Warning: --fasta was provided but is not used when processing with --methylation-calls.")
        if methylation_calls_dir:
            logging.warning("Warning: --methylation-calls-dir was provided but is not used when processing with --methylation-calls (as calls are already extracted).")

        success = mhb_workflow_from_precomputed_calls(
            precomputed_calls_files=methylation_calls,
            regions_bed_file_path=regions,
            cpg_annotations_bed_path=cpgs,
            output_dir_path=output_dir,
            custom_temp_dir_base=temp_dir,
            num_threads=threads,
            output_prefix=output_prefix,
            output_formats=parsed_output_formats, # Pass the parsed list
            mld_pairing_threshold=mld_min_reads,
            block_r2_soft_threshold=block_r2_soft,
            block_r2_hard_threshold=block_r2_hard,
            block_min_cpgs_in_block=block_min_cpgs,
            block_tolerance_gaps=block_gaps_tolerance,
            keep_temp_files=keep_temp_files
        )

    if success:
        logging.info("MHB calling workflow finished successfully.")
    else:
        logging.error("MHB calling workflow failed.")
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()