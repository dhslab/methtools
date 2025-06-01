import typer
from .snpsplit.snp_split import main as snpsplit_main
from .readmeth.read_meth import main as readmeth_main
from types import SimpleNamespace

app = typer.Typer()


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context):
    """Methtools: A toolkit for methylation analysis"""
    if ctx.invoked_subcommand is None:
        typer.echo("Welcome to Methtools! Use --help for available commands.")


@app.command()
def snpsplit(
    bam_path: str = typer.Option(..., help="Path to input BAM/CRAM file"),
    ref_fasta: str = typer.Option(
        "/storage2/fs1/dspencer/Active/spencerlab/projects/longread_methylation/data/reference/genome/hg38_mgi_patch.fa",
        help="Path to reference FASTA file",
    ),
    vcf_path: str = typer.Option(..., help="Path to input VCF file"),
    regions_path: str = typer.Option(..., help="Path to input Regions BED file"),
    results_path: str = typer.Option(
        "./results", help="Path to output directory (default: ./results)"
    ),
    bp_distance: int = typer.Option(
        2000, help="Distance in bp for merging Regions (default: 2000)"
    ),
    out_prefix: str = typer.Option(
        "test", help="Prefix for output files (default: test)"
    ),
):
    """Run the SNP split analysis"""
    typer.echo("Running SNP split analysis...")
    # Create an args object similar to argparse
    args = SimpleNamespace(
        bam_path=bam_path,
        ref_fasta=ref_fasta,
        vcf_path=vcf_path,
        regions_path=regions_path,
        results_path=results_path,
        bp_distance=bp_distance,
        out_prefix=out_prefix,
    )
    snpsplit_main(args)


@app.command()
def readmeth(
    bam: str = typer.Option(..., help="Comma-separated list of one or two BAM files."),
    haplotypes: str = typer.Option(
        "ref,alt", help="Comma-separated haplotype names for the two BAM files."
    ),
    out_file: str = typer.Option(..., help="Output file path."),
    regions_bed: str = typer.Option(..., help="Regions BED file."),
    sample: str = typer.Option(..., help="Sample name."),
):
    """Extract read level methylation data from BAM files."""
    typer.echo("Running methylation data processing...")
    bam_list = bam.split(",")  # Split the comma-separated BAM file paths into a list
    haplotype_list = haplotypes.split(
        ","
    )  # Split the comma-separated haplotype names into a list
    readmeth_main(
        bam=bam_list,
        haplotypes=haplotype_list,
        out_file=out_file,
        regions_bed=regions_bed,
        sample=sample,
    )


if __name__ == "__main__":
    app()
