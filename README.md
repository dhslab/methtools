# Methtools

Methtools is a Python toolkit for methylation analysis, providing command-line tools and API utilities for handling Whole Genome Bisulfite Sequencing (WGBS) and Oxford Nanopore Technologies (ONT) methylation data.

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/m-mahgoub/methtools.git
    cd methtools
    ```

2.  **Install dependencies:**
    It is recommended to use a virtual environment.
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows use `venv\\Scripts\\activate`
    pip install -r requirements.txt
    ```

3.  **Install Methtools:**
    ```bash
    pip install .
    ```
    Or for development (editable install):
    ```bash
    pip install -e .
    ```

After installation, the `methtools` command should be available in your terminal.

## Prerequisite Tools

Methtools requires the following command-line tools to be installed and available in your PATH:

- **Samtools** (minimum version 1.9)
- **Modkit** (minimum version 0.1.10)
- **Bedtools** (minimum version 2.31.1)

Ensure these tools are installed and meet the minimum version requirements before running workflows.

## Command-Line Interface (CLI)

Methtools offers the following command-line tools:

### 1. `snpsplit`

Splits BAM/CRAM files based on SNP information. This is useful for separating reads based on haplotype.

**Example:**

```bash
methtools snpsplit \
    --bam-path input.bam \
    --ref-fasta reference.fasta \
    --vcf-path variants.vcf \
    --regions-path regions.bed \
    --results-path ./snpsplit_results \
    --out-prefix mysnpsplit
```

### 2. `readmeth`

Extracts read-level methylation data from BAM files.

**Example:**

```bash
methtools readmeth \
    --bam hap1.bam,hap2.bam \
    --haplotypes HAP1,HAP2 \
    --out-file methylation_calls.tsv \
    --regions-bed regions_to_analyze.bed \
    --sample MySample
```

**Example (single BAM):**
```bash
methtools readmeth \
    --bam merged_sample.bam \
    --out-file methylation_calls.tsv \
    --regions-bed regions_to_analyze.bed \
    --sample MySample
```

### 3. `call-mhb`

Calls Methylation Haplotype Blocks (MHBs) from BAM files or precomputed methylation call files.

**Note:** When running `call_mhb` with BAM files (using the `--bam` option), `samtools` (>=1.9) and `modkit` (>=0.1.10) must be installed and available in your system's PATH. These tools are used for processing BAM files and extracting methylation calls.

**Example (from BAM files):**

```bash
methtools call-mhb \
    --bam input1.bam,input2.bam \
    --fasta reference.fasta \
    --regions regions_for_mhb.bed \
    --cpgs all_cpg_sites.bed \
    --output-dir ./mhb_output \
    --threads 4 \
    --output-prefix my_mhb_run \
    --output-formats tsv,bed
```

**Example (from precomputed methylation calls):**

```bash
methtools call-mhb \
    --methylation-calls precomputed_calls.tsv \
    --regions regions_for_mhb.bed \
    --cpgs all_cpg_sites.bed \
    --output-dir ./mhb_output \
    --threads 4 \
    --output-prefix my_mhb_run_precomputed
```

## API Utilities

Methtools also provides underlying Python modules and functions that can be imported and used in custom scripts for more advanced or specific methylation data analyses. These can be found within the `methtools`.
