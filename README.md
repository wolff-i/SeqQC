# SeqQC

SeqQC is a Python and Bash pipeline for validating
paired-end sequencing files and generating automated
quality-control reports.

## Features

- Paired-end FASTQ.gz input validation
- gzip integrity and FASTQ structure checks
- Read count and identifier validation
- Automated FastQC execution
- Aggregated MultiQC reports
- CSV summaries of sequencing QC metrics
- Configurable paths and thread counts
- Run-specific output directories and logs
- Automated tests using pytest

## Requirements

- Linux or WSL2
- Conda
- Python 3.15
- Biopython
- FastQC
- MultiQC
- PyYAML
- pytest

## Installation

Create the environment:

    conda env create -f environment.yml
    conda activate seqqc

## Usage

Run with default settings:

    bash scripts/run_qc.sh

Override settings:

    bash scripts/run_qc.sh \
        --input data/input \
        --output data/output/custom_run \
        --threads 4

Edit `config/config.yaml` to change the defaults.

## Input Format

Input files must use the following naming pattern:

    sample_001_R1.fastq.gz
    sample_001_R2.fastq.gz

Reads must be paired and appear in corresponding order.

## Outputs

Each run produces:

- Individual FastQC HTML and ZIP reports
- An aggregated MultiQC HTML report
- `qc_summary.csv`
- `run.log`
