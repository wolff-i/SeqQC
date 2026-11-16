# SeqQC

SeqQC is a Python and Bash pipeline for validating
paired-end sequencing files and generating automated
quality-control reports. This is an independent 
demonstration project, not a validated clinical 
diagnostic workflow.

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

## Quick Start

Generate synthetic paired-end FASTQ test data:

```bash
python generate_test_data.py
```

Run the pipeline:

```bash
bash scripts/run_qc.sh
```

SeqQC creates a new output directory containing
the FastQC reports, aggregated MultiQC report,
CSV summary, log file, and `_SUCCESS` marker.

The synthetic FASTQ files are intended for testing
the workflow, not evaluating real sequencing quality.

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

## Testing

Run unit tests:

```bash
PYTHONPATH=src python -m pytest \
    tests/ -m "not integration" -v
```

Run the end-to-end integration test:

```bash
PYTHONPATH=src python -m pytest \
    tests/ -m integration -v
```

The integration test requires FastQC and MultiQC
to be installed and available on PATH.

Unit tests use temporary fixtures and mocked external
commands to verify validation, error handling,
reporting, and pipeline orchestration.
