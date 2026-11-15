"""Tests for FastQC report parsing and summary generation."""
import csv
from pathlib import Path
from zipfile import ZipFile
import pytest
from seqqc.reporting import (
    parse_fastqc_zip,
    write_qc_summary,
)


def make_fastqc_zip(
    path,
    filename,
    total=2,
    length="16",
    gc="50",
    with_summary=True,
    with_basic=True,
):
    """Create a minimal FastQC ZIP archive for testing."""

    # Start the FastQC data file with the required header.
    data = "##FastQC\t0.12.1\n"

    # Add the Basic Statistics module when requested.
    if with_basic:
        data += (
            ">>Basic Statistics\tpass\n"
            "#Measure\tValue\n"
            f"Filename\t{filename}\n"
            f"Total Sequences\t{total}\n"
            f"Sequence length\t{length}\n"
            f"%GC\t{gc}\n"
            ">>END_MODULE\n"
        )

    # Define the FastQC summary contents.
    summary = (
        f"PASS\tBasic Statistics\t{filename}\n"
        f"WARN\tPer base sequence quality\t{filename}\n"
        f"FAIL\tAdapter Content\t{filename}\n"
    )

    # Create the FastQC ZIP archive.
    with ZipFile(path, "w") as archive:
        archive.writestr(
            f"{path.stem}/fastqc_data.txt",
            data,
        )
        # Add the summary file when requested.
        if with_summary:
            archive.writestr(
                f"{path.stem}/summary.txt",
                summary,
            )


def test_parse_fastqc_zip(tmp_path):
    """Parse a valid FastQC ZIP archive."""
    path = tmp_path / "sample_R1_fastqc.zip"

    # Create a valid FastQC ZIP archive for testing.
    make_fastqc_zip(
        path,
        "sample_R1.fastq.gz",
        length="16-20",
    )

    # Parse the FastQC ZIP archive.
    result = parse_fastqc_zip(
        path, "sample_R1.fastq.gz"
    )

    # Verify the parsed FastQC report values.
    assert result["total_sequences"] == 2
    assert result["sequence_length"] == "16-20"
    assert result["gc_percent"] == 50.0

    assert (
        result["pass_modules"],
        result["warn_modules"],
        result["fail_modules"],
    ) == (1, 1, 1)

    assert result["failed_checks"] == "Adapter Content"


def test_write_qc_summary(tmp_path):
    """Write a summary from parsed FastQC reports."""
    files = (
        Path("sample_R1.fastq.gz"),
        Path("sample_R2.fastq.gz"),
    )

    # Create FastQC ZIP archives for each FASTQ file.
    for fastq in files:
        prefix = fastq.name.removesuffix(".fastq.gz")

        make_fastqc_zip(
            tmp_path / f"{prefix}_fastqc.zip",
            fastq.name,
        )

    # Define the summary output file.
    csv_file = tmp_path / "summary.csv"

    # Write the summary from the parsed FastQC reports.
    write_qc_summary(
        {"sample": files},
        tmp_path,
        csv_file,
        expected_counts={"sample": 2},
    )

    # Read the generated summary file.
    with csv_file.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    # Verify the generated summary contents.
    assert len(rows) == 2
    assert {row["mate"] for row in rows} == {"R1", "R2"}

    assert all(
        row["total_sequences"] == "2"
        for row in rows
    )


def test_corrupted_zip(tmp_path):
    """Reject a corrupted FastQC ZIP archive."""
    # Create a corrupted FastQC ZIP archive.
    path = tmp_path / "bad_fastqc.zip"
    path.write_bytes(b"not a zip")

    # Verify that the corrupted FastQC ZIP archive is rejected.
    with pytest.raises(
        ValueError,
        match="Cannot read FastQC archive",
    ):
        parse_fastqc_zip(path, "bad.fastq.gz")


def test_missing_summary(tmp_path):
    """Reject a FastQC ZIP archive without a summary."""
    path = tmp_path / "test_fastqc.zip"

    # Create a FastQC ZIP archive without a summary.
    make_fastqc_zip(
        path,
        "test.fastq.gz",
        with_summary=False,
    )

    # Verify that the FastQC ZIP archive without a summary is rejected.
    with pytest.raises(
        ValueError,
        match="summary.txt",
    ):
        parse_fastqc_zip(path, "test.fastq.gz")


def test_read_count_mismatch(tmp_path):
    """Reject a summary when expected read counts do not match."""
    files = (
        Path("sample_R1.fastq.gz"),
        Path("sample_R2.fastq.gz"),
    )

    # Create FastQC ZIP archives for each FASTQ file.
    for fastq in files:
        prefix = fastq.name.removesuffix(".fastq.gz")

        make_fastqc_zip(
            tmp_path / f"{prefix}_fastqc.zip",
            fastq.name,
        )

    # Define the summary output file.
    path = tmp_path / "summary.csv"

    # Verify that the summary is rejected when expected read counts do not match.
    with pytest.raises(
        ValueError,
        match="read count mismatch",
    ):
        write_qc_summary(
            {"sample": files},
            tmp_path,
            path,
            expected_counts={"sample": 3},
        )

    # Verify that the summary file was not created.
    assert not path.exists()


def test_missing_basic_statistics(tmp_path):
    """Reject a FastQC ZIP archive without Basic Statistics."""
    path = tmp_path / "test_fastqc.zip"

    # Create a FastQC ZIP archive without Basic Statistics.
    make_fastqc_zip(
        path,
        "test.fastq.gz",
        with_basic=False,
    )

    # Verify that the FastQC ZIP archive without Basic Statistics is rejected.
    with pytest.raises(
        ValueError,
        match="Missing or invalid Basic Statistics",
    ):
        parse_fastqc_zip(path, "test.fastq.gz")
