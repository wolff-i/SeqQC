"""Tests for FASTQ pair discovery and validation."""
import gzip
import pytest
from seqqc.validation import (
    discover_pairs,
    validate_fastq_pair,
)


def write_fastq(path, reads):
    """Create a small gzipped FASTQ file for testing."""

    with gzip.open(path, "wt", encoding="ascii") as handle:
        for read_id, sequence in reads:
            quality = "I" * len(sequence)

            handle.write(
                f"@{read_id}\n"
                f"{sequence}\n"
                "+\n"
                f"{quality}\n"
            )


def test_discover_pairs(tmp_path):
    """Test that FASTQ pairs are discovered correctly."""
    write_fastq(
        tmp_path / "sample_A_R1.fastq.gz",
        [("r1/1", "ACGT")]
    )
    write_fastq(
        tmp_path / "sample_A_R2.fastq.gz",
        [("r1/2", "TGCA")]
    )

    pairs = discover_pairs(tmp_path)

    # Ensure that the discovered pairs match the expected sample and file names
    assert list(pairs) == ["sample_A"]
    assert (
        pairs["sample_A"][0].name
        == "sample_A_R1.fastq.gz"
    )


def test_missing_mate(tmp_path):
    """Test that missing mates are rejected."""
    write_fastq(
        tmp_path / "sample_A_R1.fastq.gz",
        [("r1/1", "ACGT")]
    )

    with pytest.raises(ValueError, match="missing R2"):
        discover_pairs(tmp_path)


def test_valid_pair(tmp_path):
    """Test that a valid FASTQ pair is validated correctly."""
    r1 = tmp_path / "sample_R1.fastq.gz"
    r2 = tmp_path / "sample_R2.fastq.gz"

    # Write valid FASTQ files for both mates
    write_fastq(
        r1,
        [("read1/1", "ACGT"), ("read2/1", "TTTT")]
    )
    write_fastq(
        r2,
        [("read1/2", "TGCA"), ("read2/2", "AAAA")]
    )

    assert validate_fastq_pair(r1, r2) == 2


def test_mismatched_read_ids(tmp_path):
    """Test that mismatched read IDs are rejected."""
    r1 = tmp_path / "R1.gz"
    r2 = tmp_path / "R2.gz"

    # Write FASTQ files with mismatched read IDs
    write_fastq(r1, [("read1/1", "ACGT")])
    write_fastq(r2, [("read2/2", "TGCA")])

    with pytest.raises(
        ValueError,
        match="Read IDs do not match"
    ):
        validate_fastq_pair(r1, r2)


def test_unequal_read_counts(tmp_path):
    """Test that unequal read counts are rejected."""
    r1 = tmp_path / "R1.gz"
    r2 = tmp_path / "R2.gz"

    # Write FASTQ files with unequal read counts
    write_fastq(
        r1,
        [("read1/1", "ACGT"), ("read2/1", "TTTT")]
    )
    write_fastq(
        r2,
        [("read1/2", "TGCA")]
    )

    with pytest.raises(
        ValueError,
        match="different read counts"
    ):
        validate_fastq_pair(r1, r2)


def test_malformed_fastq(tmp_path):
    """Test that malformed FASTQ records are rejected."""
    r1 = tmp_path / "R1.gz"
    r2 = tmp_path / "R2.gz"

    # Write a malformed FASTQ file for R1 and a valid FASTQ file for R2
    with gzip.open(r1, "wt") as handle:
        handle.write("@read1/1\nACGT\n+\nII\n")
    write_fastq(r2, [("read1/2", "TGCA")])

    with pytest.raises(
        ValueError,
        match="FASTQ validation failed"
    ):
        validate_fastq_pair(r1, r2)


def test_corrupted_gzip(tmp_path):
    """Test that corrupted gzip files are rejected."""
    r1 = tmp_path / "R1.gz"
    r2 = tmp_path / "R2.gz"

    # Write a corrupted gzip file for R1 and a valid FASTQ file for R2
    r1.write_bytes(b"not a gzip file")
    write_fastq(r2, [("read1/2", "TGCA")])

    with pytest.raises(
        ValueError,
        match="FASTQ validation failed"
    ):
        validate_fastq_pair(r1, r2)


def test_no_records(tmp_path):
    """Test that empty FASTQ files are rejected."""
    r1 = tmp_path / "R1.gz"
    r2 = tmp_path / "R2.gz"

    # Write empty FASTQ files for both mates
    write_fastq(r1, [])
    write_fastq(r2, [])

    with pytest.raises(
        ValueError,
        match="No FASTQ records"
    ):
        validate_fastq_pair(r1, r2)
