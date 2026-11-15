"""Tests for FastQC and MultiQC quality control."""
import gzip
import subprocess
import pytest
from seqqc import quality_control as qc


# Create a valid FASTQ record for testing.
def write_fastq(path, read_id):
    """Create one valid FASTQ record."""

    with gzip.open(path, "wt") as handle:
        handle.write(
            f"@{read_id}\n"
            "ACGT\n"
            "+\n"
            "IIII\n"
        )


# Test that a failed external command raises RuntimeError.
def test_run_command_failure(monkeypatch):
    """A failed external command raises RuntimeError."""
    # Define a command runner that simulates a failure.
    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(
            2,
            args[0],
            stderr="out of memory",
        )

    # Replace the subprocess runner with the failing command.
    monkeypatch.setattr(qc.subprocess, "run", fail)

    # Verify that the command failure is reported.
    with pytest.raises(
        RuntimeError,
        match="exit code 2: out of memory",
    ):
        qc.run_command(["fastqc", "some.fastq.gz"])


def test_fastqc_checks_outputs(tmp_path, monkeypatch):
    """A successful command without reports still fails."""

    # Define the input FASTQ files and output directory.
    r1 = tmp_path / "sample_R1.fastq.gz"
    r2 = tmp_path / "sample_R2.fastq.gz"

    fastqc_dir = tmp_path / "fastqc"
    fastqc_dir.mkdir()

    # Replace the command runner with a no-op.
    monkeypatch.setattr(
        qc,
        "run_command",
        lambda command: None,
    )

    # Verify that missing FastQC reports are detected.
    with pytest.raises(
        RuntimeError,
        match="did not create a valid report",
    ):
        qc.run_fastqc(r1, r2, fastqc_dir, 2)


def test_multiqc_checks_report(tmp_path, monkeypatch):
    """MultiQC must create a nonempty HTML report."""

    # Define the FastQC and MultiQC output directories.
    fastqc_dir = tmp_path / "fastqc"
    multiqc_dir = tmp_path / "multiqc"

    fastqc_dir.mkdir()
    multiqc_dir.mkdir()

    # Replace the command runner with a no-op.
    monkeypatch.setattr(
        qc,
        "run_command",
        lambda command: None,
    )

    # Verify that a missing MultiQC report is detected.
    with pytest.raises(
        RuntimeError,
        match="missing or empty",
    ):
        qc.run_multiqc(fastqc_dir, multiqc_dir)


def test_pipeline_with_mocked_tools(tmp_path, monkeypatch):
    """Test the full orchestration without external tools."""

    # Define the input directory and create test FASTQ files.
    inputs = tmp_path / "input"
    inputs.mkdir()

    write_fastq(
        inputs / "sample_001_R1.fastq.gz",
        "read1/1",
    )
    write_fastq(
        inputs / "sample_001_R2.fastq.gz",
        "read1/2",
    )

    # Define the output directory and record executed commands.
    output = tmp_path / "run_001"
    calls = []

    # Replace the executable lookup with a mock path.
    monkeypatch.setattr(
        qc.shutil,
        "which",
        lambda tool: "/mock/" + tool,
    )

    # Replace the command runner with a mock that creates reports.
    def fake_execute(command):
        calls.append(command[0])

        if command[0] == "fastqc":
            fastqc_dir = output / "fastqc"

            for mate in ("R1", "R2"):
                prefix = f"sample_001_{mate}_fastqc"

                (fastqc_dir / f"{prefix}.html").write_text(
                    "html"
                )
                (fastqc_dir / f"{prefix}.zip").write_bytes(
                    b"zip"
                )

        else:
            report = (
                output / "multiqc" / "multiqc_report.html"
            )
            report.write_text("html")

    monkeypatch.setattr(
        qc,
        "run_command",
        fake_execute,
    )

    # Replace the QC summary writer with a mock.
    def fake_write_summary(
        pairs,
        fastqc_dir,
        output_csv,
        expected_counts=None,
    ):
        output_csv.write_text(
            "sample_id,mate\nsample_001,R1\n"
        )
        return output_csv

    monkeypatch.setattr(
        qc,
        "write_qc_summary",
        fake_write_summary,
    )
    # Run the full quality control pipeline.
    report = qc.run_pipeline(
        inputs,
        output,
        threads=2,
    )

    # Verify that the report and expected commands were created.
    assert report.exists()
    assert calls == ["fastqc", "multiqc"]
    assert (output / "run.log").exists()
    assert (output / "qc_summary.csv").exists()


def test_existing_output_rejected(tmp_path):
    """Prevent mixing reports from different runs."""

    # Define the existing output directory.
    output = tmp_path / "run_001"
    output.mkdir()

    # Verify that an existing output directory is rejected.
    with pytest.raises(
        FileExistsError,
        match="already exists",
    ):
        qc.run_pipeline(tmp_path, output)
