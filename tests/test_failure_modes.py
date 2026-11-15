
"""Test SeqQC behavior when pipeline steps fail."""

import gzip

import pytest

from seqqc import quality_control as qc


def make_input_pair(input_dir, r1_id="read1/1", r2_id="read1/2",):
    """Create a minimal paired-end FASTQ dataset."""

    input_dir.mkdir()

    for mate, read_id in (
        ("R1", r1_id),
        ("R2", r2_id),
    ):
        path = input_dir / f"sample_R{mate[-1]}.fastq.gz"

        with gzip.open(path, "wt") as handle:
            handle.write(
                f"@{read_id}\n"
                "ACGTACGT\n"
                "+\n"
                "IIIIIIII\n"
            )


def test_invalid_input_stops_pipeline(tmp_path, monkeypatch):
    """Stop the pipeline when input is invalid."""

    input_dir = tmp_path / "input"

    make_input_pair(
        input_dir,
        r1_id="read1/1",
        r2_id="wrong_read/2",
    )

    # Replace the QC tool runners with mocks that fail if called.
    def unexpected_execution(*args, **kwargs):
        pytest.fail(
            "QC tools should not run on invalid input"
        )
    
    monkeypatch.setattr(
        qc, "run_fastqc", unexpected_execution
    )
    monkeypatch.setattr(
        qc, "run_multiqc", unexpected_execution
    )

    output_dir = tmp_path / "results"

    # Verify that the pipeline stops before creating output.
    with pytest.raises(
        ValueError,
        match="Read IDs do not match",
    ):
        qc.run_pipeline(input_dir, output_dir)

    assert not output_dir.exists()


def test_fastqc_failure_stops_pipeline(tmp_path, monkeypatch):
    """Stop execution and log the failure if FastQC fails."""

    input_dir = tmp_path / "input"
    make_input_pair(input_dir)

    # Mock the tool lookup so the pipeline reaches FastQC.
    monkeypatch.setattr(
        qc.shutil, "which", lambda tool: "/mock/tool"
    )
    # Replace the QC tool runners with mocks that simulate failure.
    def fail_fastqc(*args, **kwargs):
        raise RuntimeError("simulated FastQC failure")
    # Ensure MultiQC does not run after FastQC fails.
    def unexpected_multiqc(*args, **kwargs):
        pytest.fail(
            "MultiQC must not run after FastQC failure"
        )

    monkeypatch.setattr(
        qc, "run_fastqc", fail_fastqc
    )
    monkeypatch.setattr(
        qc, "run_multiqc", unexpected_multiqc
    )

    output_dir = tmp_path / "results"

    # Verify that the pipeline stops when FastQC fails.
    with pytest.raises(
        RuntimeError,
        match="simulated FastQC failure",
    ):
        qc.run_pipeline(input_dir, output_dir)

    # Verify that the failure is recorded in the run log.
    log_path = output_dir / "run.log"

    assert log_path.exists()
    assert "simulated FastQC failure" in (
        log_path.read_text()
    )
    assert not (output_dir / "_SUCCESS").exists()
    assert not (
        output_dir / "qc_summary.csv"
    ).exists()


def test_multiqc_failure_stops_reporting(tmp_path, monkeypatch):
    """Stop reporting and log the failure if MultiQC fails."""

    input_dir = tmp_path / "input"
    make_input_pair(input_dir)

    # Mock the tool lookup so the pipeline reaches MultiQC.
    monkeypatch.setattr(
        qc.shutil, "which", lambda tool: "/mock/tool"
    )

    # Replace the QC tool runners with mocks that simulate failure.
    monkeypatch.setattr(
        qc, "run_fastqc", lambda *args: None
    )

    def fail_multiqc(*args):
        raise RuntimeError(
            "simulated MultiQC failure"
        )

    # Ensure CSV reporting does not run after MultiQC fails.
    def unexpected_reporting(*args, **kwargs):
        pytest.fail(
            "CSV reporting must not run after "
            "MultiQC failure"
        )

    monkeypatch.setattr(
        qc, "run_multiqc", fail_multiqc
    )
    monkeypatch.setattr(
        qc, "write_qc_summary", unexpected_reporting
    )

    output_dir = tmp_path / "results"

    # Verify that the pipeline stops when MultiQC fails.
    with pytest.raises(
        RuntimeError,
        match="simulated MultiQC failure",
    ):
        qc.run_pipeline(input_dir, output_dir)

    # Verify that the failure is recorded in the run log.
    assert not (output_dir / "_SUCCESS").exists()
    assert not (
        output_dir / "qc_summary.csv"
    ).exists()

    assert "simulated MultiQC failure" in (
        output_dir / "run.log"
    ).read_text()
