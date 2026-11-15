
"""End-to-end testing with actual FastQC and MultiQC."""
import csv
import gzip
import shutil
import subprocess
from pathlib import Path
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def generate_fastq_pair(input_dir: Path) -> None:
    """Generate small paired-end FASTQ fixtures."""

    input_dir.mkdir()

    for mate, sequence in (
        ("R1", "ACGT" * 20),
        ("R2", "TGCA" * 20),
    ):
        path = input_dir / f"demo_{mate}.fastq.gz"

        with gzip.open(path, "wt") as handle:
            for index in range(1, 11):
                handle.write(
                    f"@read_{index:04d}/{mate[-1]}\n"
                    f"{sequence}\n"
                    "+\n"
                    f"{'I' * len(sequence)}\n"
                )


@pytest.mark.integration
def test_end_to_end_pipeline(tmp_path):
    """Run real QC tools and verify final outputs."""

    # Skip the test if the required QC tools are not installed.
    for tool in ("fastqc", "multiqc"):
        if shutil.which(tool) is None:
            pytest.skip(
                f"{tool} is not installed in this environment"
            )

    input_dir = tmp_path / "input reads"
    output_dir = tmp_path / "results"
    
    # Generate the paired-end FASTQ fixtures.
    generate_fastq_pair(input_dir)

    # Build the QC pipeline command.
    command = [
        "bash",
        str(PROJECT_ROOT / "scripts/run_qc.sh"),
        "--input", str(input_dir),
        "--output", str(output_dir),
        "--threads", "1",
    ]

    # Run the QC pipeline.
    result = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        timeout=180,
    )

    # Verify that the pipeline completed successfully.
    assert result.returncode == 0, (
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )

    # Verify that the success marker and run log were created.
    assert (output_dir / "_SUCCESS").is_file()
    assert (output_dir / "run.log").is_file()

    # Verify that the FastQC reports were created.
    for mate in ("R1", "R2"):
        prefix = f"demo_{mate}_fastqc"

        for extension in (".html", ".zip"):
            path = (
                output_dir / "fastqc"
                / f"{prefix}{extension}"
            )
            assert path.is_file()
            assert path.stat().st_size > 0

    # Verify that the MultiQC report was created.
    multiqc_report = (
        output_dir / "multiqc"
        / "multiqc_report.html"
    )

    assert multiqc_report.is_file()
    assert multiqc_report.stat().st_size > 0

    # Verify that the QC summary contains both mates.
    summary_path = output_dir / "qc_summary.csv"

    with summary_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 2
    assert {row["mate"] for row in rows} == {
        "R1", "R2"
    }
    assert all(
        int(row["total_sequences"]) == 10
        for row in rows
    )
