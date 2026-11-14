"""Execute FastQC and aggregate sequencing reports with MultiQC."""
import argparse
import logging
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from seqqc.validation import discover_pairs, validate_fastq_pair

# Configure logging for the module.
logger = logging.getLogger(__name__)


def configure_logging(log_path: Path) -> None:
    """Write progress to terminal and detailed logs to file."""

    # Set the logging level and prevent propagation to the root logger.
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    # Remove any existing handlers to avoid duplicate log entries.
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()

    # Define the log message format.
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s %(message)s"
    )

    # Add a console handler for terminal output.
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    console.setFormatter(formatter)

    # Add a file handler for detailed log output.
    logfile = logging.FileHandler(log_path, encoding="utf-8")
    logfile.setLevel(logging.DEBUG)
    logfile.setFormatter(formatter)

    logger.addHandler(console)
    logger.addHandler(logfile)


def run_command(command: list[str]) -> None:
    """Execute a command and raise an error if it fails."""
    
    # Log the command before execution.
    logger.info("Running: %s", shlex.join(command))

    # Execute the command and capture its output.
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )

    # Handle a missing executable.
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Executable not found: {command[0]}"
        ) from exc

    # Handle a command that exits with an error.
    except subprocess.CalledProcessError as exc:
        detail = (
            exc.stderr or exc.stdout or
            "No error details supplied"
        ).strip()

        raise RuntimeError(
            f"{command[0]} failed with exit code "
            f"{exc.returncode}: {detail}"
        ) from exc

    # Log the command's standard output.
    if result.stdout:
        logger.debug("stdout:\n%s", result.stdout.rstrip())

    # Log the command's standard error.
    if result.stderr:
        logger.debug("stderr:\n%s", result.stderr.rstrip())


def run_fastqc(r1: Path, r2: Path, fastqc_dir: Path, threads: int,) -> None:
    """Run FastQC and verify reports for both mates."""

    # Build the FastQC command.
    command = [
        "fastqc",
        "--threads", str(threads),
        "--outdir", str(fastqc_dir),
        str(r1),
        str(r2),
    ]

    # Execute the FastQC command.
    run_command(command)

    # Verify that FastQC created a report for each mate.
    for fastq in (r1, r2):
        prefix = fastq.name.removesuffix(".fastq.gz")
        
        # Verify each expected report file.
        for suffix in (".html", ".zip"):
            report = (
                fastqc_dir / f"{prefix}_fastqc{suffix}"
            )
            
            # Verify that the report exists and is not empty.
            if (
                not report.is_file()
                or report.stat().st_size == 0
            ):
                raise RuntimeError(
                    "FastQC did not create a valid report: "
                    f"{report}"
                )


def run_multiqc(fastqc_dir: Path,multiqc_dir: Path,) -> Path:
    """Aggregate FastQC reports and return HTML path."""

    # Define the expected MultiQC report path.
    report = multiqc_dir / "multiqc_report.html"

    # Build the MultiQC command.
    command = [
        "multiqc",
        str(fastqc_dir),
        "--module", "fastqc",
        "--require-logs",
        "--outdir", str(multiqc_dir),
        "--filename", report.name,
    ]

    # Execute the MultiQC command.
    run_command(command)

    # Verify that the MultiQC report exists and is not empty.
    if (
        not report.is_file()
        or report.stat().st_size == 0
    ):
        raise RuntimeError(
            f"MultiQC report is missing or empty: {report}"
        )

    return report


def run_pipeline(input_dir: Path, output_dir: Path, threads: int = 2,) -> Path:
    """Validate FASTQ files, run FastQC, then MultiQC."""
    # Convert input and output paths to Path objects.
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    # Verify that the requested thread count is valid.
    if threads < 1:
        raise ValueError("--threads must be at least 1")

    # Verify that the output directory does not already exist.
    if output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists: {output_dir}. "
            "Choose a new directory for each run."
        )

    # Validate every sample before running QC tools.
    pairs = discover_pairs(input_dir)
    counts = {}

    # Validate each sample's FASTQ pair.
    for sample, (r1, r2) in pairs.items():
        counts[sample] = validate_fastq_pair(r1, r2)

    # Check that external tools are installed.
    for tool in ("fastqc", "multiqc"):
        if shutil.which(tool) is None:
            raise RuntimeError(
                f"Required tool not found on PATH: {tool}"
            )

    # Create fresh directories for this run.
    output_dir.mkdir(parents=True, exist_ok=False)

    # Define the QC output directories.
    fastqc_dir = output_dir / "fastqc"
    multiqc_dir = output_dir / "multiqc"

    # Create the QC output directories.
    fastqc_dir.mkdir()
    multiqc_dir.mkdir()

    # Configure logging for the run.
    configure_logging(output_dir / "run.log")

    # Log the number of validated samples.
    logger.info("Validated %d sample(s)", len(pairs))

    # Process each sample.
    for sample, (r1, r2) in pairs.items():
        logger.info(
            "Processing %s (%d read pairs)",
            sample,
            counts[sample],
        )

        # Run FastQC on the sample's FASTQ pair.
        run_fastqc(r1, r2, fastqc_dir, threads)

    # Combine results.
    report = run_multiqc(fastqc_dir, multiqc_dir)

    # Log the completed pipeline.
    logger.info("Completed: %s", report)

    return report


def main() -> int:
    # Configure the command-line interface.
    parser = argparse.ArgumentParser(
        description=(
            "Validate paired FASTQ files and run "
            "FastQC/MultiQC"
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Directory containing FASTQ.gz files",
    )

    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="New output directory for this run",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=2,
        help="FastQC threads (default: 2)",
    )

    args = parser.parse_args()

    # Run the pipeline.
    try:
        report = run_pipeline(
            args.input,
            args.output,
            args.threads,
        )

    # Handle pipeline errors.
    except (OSError, ValueError, RuntimeError) as exc:
        if logger.handlers:
            logger.error("Pipeline failed: %s", exc)

        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    # Print the MultiQC report path.
    print(f"MultiQC report: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
