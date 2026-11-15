"""Summarize FastQC ZIP results as a machine-readable CSV."""
import argparse
import csv
from pathlib import Path
from zipfile import BadZipFile, ZipFile
from seqqc.validation import discover_pairs


# Define the columns for the quality control summary CSV.
COLUMNS = [
    "sample_id",
    "mate",
    "fastq_file",
    "total_sequences",
    "sequence_length",
    "gc_percent",
    "pass_modules",
    "warn_modules",
    "fail_modules",
    "warned_checks",
    "failed_checks",
]


def _read_member(archive: ZipFile, member: str) -> str:
    """Read a uniquely named FastQC report file from a ZIP."""
    # Find the uniquely named report member.
    matches = [
        name for name in archive.namelist()
        if name.endswith(f"/{member}")
    ]
    # Verify that exactly one matching member was found.
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {member}; found {len(matches)}"
        )
    # Read the matching report member as UTF-8 text.
    return archive.read(matches[0]).decode("utf-8")


def parse_fastqc_zip(zip_path: Path, expected_filename: str,) -> dict:
    """Extract basic statistics and module status counts."""

    zip_path = Path(zip_path)

    # Read the FastQC report files from the ZIP archive.
    try:
        with ZipFile(zip_path) as archive:
            data_text = _read_member(
                archive, "fastqc_data.txt"
            )
            summary_text = _read_member(
                archive, "summary.txt"
            )

    except (
        OSError,
        BadZipFile,
        UnicodeError,
        RuntimeError,
        ValueError,
    ) as exc:
        raise ValueError(
            f"Cannot read FastQC archive {zip_path}: {exc}"
        ) from exc

    # Extract values from the Basic Statistics section.
    statistics = {}
    in_basic = False

    # Parse the Basic Statistics section of the FastQC data file.
    for line in data_text.splitlines():

        if line.startswith(">>Basic Statistics\t"):
            in_basic = True
            continue

        if in_basic and line == ">>END_MODULE":
            break

        if in_basic and line and not line.startswith("#"):
            key, separator, value = line.partition("\t")

            if not separator:
                raise ValueError(
                    f"Malformed Basic Statistics in {zip_path}"
                )

            statistics[key] = value

    # Extract the required values from the Basic Statistics section.
    try:
        filename = statistics["Filename"]
        total = int(statistics["Total Sequences"])
        length = statistics["Sequence length"]
        gc_percent = float(statistics["%GC"])

    except (KeyError, ValueError) as exc:
        raise ValueError(
            f"Missing or invalid Basic Statistics in {zip_path}"
        ) from exc

    # Verify that the reported filename matches the expected filename.
    if Path(filename).name != expected_filename:
        raise ValueError(
            f"Filename mismatch in {zip_path}: "
            f"{filename} != {expected_filename}"
        )

    # Verify that the reported metrics are within valid ranges.
    if total < 0 or not 0 <= gc_percent <= 100:
        raise ValueError(
            f"Invalid FastQC metrics in {zip_path}"
        )

    # Parse the PASS / WARN / FAIL summary.
    statuses = {
        "PASS": [],
        "WARN": [],
        "FAIL": [],
    }

    # Parse each summary line and group modules by status.
    for line in summary_text.splitlines():

        if not line.strip():
            continue

        parts = line.split("\t")

        if len(parts) != 3:
            raise ValueError(
                f"Malformed summary.txt in {zip_path}"
            )

        status, module, reported_filename = parts

        if status not in statuses:
            raise ValueError(
                f"Invalid module status {status} in {zip_path}"
            )

        if Path(reported_filename).name != expected_filename:
            raise ValueError(
                f"Summary filename mismatch in {zip_path}"
            )

        statuses[status].append(module)

    # Verify that at least one module result was found.
    if not any(statuses.values()):
        raise ValueError(
            f"No FastQC module results in {zip_path}"
        )

    # Return the parsed FastQC metrics.
    return {
        "fastq_file": expected_filename,
        "total_sequences": total,
        "sequence_length": length,
        "gc_percent": gc_percent,
        "pass_modules": len(statuses["PASS"]),
        "warn_modules": len(statuses["WARN"]),
        "fail_modules": len(statuses["FAIL"]),
        "warned_checks": "; ".join(statuses["WARN"]),
        "failed_checks": "; ".join(statuses["FAIL"]),
    }


def write_qc_summary(
    pairs: dict[str, tuple[Path, Path]],
    fastqc_dir: Path,
    output_csv: Path,
    expected_counts: dict[str, int] | None = None,
) -> Path:
    """Write one CSV row per FASTQ file."""

    fastqc_dir = Path(fastqc_dir)
    output_csv = Path(output_csv)
    rows = []

    # Verify that at least one paired FASTQ file was provided.
    if not pairs:
        raise ValueError(
            "No paired FASTQ files to summarize"
        )

    # Process each paired FASTQ file in sorted order.
    for sample, (r1, r2) in sorted(pairs.items()):

        for mate, fastq in (("R1", r1), ("R2", r2)):

            # Derive the expected FastQC report filename.
            filename = fastq.name
            prefix = filename.removesuffix(".fastq.gz")
            zip_path = (
                fastqc_dir / f"{prefix}_fastqc.zip"
            )

            # Parse the FastQC metrics from the ZIP report.
            metrics = parse_fastqc_zip(
                zip_path, filename
            )

            # Compare FastQC counts with input validation.
            if (
                expected_counts is not None
                and metrics["total_sequences"]
                != expected_counts[sample]
            ):
                raise ValueError(
                    f"FastQC read count mismatch for {filename}: "
                    f"expected {expected_counts[sample]}, "
                    f"found {metrics['total_sequences']}"
                )

            rows.append({
                "sample_id": sample,
                "mate": mate,
                **metrics,
            })

    # Create the output directory if it does not exist.
    output_csv.parent.mkdir(
        parents=True, exist_ok=True
    )

    # Write the quality control summary to the CSV file.
    with output_csv.open(
        "w", newline="", encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle, fieldnames=COLUMNS
        )

        writer.writeheader()
        writer.writerows(rows)

    return output_csv


def main() -> int:
    # Configure the command-line interface.
    parser = argparse.ArgumentParser(
        description="Create a CSV summary from FastQC ZIP reports"
    )

    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--fastqc", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)

    args = parser.parse_args()

    # Discover paired FASTQ files and write the quality control summary.
    try:
        pairs = discover_pairs(args.input)

        result = write_qc_summary(pairs, args.fastqc, args.output)

    except (OSError, ValueError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")

    print(f"Created QC summary: {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
