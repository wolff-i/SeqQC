"""FASTQ pair discovery and validation."""
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
from itertools import zip_longest
from pathlib import Path
import re
import sys
import zlib


# Regular expression to parse FASTQ file names and extract sample and mate information
FASTQ_NAME = re.compile(r"^(?P<sample>.+)_R(?P<mate>[12])\.fastq\.gz$")


def discover_pairs(input_dir: Path) -> dict[str, tuple[Path, Path]]:
    """Return sample IDs mapped to (R1, R2) FASTQ paths."""

    # Ensure the input directory is a Path object
    input_dir = Path(input_dir)

    # Validate the input directory
    if not input_dir.is_dir():
        raise ValueError(f"Input directory does not exist: {input_dir}")

    # Discover FASTQ files matching the expected naming pattern
    files = sorted(input_dir.glob("*_R[12].fastq.gz"))

    # Ensure that there are FASTQ files to process
    if not files:
        raise ValueError(f"No FASTQ files found in input directory: {input_dir}")

    # Group FASTQ files by sample and mate
    samples: dict[str, dict[int, Path]] = {}

    # Process each FASTQ file and add it to the appropriate sample and mate group
    for file in files:
        match = FASTQ_NAME.match(file.name)
        if match is None:
            raise ValueError(f"Invalid FASTQ file name: {file.name}")

        # Extract the sample and mate information from the file name
        sample = match.group("sample")
        mate = int(match.group("mate"))

        # Ensure that the sample and mate combination is unique
        if sample in samples and mate in samples[sample]:
            raise ValueError(f"Duplicate FASTQ file for sample {sample}, mate {mate}")

        # Add the FASTQ file to the appropriate sample and mate group
        samples.setdefault(sample, {})[mate] = file

    # Validate and create FASTQ pairs for each sample
    pairs: dict[str, tuple[Path, Path]] = {}
    for sample, mates in sorted(samples.items()):
        if set(mates) != {1, 2}:
            missing = 1 if 1 not in mates else 2
            raise ValueError(f"Sample {sample} is missing R{missing}")

        # Add the validated FASTQ pair to the results
        pairs[sample] = (mates[1], mates[2])

    return pairs


def normalize_read_id(title: str) -> str:
    """Extract the normalized read ID from a FASTQ title."""

    read_id = title.split()[0]

    # Remove the trailing /1 or /2 if present
    if read_id.endswith(("/1", "/2")):
        return read_id[:-2]
    
    return read_id


def validate_fastq_pair(r1: Path, r2: Path) -> int:
    """Validate paired FASTQ files and return read-pair count."""

    # Ensure the input paths are Path objects
    r1, r2 = Path(r1), Path(r2)

    # Validate that the input FASTQ files exist
    for path in (r1, r2):
        if not path.is_file():
            raise ValueError(
                f"Missing FASTQ file: {path}"
            )
        # Ensure the input FASTQ files are not empty
        if path.stat().st_size == 0:
            raise ValueError(
                f"Empty compressed file: {path}"
            )

    # Initialize the read-pair count
    count = 0

    # Validate the paired FASTQ records
    try:
        with gzip.open(r1, "rt", encoding="ascii") as handle1, \
             gzip.open(r2, "rt", encoding="ascii") as handle2:

            reads1 = FastqGeneralIterator(handle1)
            reads2 = FastqGeneralIterator(handle2)

            # Validate each read pair
            for index, (record1, record2) in enumerate(
                zip_longest(reads1, reads2),
                start=1
            ):
                # Ensure both read files contain the same number of records
                if record1 is None or record2 is None:
                    raise ValueError(
                        "R1/R2 have different read counts "
                        f"(at record {index})"
                    )
                # Normalize the read IDs
                id1 = normalize_read_id(record1[0])
                id2 = normalize_read_id(record2[0])
                # Ensure the read IDs match between R1 and R2
                if id1 != id2:
                    raise ValueError(
                        "Read IDs do not match at "
                        f"record {index}: {id1} != {id2}"
                    )

                count += 1

    # Handle errors during FASTQ validation
    except (
        OSError,
        EOFError,
        UnicodeError,
        ValueError,
        zlib.error,
    ) as exc:
        raise ValueError(
            f"FASTQ validation failed for "
            f"{r1.name} and {r2.name}: {exc}"
        ) from exc

    # Ensure the FASTQ files contain at least one read pair
    if count == 0:
        raise ValueError(
            f"No FASTQ records in {r1.name} and {r2.name}"
        )

    return count


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate paired-end FASTQ.gz files"
    )

    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Directory containing paired FASTQ.gz files"
    )

    args = parser.parse_args()

    # Discover and validate paired FASTQ files
    try:
        pairs = discover_pairs(args.input)

        for sample, (r1, r2) in pairs.items():
            count = validate_fastq_pair(r1, r2)

            print(
                f"PASS {sample}: {count} read pairs"
            )
    # Handle errors during FASTQ file discovery or validation
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"Validated {len(pairs)} sample(s) successfully"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())