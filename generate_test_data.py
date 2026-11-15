"""Generate test FASTQ files for quality control."""
from pathlib import Path
import gzip

# Define the output directory for the generated test FASTQ files.
output_dir = Path("data/input")
output_dir.mkdir(parents=True, exist_ok=True)

# Define the test read pairs for each sample.
samples = {
    "sample_001": [
        ("ACGTACGTACGTACGT", "TGCATGCATGCATGCA"),
        ("GGGAAACCCGGGAAAC", "CCCTTTGGGCCCTTTG"),
    ],
    "sample_002": [
        ("ATATCGCGATATCGCG", "CGCGATATCGCGATAT"),
        ("AAAACCCCGGGGTTTT", "TTTTGGGGCCCCAAAA"),
    ],
}

# Generate the test FASTQ files for each sample.
for sample_id, read_pairs in samples.items():
    r1_path = output_dir / f"{sample_id}_R1.fastq.gz"
    r2_path = output_dir / f"{sample_id}_R2.fastq.gz"

    with gzip.open(r1_path, "wt") as r1, \
         gzip.open(r2_path, "wt") as r2:

        for i, (sequence_1, sequence_2) in enumerate(read_pairs, start=1):
            quality_1 = "I" * len(sequence_1)
            quality_2 = "I" * len(sequence_2)

            r1.write(
                f"@{sample_id}_{i}/1\n"
                f"{sequence_1}\n"
                "+\n"
                f"{quality_1}\n"
            )

            r2.write(
                f"@{sample_id}_{i}/2\n"
                f"{sequence_2}\n"
                "+\n"
                f"{quality_2}\n"
            )

print(f"Generated test FASTQ files in {output_dir}")