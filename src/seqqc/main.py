"""Configurable command-line entry point for SeqQC."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import sys
from uuid import uuid4
import yaml
from seqqc.quality_control import run_pipeline


# Project root directory for resolving relative paths
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def resolve_path(value: str) -> Path:
    """Resolve paths relative to the SeqQC project root."""
    # Ensure the provided path is a nonempty string.
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Path must be a nonempty string")
    # Expand any user home directory references.
    path = Path(value).expanduser()
    # Resolve relative paths against the project root.
    if not path.is_absolute():
        path = PROJECT_ROOT / path

    return path.resolve()


def load_config(config_path: Path) -> dict:
    """Load and validate the YAML configuration."""

    # Load the YAML configuration.
    try:
        with config_path.open(encoding="utf-8") as handle:
            config = yaml.safe_load(handle)

    except yaml.YAMLError as exc:
        raise ValueError(
            f"Invalid YAML in {config_path}: {exc}"
        ) from exc

    # Ensure the configuration is a mapping.
    if not isinstance(config, dict):
        raise ValueError(
            "Configuration must be a YAML mapping"
        )

    # Validate the configuration keys.
    keys = {"input_dir", "output_root", "threads"}

    missing = keys - config.keys()
    unknown = config.keys() - keys

    if missing:
        raise ValueError(
            "Missing configuration keys: "
            + ", ".join(sorted(missing))
        )

    if unknown:
        raise ValueError(
            "Unknown configuration keys: "
            + ", ".join(sorted(unknown))
        )

    # Validate the configured directory paths.
    for key in ("input_dir", "output_root"):
        if (
            not isinstance(config[key], str)
            or not config[key].strip()
        ):
            raise ValueError(
                f"{key} must be a nonempty path string"
            )

    # Validate the configured thread count.
    threads = config["threads"]

    if (
        isinstance(threads, bool)
        or not isinstance(threads, int)
        or threads < 1
    ):
        raise ValueError(
            "threads must be a positive integer"
        )

    return config


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and execute the SeqQC pipeline."""

    parser = argparse.ArgumentParser(
        description=(
            "Run SeqQC: validate FASTQs, run "
            "FastQC/MultiQC, create CSV"
        )
    )

    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="YAML configuration path",
    )

    parser.add_argument(
        "--input",
        help="Override the configured input directory",
    )

    parser.add_argument(
        "--output",
        help=(
            "New output directory; if omitted, "
            "auto-generate one"
        ),
    )

    parser.add_argument(
        "--threads",
        type=int,
        help="Override FastQC thread count",
    )

    args = parser.parse_args(argv)

    # Load the configuration and resolve the pipeline paths.
    try:
        config = load_config(resolve_path(args.config))

        input_path = resolve_path(
            args.input
            if args.input is not None
            else config["input_dir"]
        )

        threads = (
            args.threads
            if args.threads is not None
            else config["threads"]
        )

        if threads < 1:
            raise ValueError(
                "--threads must be a positive integer"
            )

        if args.output is not None:
            output_path = resolve_path(args.output)

        else:
            timestamp = datetime.now(
                timezone.utc
            ).strftime("%Y%m%dT%H%M%SZ")

            run_id = (
                f"run_{timestamp}_{uuid4().hex[:8]}"
            )

            output_path = (
                resolve_path(config["output_root"])
                / run_id
            )

        report = run_pipeline(
            input_path,
            output_path,
            threads=threads,
        )

    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Output directory: {output_path}")
    print(f"MultiQC report: {report}")
    print(
        f"QC summary: {output_path / 'qc_summary.csv'}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
