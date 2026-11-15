"""Tests for the SeqQC command-line interface."""
from pathlib import Path

from seqqc import main as cli


def write_config(
    path: Path,
    body: str = (
        "input_dir: data/input\n"
        "output_root: data/output\n"
        "threads: 2\n"
    ),
):
    """Write a YAML configuration file."""
    path.write_text(body, encoding="utf-8")


def test_defaults_generate_unique_run_directory(tmp_path, monkeypatch):
    """Test that default runs create unique output directories."""
    config = tmp_path / "config.yaml"
    write_config(config)

    monkeypatch.setattr(cli, "PROJECT_ROOT", tmp_path)

    calls = []
    # Replace the pipeline with a fake implementation that records its arguments.
    def fake_pipeline(input_dir, output_dir, threads):
        calls.append((input_dir, output_dir, threads))
        return (
            output_dir / "multiqc" / "multiqc_report.html"
        )
    # Replace the pipeline with the fake implementation.
    monkeypatch.setattr(
        cli, "run_pipeline", fake_pipeline
    )
    # Run the CLI twice to ensure unique output directories are created.
    assert cli.main(["--config", str(config)]) == 0
    assert cli.main(["--config", str(config)]) == 0

    assert calls[0][0] == tmp_path / "data/input"
    assert calls[0][2] == 2
    assert calls[0][1].parent == tmp_path / "data/output"
    assert calls[0][1].name.startswith("run_")
    assert calls[0][1] != calls[1][1]


def test_cli_overrides_config(tmp_path, monkeypatch):
    """Test that command-line arguments override configuration values."""
    config = tmp_path / "config.yaml"
    write_config(config)

    monkeypatch.setattr(cli, "PROJECT_ROOT", tmp_path)
    
    captured = {}
    # Replace the pipeline with a fake implementation that records its arguments.
    def fake_pipeline(input_dir, output_dir, threads):
        captured.update(
            input_dir=input_dir,
            output_dir=output_dir,
            threads=threads,
        )
        return (
            output_dir / "multiqc" / "multiqc_report.html"
        )

    # Replace the pipeline with the fake implementation.
    monkeypatch.setattr(cli, "run_pipeline", fake_pipeline)

    # Run the CLI with command-line overrides.
    status = cli.main([
        "--config", str(config),
        "--input", "my-input",
        "--output", "my-results/run_01",
        "--threads", "4",
    ])

    assert status == 0

    assert captured == {
        "input_dir": tmp_path / "my-input",
        "output_dir": tmp_path / "my-results/run_01",
        "threads": 4,
    }


def test_malformed_config_rejected(tmp_path, monkeypatch, capsys):
    """Test that malformed configuration files are rejected."""
    config = tmp_path / "config.yaml"
    # Write a malformed configuration file.
    write_config(config, "- not_a_mapping\n")

    monkeypatch.setattr(
        cli, "run_pipeline", lambda *args, **kwargs: None
    )

    # Run the CLI with the malformed configuration.
    assert cli.main(["--config", str(config)]) == 1
    # Verify that the error message is displayed.
    assert (
        "Configuration must be a YAML mapping"
        in capsys.readouterr().err
    )


def test_nonpositive_threads_rejected(tmp_path, monkeypatch, capsys):
    """Test that nonpositive thread counts are rejected."""
    config = tmp_path / "config.yaml"
    write_config(config)

    # Replace the pipeline with a fake implementation.
    monkeypatch.setattr(
        cli, "run_pipeline", lambda *args, **kwargs: None
    )

    # Run the CLI with a nonpositive thread count.
    assert cli.main([
        "--config", str(config),
        "--threads", "0",
    ]) == 1

    # Verify that the error message is displayed.
    assert (
        "--threads must be a positive integer"
        in capsys.readouterr().err
    )


def test_missing_config_rejected(tmp_path, monkeypatch, capsys):
    """Test that missing configuration files are rejected."""
    # Replace the pipeline with a fake implementation.
    monkeypatch.setattr(
        cli, "run_pipeline", lambda *args, **kwargs: None
    )

    # Run the CLI with a missing configuration file.
    assert cli.main([
        "--config", str(tmp_path / "missing.yaml")
    ]) == 1

    # Verify that the error message is displayed.
    assert "ERROR:" in capsys.readouterr().err
