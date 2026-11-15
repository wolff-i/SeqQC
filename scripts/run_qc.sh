#!/usr/bin/env bash

set -Eeuo pipefail

# Find the project root regardless of where the script is executed.
PROJECT_ROOT="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd
)"

# Ensure the active environment has Python.
if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: python not found; activate the seqqc Conda environment." >&2
    exit 127
fi

# Make the project's Python package importable.
export PYTHONPATH="${PROJECT_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"

# Forward every argument to the Python entry point.
exec python -m seqqc.main "$@"
