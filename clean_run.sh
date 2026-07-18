#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
rm -rf .venv-repro
python3 -m venv .venv-repro
. .venv-repro/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-lock.txt
python scripts/06_validate_file_manifest.py
python scripts/run_pipeline.py
echo "Clean run completed. Open generated/reproducibility_validation_report.md"
