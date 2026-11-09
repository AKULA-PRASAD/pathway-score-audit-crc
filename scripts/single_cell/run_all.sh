#!/usr/bin/env bash
# Re-run the registered single-cell extension S1-S5 (docs/SINGLE_CELL_EXTENSION_PLAN.md). Outcome-free.
# Needs the public GEO files in data/external/single_cell/ (see their PROVENANCE.json) and the .venv-sc environment
# (versions in results/single_cell/ENVIRONMENT.txt).
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=.venv-sc/bin/python
$PY scripts/single_cell/s0_build.py
$PY scripts/single_cell/s1_tests.py
$PY scripts/single_cell/s_validate.py
$PY scripts/single_cell/s_figures.py
