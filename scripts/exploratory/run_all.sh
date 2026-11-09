#!/usr/bin/env bash
# Re-run every registered exploratory analysis (docs/EXPLORATORY_ANALYSIS_PLAN.md) in order.
# Requires the project .venv and Rscript with the base `survival` package. Reads frozen inputs only.
set -euo pipefail
cd "$(dirname "$0")"
PY=../../.venv/bin/python
for s in x1_molecular_anchoring x2_competing_risks x3_batch x4_mmr_provenance x5_composition \
         x6_scores_vs_composition x7_composition_adjusted_cox x8_external_audit x_figures x_validate; do
  echo "== $s"; $PY "$s.py" > /dev/null
done
echo "done: outputs in reports/exploratory/"
