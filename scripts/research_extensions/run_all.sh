#!/usr/bin/env bash
# Re-run registered research extensions E1–E3 (docs/RESEARCH_EXTENSION_REGISTER.md). Outcome-free; frozen inputs read-only.
set -euo pipefail
cd "$(dirname "$0")"
PY=../../.venv/bin/python
for s in e1_gene_architecture e1b_clusterwise_stability e2_conditioned_structure e3_replication e3b_positive_control_qc e3c_tcga e3d_tcga_challenge e_challenge e_figures e_validate; do
  echo "== $s"; $PY "$s.py" > /dev/null
done
