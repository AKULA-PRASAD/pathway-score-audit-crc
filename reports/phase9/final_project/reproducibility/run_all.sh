#!/bin/zsh
# Reproduce every Phase 9 final-project output from the five frozen checkpoints.
# Writes ONLY under reports/phase9/final_project/. Never re-runs the bootstrap.
# Usage (from the repository root):  zsh reports/phase9/final_project/reproducibility/run_all.sh
set -euo pipefail
cd "$(dirname "$0")/../../../.."           # repository root
PY=.venv/bin/python
FP=reports/phase9/final_project
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/p9_mpl"
if pgrep -f '40_phase9_dc_interval' >/dev/null; then echo "a Phase 9 driver is running -- refusing to start"; exit 2; fi

echo "== S0 starting-state verification";  $PY -B -u $FP/analysis/s0_verify_starting_state.py > $FP/audit/s0_run.log
echo "== S1 section 9.1 analysis";          $PY -B -u $FP/analysis/s1_interval_analysis.py > $FP/audit/s1_run.log
echo "== independent re-computation";      $PY -B -u $FP/validation/independent_check.py | tee $FP/audit/independent_check_run.log
echo "== CV seed-state evidence";           $PY -B -I $FP/validation/cv_seed_state_evidence.py > $FP/audit/cv_seed_state_evidence_run.log
echo "== S2 figures";                       $PY -B -u $FP/analysis/s2_figures.py > $FP/audit/s2_run.log
echo "== S3 documents";                     $PY -B -u $FP/analysis/s3_build_docs.py | tee $FP/audit/s3_run.log
echo "== tests (unit, integration, documents)"
$PY -B -m pytest -p no:cacheprovider -c /dev/null $FP/validation/tests -q 2>&1 | tee $FP/audit/pytest_final_project.log
echo "== final QC (pass 1)";                $PY -B -u $FP/analysis/s9_final_qc.py > $FP/audit/final_qc_pass1.log || true
echo "== S3 re-render with QC result";      $PY -B -u $FP/analysis/s3_build_docs.py >> $FP/audit/s3_run.log
echo "== document tests on final render";   $PY -B -m pytest -p no:cacheprovider -c /dev/null $FP/validation/tests -q 2>&1 | tee $FP/audit/pytest_final_project.log
echo "== final QC (pass 2)";                $PY -B -u $FP/analysis/s9_final_qc.py | tee $FP/audit/final_qc_run.log
