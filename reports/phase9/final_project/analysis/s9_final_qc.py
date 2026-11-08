"""Final independent quality-control audit (Section 14 of the work order). Read-only except for
audit/final_qc.json and reproducibility/MANIFEST_SHA256.txt.

Checks
  computational  -- 17 frozen/audit hashes, 5 checkpoint hashes, driver check_frozen(), H5 files,
                    v1 deviation draft; no file outside final_project/ modified after work began
  statistical    -- results computed on the verified checkpoints; independent check all-agree
  documentation  -- every deliverable present; no unrendered placeholders; PROVISIONAL labels
  testing        -- pytest logs show 0 failures
Usage: .venv/bin/python -B -u reports/phase9/final_project/analysis/s9_final_qc.py
"""
from __future__ import annotations

import importlib.util, json, os, sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FP, ROOT, CKPT_DIR, EXPECTED_FILES, EXPECTED_CKPT, sha256  # noqa: E402

H5_DIR = ROOT / "reports/phase9/h5_spotcheck"
H5_FILES = {"h5_spotcheck.py": "41bb774de3cec0a6", "recomputed_records.jsonl": "59266c194947158b",
            "h5_comparison.json": "35e2ed1e8c271da5", "CV_SEED_DEVIATION_DRAFT.md": "2866571cecf471af",
            "pre_verification.json": None, "post_verification.json": None, "h5_run.log": None}
DELIVERABLES = [
    "README.md", "PROJECT_COMPLETION_REPORT.md", "audit/EXECUTION_LEDGER.md",
    "documentation/PHASE9_FINAL_REPORT.md", "documentation/CV_SEED_DEVIATION_DRAFT_v2.md",
    "documentation/STATISTICAL_CONVENTIONS.md", "presentation/PHASE9_TECHNICAL_PRESENTATION.html",
    "results/s0_starting_state.json", "results/s1_interval_results.json",
    "validation/independent_check_results.json", "validation/cv_seed_state_evidence.json",
    "reproducibility/README_REPRODUCE.md", "reproducibility/run_all.sh", "reproducibility/environment.json",
    "figures/FIGURE_CAPTIONS.md", *[f"tables/{t}" for t in (
        "T1_intervals_per_run.csv", "T2_kernel_MCSE.csv", "T3_between_run.csv", "T4_EP1_grid.csv",
        "T5_warnings_hardfailures.csv", "T6_W4_kernel_vs_between.csv")],
    *[f"figures/F{i}_{n}.png" for i, n in enumerate(("outer_distributions", "interval_comparison",
      "endpoint_variability", "MC1_running_trace", "EP1_grid", "shift_decomposition"), 1)]]
EXCLUDE_DIRS = {".venv", ".git", "__pycache__", ".pytest_cache", "final_project"}


def main() -> int:
    qc = {"captured_utc": datetime.now(timezone.utc).isoformat(), "checks": {}}
    c = qc["checks"]
    spec = importlib.util.spec_from_file_location("p9_driver_qc", ROOT / "scripts/40_phase9_dc_interval.py")
    drv = importlib.util.module_from_spec(spec); spec.loader.exec_module(drv)
    try:
        drv.check_frozen(); c["driver_check_frozen"] = True
    except SystemExit:
        c["driver_check_frozen"] = False
    c["files_17_unchanged"] = all(sha256(ROOT / k).startswith(v) for k, v in EXPECTED_FILES.items())
    c["checkpoints_5_unchanged"] = all(sha256(CKPT_DIR / f"run{r}_outer.jsonl").startswith(h) for r, h in EXPECTED_CKPT.items())
    c["h5_evidence_unchanged"] = all(sha256(H5_DIR / n).startswith(h) for n, h in H5_FILES.items() if h)
    c["h5_directory_files"] = sorted(p.name for p in H5_DIR.iterdir())
    c["dc_interval_directory_files"] = sorted(p.name for p in CKPT_DIR.iterdir())
    c["dc_interval_no_new_files"] = c["dc_interval_directory_files"] == sorted(
        ["HALT_RECORD.json", "RESUME_AUTHORIZATION.json", "resume_driver.log"] + [f"run{r}_outer.jsonl" for r in range(5)])

    # nothing outside final_project/ modified since this work order started (S0 timestamp)
    # write-once timestamp of the FIRST S0 run of this work order (re-running S0 does not move it)
    start = datetime.fromisoformat((FP / "audit/WORK_ORDER_START_UTC.txt").read_text().strip()).timestamp() - 60
    c["work_order_start_utc"] = (FP / "audit/WORK_ORDER_START_UTC.txt").read_text().strip()
    touched = []
    for dp, dns, fns in os.walk(ROOT):
        dns[:] = [d for d in dns if d not in EXCLUDE_DIRS]
        for fn in fns:
            p = Path(dp) / fn
            if p.stat().st_mtime > start:
                touched.append(str(p.relative_to(ROOT)))
    c["files_outside_final_project_modified_since_start"] = touched
    c["no_outside_modification"] = not touched

    res = json.loads((FP / "results/s1_interval_results.json").read_text())
    c["results_on_verified_checkpoints"] = {int(k): v for k, v in res["checkpoint_hashes_16"].items()} == EXPECTED_CKPT
    ic = json.loads((FP / "validation/independent_check_results.json").read_text())
    c["independent_check"] = f"{ic['n_ok']}/{ic['n_checks']}"
    c["independent_check_all_ok"] = ic["all_ok"]

    missing = [d for d in DELIVERABLES if not (FP / d).exists()]
    c["deliverables_missing"] = missing
    docs = [FP / "documentation/PHASE9_FINAL_REPORT.md", FP / "presentation/PHASE9_TECHNICAL_PRESENTATION.html",
            FP / "PROJECT_COMPLETION_REPORT.md"]
    c["docs_no_placeholders"] = all("{{" not in d.read_text() for d in docs if d.exists())
    c["docs_labelled_provisional"] = all("PROVISIONAL" in d.read_text() for d in docs if d.exists())
    log = (FP / "audit/pytest_final_project.log").read_text() if (FP / "audit/pytest_final_project.log").exists() else ""
    c["pytest_final_project_log"] = log.strip().splitlines()[-1] if log.strip() else "missing"
    c["pytest_no_failures"] = "passed" in log and "failed" not in log and "error" not in log.lower()
    reg = (FP / "audit/pytest_project_regression.log").read_text()
    c["project_regression_log"] = reg.strip().splitlines()[-1]

    keys_bool = [k for k, v in c.items() if isinstance(v, bool)]
    qc["all_ok"] = all(c[k] for k in keys_bool) and not missing
    qc["failed"] = [k for k in keys_bool if not c[k]] + (["deliverables_missing"] if missing else [])
    (FP / "audit/final_qc.json").write_text(json.dumps(qc, indent=1))

    # manifest: every final_project file + referenced originals
    lines = ["# sha256  path   (final_project files, then referenced original artifacts; read-only)"]
    for p in sorted(FP.rglob("*")):
        if p.is_file() and "__pycache__" not in p.parts and p.name != "MANIFEST_SHA256.txt":
            lines.append(f"{sha256(p)}  {p.relative_to(ROOT)}")
    for rel in [*EXPECTED_FILES, *[f"reports/phase9/dc_interval/run{r}_outer.jsonl" for r in range(5)],
                "reports/phase9/dc_interval/resume_driver.log", "scripts/38_phase9_execute.py",
                *[f"reports/phase9/h5_spotcheck/{n}" for n in H5_FILES]]:
        lines.append(f"{sha256(ROOT / rel)}  {rel}")
    (FP / "reproducibility/MANIFEST_SHA256.txt").write_text("\n".join(lines) + "\n")
    print(json.dumps({"all_ok": qc["all_ok"], "failed": qc["failed"], **{k: c[k] for k in keys_bool},
                      "touched": touched, "missing": missing, "independent_check": c["independent_check"],
                      "pytest": c["pytest_final_project_log"], "regression": c["project_regression_log"]}, indent=1))
    return 0 if qc["all_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
