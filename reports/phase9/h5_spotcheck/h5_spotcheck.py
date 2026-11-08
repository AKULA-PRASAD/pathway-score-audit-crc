"""Phase 9 -- H5 reproducibility spot-check (plan v2.1 section 8.4, hard failure H5).

Option A (user-approved 2026-11-04): recompute run-0 outer replicates
k = 0, 25, ..., 475 with the FROZEN driver scripts/40_phase9_dc_interval.py exactly
as used for the original five runs (its own seeding, its own outer_replicate), and
require bit-identical dC_corrected (float.hex equality). Any mismatch -> H5 FAIL.

The driver is loaded with importlib; its main() is NOT invoked, so the driver writes
nothing. This script reads the original checkpoints and writes ONLY into
reports/phase9/h5_spotcheck/.

Usage:  .venv/bin/python -u reports/phase9/h5_spotcheck/h5_spotcheck.py
"""
from __future__ import annotations

import hashlib, importlib.util, json, os, sys, time
import multiprocessing as mp
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
DRIVER = ROOT / "scripts" / "40_phase9_dc_interval.py"
CKPT_DIR = ROOT / "reports" / "phase9" / "dc_interval"

SELECTED_K = list(range(0, 500, 25))          # section 8.4 fixed rule: 0, 25, ..., 475
RUN_ID, B_IN = 0, 100

# 17-file verification (driver FROZEN set + cohort membership + driver + plan + audit records)
EXPECTED_FILES = {
    "docs/PHASE9_MODELING_PROTOCOL_DRAFT.md": "ab56708919f3746f",
    "config/pathways/m2_panel_4set.json": "43ce15b2f0e06f75",
    "config/pathways/M2_FREEZE_RECORD.json": "64a86fa8ed0d6a8a",
    "config/analysis/A19_FREEZE_RECORD.json": "a77de7b37fe09eeb",
    "src/pathway_scoring.py": "38cac035c0b2d450",
    "src/score_standardisation.py": "41b5c091125e1b11",
    "src/phase9_execution.py": "571d0c233f367be6",
    "data/processed/analysis/expression_median_S0.parquet": "a1ca41807c534584",
    "data/processed/analysis/covariates.csv": "47b2d9829a4fccec",
    "data/processed/analysis/endpoint.csv": "802a68575a34a3b6",
    "reports/phase9/execution/phase9_execution_results.json": "05d13b8b339fa7d4",
    "reports/phase9/execution_s7/phase9_s7_results.json": "b9202305bb6416e1",
    "data/processed/analysis/cohort_membership.csv": "bd3671f8df64c594",
    "scripts/40_phase9_dc_interval.py": "4b1738a4d0ee0e39",
    "docs/PHASE9_DELTAC_INTERVAL_PLAN.md": "d27c52a29930a338",
    "reports/phase9/dc_interval/HALT_RECORD.json": "c1b404309bd63744",
    "reports/phase9/dc_interval/RESUME_AUTHORIZATION.json": "0571c10e02e9059d",
}
# checkpoint hashes established at completion of each run (2026-11-03/04)
EXPECTED_CKPT = {
    "run0_outer.jsonl": "885bbd1bd891b570",
    "run1_outer.jsonl": "1bb27c56612a5f9d",
    "run2_outer.jsonl": "0f12dce12bfce95e",
    "run3_outer.jsonl": "4115bd4cd014d345",
    "run4_outer.jsonl": "2202d559b0418392",
}
THREAD_VARS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "P9_WORKERS")


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify(drv, stage: str) -> dict:
    out = {"stage": stage, "captured_utc": now()}
    try:
        drv.check_frozen()
        out["driver_check_frozen"] = "PASS"
    except SystemExit as ex:
        out["driver_check_frozen"] = f"FAIL: {ex}"
    files = {}
    for rel, exp in EXPECTED_FILES.items():
        full = sha(ROOT / rel)
        files[rel] = {"expected_16": exp, "sha256": full, "ok": full.startswith(exp)}
    out["files_17"] = files
    ck = {}
    for name, exp in EXPECTED_CKPT.items():
        full = sha(CKPT_DIR / name)
        n = sum(1 for l in (CKPT_DIR / name).read_text().splitlines() if l.strip())
        ck[name] = {"expected_16": exp, "sha256": full, "records": n, "ok": full.startswith(exp)}
    out["checkpoints_5"] = ck
    out["all_ok"] = (out["driver_check_frozen"] == "PASS"
                     and all(v["ok"] for v in files.values())
                     and all(v["ok"] for v in ck.values()))
    return out


def load_driver():
    spec = importlib.util.spec_from_file_location("p9_dc_driver_frozen", DRIVER)
    drv = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = drv           # so forked workers can unpickle drv.outer_replicate
    spec.loader.exec_module(drv)           # defines functions only; main() is NOT called
    return drv


def main() -> int:
    t0 = time.time()
    drv = load_driver()
    from src.reproducibility import environment_report  # driver put ROOT on sys.path

    env = {"thread_env": {v: os.environ.get(v) for v in THREAD_VARS},
           "driver_N_WORKERS": drv.N_WORKERS, "driver_BASE": drv.BASE,
           "driver_path": str(DRIVER), "environment": environment_report()}

    # ---------- pre-execution verification ----------
    pre = verify(drv, "pre")
    orig = {}
    for line in (CKPT_DIR / "run0_outer.jsonl").read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            if r["k"] in SELECTED_K:
                orig[r["k"]] = r
    pre["selected_k"] = SELECTED_K
    pre["selected_present_in_run0"] = sorted(orig) == SELECTED_K
    pre["selected_status"] = {k: orig[k]["status"] for k in sorted(orig)}
    pre["all_ok"] = pre["all_ok"] and pre["selected_present_in_run0"]
    pre["execution_config"] = env
    (HERE / "pre_verification.json").write_text(json.dumps(pre, indent=1))
    print(f"[h5] pre-verification all_ok={pre['all_ok']}  "
          f"check_frozen={pre['driver_check_frozen']}  "
          f"selected_present={pre['selected_present_in_run0']}", flush=True)
    if not pre["all_ok"]:
        print("[h5] PRE-VERIFICATION FAILED -- not recomputing", flush=True)
        return 2

    # ---------- recomputation: same configuration as the original driver ----------
    todo = [(RUN_ID, k, B_IN) for k in SELECTED_K]
    drv.init_worker()                                  # load once before fork, as the driver does
    ctx = mp.get_context("fork")
    recs = {}
    with ctx.Pool(drv.N_WORKERS) as pool, open(HERE / "recomputed_records.jsonl", "w") as fh:
        for i, rec in enumerate(pool.imap_unordered(drv.outer_replicate, todo, chunksize=1), 1):
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            recs[rec["k"]] = rec
            print(f"[h5] {i}/{len(todo)} recomputed k={rec['k']}  {(time.time()-t0)/60:.1f} min",
                  flush=True)

    # ---------- comparison ----------
    rows = []
    for k in SELECTED_K:
        o, n = orig[k], json.loads(json.dumps(recs[k]))   # same JSON round-trip as the checkpoint
        ss = drv.np.random.SeedSequence(entropy=[drv.BASE, RUN_ID]).spawn(k + 1)[k]
        oh = float(o["dC_corrected"]).hex() if "dC_corrected" in o else None
        nh = float(n["dC_corrected"]).hex() if "dC_corrected" in n else None
        diff_keys = sorted(x for x in set(o) | set(n) if o.get(x) != n.get(x))
        rows.append({
            "k": k, "seed_entropy": [drv.BASE, RUN_ID], "spawn_key": list(ss.spawn_key),
            "dC_corrected_original": o.get("dC_corrected"),
            "dC_corrected_recomputed": n.get("dC_corrected"),
            "dC_corrected_original_hex": oh, "dC_corrected_recomputed_hex": nh,
            "dC_corrected_bit_identical": oh is not None and oh == nh,
            "multiplicity_sha_match": o["multiplicity_sha"] == n["multiplicity_sha"],
            "full_record_identical": not diff_keys, "differing_fields": diff_keys,
            "status_original": o["status"], "status_recomputed": n["status"],
        })
    n_match = sum(r["dC_corrected_bit_identical"] for r in rows)

    # ---------- post-execution verification ----------
    post = verify(drv, "post")
    (HERE / "post_verification.json").write_text(json.dumps(post, indent=1))

    verdict = "PASS" if (n_match == len(SELECTED_K) and pre["all_ok"] and post["all_ok"]) else "FAIL"
    result = {
        "check": "H5 reproducibility spot-check (plan v2.1 section 8.4 / section 11 H5)",
        "option": "A -- frozen driver, original seeding, outer_replicate((0, k, 100))",
        "criterion": "dC_corrected bit-identical (float.hex equality) for all 20 selected k",
        "completed_utc": now(), "runtime_min": round((time.time() - t0) / 60, 2),
        "n_selected": len(SELECTED_K), "n_bit_identical": n_match,
        "n_mismatch": len(SELECTED_K) - n_match,
        "n_full_record_identical": sum(r["full_record_identical"] for r in rows),
        "n_multiplicity_sha_match": sum(r["multiplicity_sha_match"] for r in rows),
        "pre_all_ok": pre["all_ok"], "post_all_ok": post["all_ok"],
        "script_sha256": sha(Path(__file__).resolve()),
        "recomputed_records_sha256": sha(HERE / "recomputed_records.jsonl"),
        "H5_verdict": verdict, "comparisons": rows,
    }
    (HERE / "h5_comparison.json").write_text(json.dumps(result, indent=1))
    print(f"[h5] bit-identical {n_match}/{len(SELECTED_K)}  post all_ok={post['all_ok']}  "
          f"H5 {verdict}", flush=True)
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
