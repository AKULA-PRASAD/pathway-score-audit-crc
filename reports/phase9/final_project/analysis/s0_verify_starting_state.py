"""Stage 0 -- independent verification of the starting state (read-only).

Checks, without trusting any earlier report:
  * the 17 frozen/audit file hashes and the driver's own check_frozen()
  * the five run checkpoints: hashes, record structure, k coverage, status,
    seeds, paired-complete counts, inner failures, H1 and H6
  * internal arithmetic of every record (dC_app, dOpt, dC_corrected, inner_C)
  * the H5 evidence files: hashes and internal consistency with run0
Writes ONLY reports/phase9/final_project/results/s0_starting_state.json.

Usage: .venv/bin/python -u reports/phase9/final_project/analysis/s0_verify_starting_state.py
"""
from __future__ import annotations

import hashlib, importlib.util, json, math, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
FP = ROOT / "reports" / "phase9" / "final_project"
CK = ROOT / "reports" / "phase9" / "dc_interval"
H5 = ROOT / "reports" / "phase9" / "h5_spotcheck"
sys.path.insert(0, str(FP / "analysis"))
from common import EXPECTED_FILES, EXPECTED_CKPT, BASE, B_OUT, B_IN, sha256  # noqa: E402

H5_FILES = {  # hashes reported at H5 completion, 2026-11-04
    "h5_spotcheck.py": "41bb774de3cec0a6",
    "recomputed_records.jsonl": "59266c194947158b",
    "h5_comparison.json": "35e2ed1e8c271da5",
    "CV_SEED_DEVIATION_DRAFT.md": "2866571cecf471af",
}


def main() -> int:
    out = {"stage": "S0 starting-state verification",
           "captured_utc": datetime.now(timezone.utc).isoformat()}

    # --- driver's own H4 check (driver loaded without running main) ---
    spec = importlib.util.spec_from_file_location("p9_driver_ro", ROOT / "scripts/40_phase9_dc_interval.py")
    drv = importlib.util.module_from_spec(spec); spec.loader.exec_module(drv)
    try:
        drv.check_frozen(); out["driver_check_frozen"] = "PASS"
    except SystemExit as ex:
        out["driver_check_frozen"] = f"FAIL {ex}"

    out["files_17"] = {r: {"expected": e, "observed": sha256(ROOT / r)[:16],
                           "ok": sha256(ROOT / r).startswith(e)} for r, e in EXPECTED_FILES.items()}

    # --- checkpoints ---
    runs = {}
    max_arith = {"dC_app": 0.0, "dOpt_from_inner": 0.0, "dC_corrected": 0.0}
    for r in range(5):
        p = CK / f"run{r}_outer.jsonl"
        raw = p.read_bytes()
        rows = [json.loads(l) for l in raw.decode().splitlines() if l.strip()]
        ks = [x["k"] for x in rows]
        h6_bad, keyset = 0, set()
        for x in rows:
            keyset |= set(x)
            n1, n2 = len(x["inner_C"]["M1"]), len(x["inner_C"]["M2"])
            if not (n1 == n2 == x["B_in_paired_complete"]):
                h6_bad += 1
            # internal arithmetic, recomputed independently
            max_arith["dC_app"] = max(max_arith["dC_app"],
                                      abs((x["C_app"]["M2"] - x["C_app"]["M1"]) - x["dC_app"]))
            o1 = np.array(x["inner_C"]["M1"]); o2 = np.array(x["inner_C"]["M2"])
            dopt = (o2[:, 0] - o2[:, 1]).mean() - (o1[:, 0] - o1[:, 1]).mean()
            max_arith["dOpt_from_inner"] = max(max_arith["dOpt_from_inner"], abs(dopt - x["dOpt"]))
            max_arith["dC_corrected"] = max(max_arith["dC_corrected"],
                                            abs((x["dC_app"] - x["dOpt"]) - x["dC_corrected"]))
        vals = [x["dC_corrected"] for x in rows]
        runs[r] = {
            "sha256_16": hashlib.sha256(raw).hexdigest()[:16],
            "sha_ok": hashlib.sha256(raw).hexdigest().startswith(EXPECTED_CKPT[r]),
            "ends_with_newline": raw.endswith(b"\n"),
            "records": len(rows), "unique_k": len(set(ks)),
            "k_complete_0_to_499": set(ks) == set(range(B_OUT)),
            "run_id": sorted({x["run_id"] for x in rows}),
            "entropy": sorted({tuple(x["entropy"]) for x in rows}),
            "entropy_ok": {tuple(x["entropy"]) for x in rows} == {(BASE, r)},
            "B_in_attempted": sorted({x["B_in_attempted"] for x in rows}),
            "status": {s: sum(x["status"] == s for x in rows) for s in {x["status"] for x in rows}},
            "paired_complete_min": min(x["B_in_paired_complete"] for x in rows),
            "inner_failures": sum(len(x["failures"]) for x in rows),
            "H1_unreliable": sum(x["status"] != "ok" for x in rows),
            "H6_pair_count_mismatch": h6_bad,
            "events_outer": sorted({x["events_outer"] for x in rows}),
            "censored_outer": sorted({x["censored_outer"] for x in rows}),
            "dC_corrected_all_finite": all(math.isfinite(v) for v in vals),
            "fields": sorted(keyset),
        }
    out["checkpoints"] = runs
    out["record_arithmetic_max_abs_error"] = max_arith

    # --- H5 evidence ---
    h5 = {n: {"expected": e, "observed": sha256(H5 / n)[:16], "ok": sha256(H5 / n).startswith(e)}
          for n, e in H5_FILES.items()}
    comp = json.loads((H5 / "h5_comparison.json").read_text())
    run0 = {json.loads(l)["k"]: json.loads(l) for l in (CK / "run0_outer.jsonl").read_text().splitlines() if l.strip()}
    rec = [json.loads(l) for l in (H5 / "recomputed_records.jsonl").read_text().splitlines() if l.strip()]
    h5["verdict"] = comp["H5_verdict"]
    h5["n_bit_identical"] = comp["n_bit_identical"]
    h5["recomputed_hash_in_comparison_ok"] = sha256(H5 / "recomputed_records.jsonl") == comp["recomputed_records_sha256"]
    h5["script_hash_in_comparison_ok"] = sha256(H5 / "h5_spotcheck.py") == comp["script_sha256"]
    h5["independent_recheck_bit_identical"] = sum(
        float(x["dC_corrected"]).hex() == float(run0[x["k"]]["dC_corrected"]).hex() for x in rec)
    h5["independent_recheck_full_record"] = sum(x == run0[x["k"]] for x in rec)
    h5["selected_k"] = sorted(x["k"] for x in rec)
    out["h5_evidence"] = h5

    ok = (out["driver_check_frozen"] == "PASS"
          and all(v["ok"] for v in out["files_17"].values())
          and all(v["sha_ok"] and v["k_complete_0_to_499"] and v["records"] == 500 and v["unique_k"] == 500
                  and v["entropy_ok"] and v["H1_unreliable"] == 0 and v["H6_pair_count_mismatch"] == 0
                  and v["ends_with_newline"] and v["dC_corrected_all_finite"] for v in runs.values())
          and all(h5[n]["ok"] for n in H5_FILES)
          and h5["verdict"] == "PASS" and h5["independent_recheck_bit_identical"] == 20
          and h5["independent_recheck_full_record"] == 20 and h5["selected_k"] == list(range(0, 500, 25)))
    out["all_ok"] = bool(ok)
    (FP / "results" / "s0_starting_state.json").write_text(json.dumps(out, indent=1, default=list))
    print(json.dumps({"all_ok": out["all_ok"], "driver_check_frozen": out["driver_check_frozen"],
                      "files_ok": sum(v["ok"] for v in out["files_17"].values()),
                      "runs": {r: {k: v[k] for k in ("sha256_16", "records", "k_complete_0_to_499", "status",
                                                     "paired_complete_min", "inner_failures",
                                                     "H6_pair_count_mismatch", "events_outer", "censored_outer")}
                               for r, v in runs.items()},
                      "arith": max_arith, "h5": {k: v for k, v in h5.items() if k != "selected_k"}},
                     indent=1, default=list))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
