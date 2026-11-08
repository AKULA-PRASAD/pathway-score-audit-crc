#!/usr/bin/env python3
"""Phase 9 corrected campaign -- independent checkpoint-level audit (read-only).

Recomputes every reported interval quantity directly from the raw checkpoint bytes, without
importing src.phase9_corrected or scripts/46. The calculation path is independent:
  * checkpoint parsing, record hashing and prefix checks use only json/hashlib;
  * DeltaC_corrected(k) is rebuilt from the per-inner-replicate C_fit / C_eval values and
    the apparent fits, not taken from the stored aggregate;
  * quantiles use a hand-written Hyndman-Fan type-7 rule (not numpy.quantile);
  * the normal CDF / inverse CDF come from statistics.NormalDist (not scipy);
  * the kernel quantile SE follows plan v2.1 section 2.2 with Silverman (1986) eq. 3.31.

Usage (repository root):
    python3 scripts/48_phase9_independent_audit.py [output_dir]
Default output_dir: reports/phase9/independent_audit. Files are created exclusively; an
existing file is never replaced. Nothing else is written.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from statistics import NormalDist


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "reports/phase9/corrected_rerun/runs"
HALTS = ROOT / "reports/phase9/corrected_rerun/halts"
OFFICIAL = ROOT / "reports/phase9/corrected_rerun_report/corrected_interval_results.json"
VERIFY = ROOT / "reports/phase9/corrected_rerun/verification/L00750.report.json"
EXECUTION = ROOT / "reports/phase9/execution/phase9_execution_results.json"
CONFIG = ROOT / "config/analysis/PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json"
RUN_IDS = (0, 1, 2, 3, 4)
B_OUT, B_IN, ALPHA = 750, 100, 0.05
N = NormalDist()


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def canonical(obj) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=True).encode()


def q7(x, p):
    """Hyndman-Fan type 7: h = (n-1)p, linear interpolation between order statistics."""
    s = sorted(float(v) for v in x)
    h = (len(s) - 1) * p
    lo = math.floor(h)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (h - lo) * (s[hi] - s[lo])


def sd1(x):
    m = sum(x) / len(x)
    return math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1))


def silverman(x):
    return 0.9 * min(sd1(x), (q7(x, 0.75) - q7(x, 0.25)) / 1.34) * len(x) ** -0.2


def kde(x, t, h):
    return sum(math.exp(-0.5 * ((t - v) / h) ** 2) for v in x) / (len(x) * h * math.sqrt(2 * math.pi))


def kernel_se(x, p, h):
    return math.sqrt(p * (1 - p) / len(x)) / kde(x, q7(x, p), h)


def rebuild(rec):
    """DeltaC_corrected from first principles (plan v2.1 section 3.1)."""
    capp = {m: rec["apparent"][m]["C_fit"] for m in ("M1", "M2")}
    opt = {"M1": [], "M2": []}
    n_pairs = 0
    for row in rec["inner"]:
        f = row["fits"]
        ok = all(f[m]["ok"] for m in ("M1", "M2"))
        if ok != row["paired_complete"]:
            raise AssertionError(f"paired_complete flag inconsistent at k={rec['k']} b={row['b']}")
        if ok:
            n_pairs += 1
            for m in ("M1", "M2"):
                opt[m].append(f[m]["C_fit"] - f[m]["C_eval"])
    om = {m: sum(opt[m]) / len(opt[m]) for m in opt}
    dapp = capp["M2"] - capp["M1"]
    dopt = om["M2"] - om["M1"]
    return {"dC_app": dapp, "dOpt": dopt, "dC_corrected": dapp - dopt, "n_pairs": n_pairs,
            "n_inner": len(rec["inner"]), "C_corr": {m: capp[m] - om[m] for m in om}}


def audit_run(r, cfg_sha, man_sha):
    raw = (RUNS / f"run{r}.jsonl").read_bytes()
    lines = raw.split(b"\n")
    assert lines[-1] == b"", "checkpoint does not end with a newline"
    lines = lines[:-1]
    head = json.loads(lines[0])
    e1 = json.loads((HALTS / f"run{r}.inv0001.end.json").read_bytes())
    s2 = json.loads((HALTS / f"run{r}.inv0002.start.json").read_bytes())
    e2 = json.loads((HALTS / f"run{r}.inv0002.end.json").read_bytes())
    prefix = b"\n".join(lines[:501]) + b"\n"
    out = {"file_bytes": len(raw), "file_sha256": sha(raw),
           "header_ok": head.get("kind") == "header" and head.get("run_id") == r
           and head.get("config_sha256") == cfg_sha and head.get("manifest_sha256") == man_sha,
           "end_record_sha_match": sha(raw) == e2["checkpoint_sha256"],
           "end_record_n_records": e2.get("n_records"), "end_record_outcome": e2.get("outcome"),
           "prefix500_sha256": sha(prefix), "prefix500_equals_inv1_end": sha(prefix) == e1["checkpoint_sha256"],
           "topup_started_from_inv1": s2["checkpoint_sha256"] == e1["checkpoint_sha256"]
           and s2["path"] == "topup" and int(s2["n_start"]) == 500 and int(s2["n_new"]) == 250}
    vals, recs, bad_hash, stored_diff, inner_counts, pair_counts = [], [], 0, 0.0, set(), set()
    for i, ln in enumerate(lines[1:]):
        o = json.loads(ln)
        rec = o["record"]
        assert o["kind"] == "record" and o["run_id"] == r and o["k"] == i == rec["k"], f"k order at {i}"
        if sha(canonical(rec)) != o["record_sha256"]:
            bad_hash += 1
        assert rec["status"] == "ok" and not rec["failures"], f"non-ok record k={i}"
        x = rebuild(rec)
        inner_counts.add(x["n_inner"]); pair_counts.add(x["n_pairs"])
        stored_diff = max(stored_diff, abs(x["dC_corrected"] - rec["dC_corrected"]),
                          abs(x["dC_app"] - rec["dC_app"]), abs(x["dOpt"] - rec["dOpt"]))
        vals.append(x["dC_corrected"]); recs.append(x)
    out.update({"n_records": len(vals), "record_hash_mismatches": bad_hash,
                "inner_per_outer": sorted(inner_counts), "paired_complete_per_outer": sorted(pair_counts),
                "max_abs_diff_rebuilt_vs_stored": stored_diff})
    return out, vals, recs


def interval_stats(v, theta):
    lo, hi = q7(v, ALPHA / 2), q7(v, 1 - ALPHA / 2)
    w = hi - lo
    h = silverman(v)
    se = {e: {str(m): kernel_se(v, p, m * h) for m in (0.5, 1.0, 2.0)}
          for e, p in (("lower", ALPHA / 2), ("upper", 1 - ALPHA / 2))}
    n_less = sum(1 for t in v if t < theta)
    n_tie = sum(1 for t in v if t == theta)
    z0 = N.inv_cdf(n_less / len(v))
    a1, a2 = N.cdf(2 * z0 + N.inv_cdf(ALPHA / 2)), N.cdf(2 * z0 + N.inv_cdf(1 - ALPHA / 2))
    mid, med = 0.5 * (lo + hi), q7(v, 0.5)
    return {"percentile": [lo, hi], "width": w, "basic": [2 * theta - hi, 2 * theta - lo],
            "z0": z0, "n_less": n_less, "n_tie": n_tie, "bc": [q7(v, a1), q7(v, a2)],
            "A7": 2 * abs(mid - theta) / w, "A8": abs(z0),
            "A7_tail": mid - med, "A7_median_offset": med - theta,
            "h": h, "kernel_SE": se,
            "ep1_ratio": {e: max(se[e].values()) / w for e in se},
            "mean": sum(v) / len(v), "median": med, "sd": sd1(v), "min": min(v), "max": max(v),
            "share_below_zero": sum(1 for t in v if t < 0) / len(v)}


def main():
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "reports/phase9/independent_audit"
    cfg = json.loads(CONFIG.read_bytes())
    official = json.loads(OFFICIAL.read_bytes())
    ver = json.loads(VERIFY.read_bytes())
    ex = json.loads(EXECUTION.read_bytes())
    cfg_sha, man_sha = ver["config_sha256"], ver["manifest_sha256"]
    mc = ex["monte_carlo"]
    theta = mc["M2_primary"]["corrected_C_mean"] - mc["M1"]["corrected_C_mean"]
    dapp_D = ex["apparent"]["M2_primary"]["C"] - ex["apparent"]["M1"]["C"]
    dopt_D = (sum(s["mean_optimism_C"] for s in mc["M2_primary"]["per_seed"])
              - sum(s["mean_optimism_C"] for s in mc["M1"]["per_seed"])) / 5

    integrity, stats, allv, comp = {}, {}, [], {}
    for r in RUN_IDS:
        integrity[r], v, recs = audit_run(r, cfg_sha, man_sha)
        stats[r] = interval_stats(v, theta)
        comp[r] = {"mean_dC_app": sum(x["dC_app"] for x in recs) / len(recs),
                   "mean_dOpt": sum(x["dOpt"] for x in recs) / len(recs),
                   "mean_C_corrected": {m: sum(x["C_corr"][m] for x in recs) / len(recs) for m in ("M1", "M2")}}
        allv += v

    between = {e: {"values": [stats[r]["percentile"][i] for r in RUN_IDS]} for i, e in enumerate(("lower", "upper"))}
    for e in between:
        vv = between[e]["values"]
        between[e].update({"sd_ddof1": sd1(vv), "range": max(vv) - min(vv)})
    pooled = {"n": len(allv), "percentile": [q7(allv, 0.025), q7(allv, 0.975)],
              "mean": sum(allv) / len(allv), "median": q7(allv, 0.5),
              "share_below_zero": sum(1 for t in allv if t < 0) / len(allv)}
    rho_max = max(max(stats[r]["ep1_ratio"].values()) for r in RUN_IDS)

    # ---- reconciliation with the official results file
    diffs = []
    def cmp(name, mine, theirs):
        d = abs(mine - theirs)
        diffs.append({"quantity": name, "independent": mine, "official": theirs, "abs_diff": d})
    cmp("theta_hat", theta, official["theta_hat"])
    for r in RUN_IDS:
        s, o = stats[r], official["runs"][str(r)]
        cmp(f"run{r}.percentile.lower", s["percentile"][0], o["percentile_primary"]["lower"])
        cmp(f"run{r}.percentile.upper", s["percentile"][1], o["percentile_primary"]["upper"])
        cmp(f"run{r}.width", s["width"], o["percentile_primary"]["width"])
        cmp(f"run{r}.basic.lower", s["basic"][0], o["basic_diagnostic"]["lower"])
        cmp(f"run{r}.basic.upper", s["basic"][1], o["basic_diagnostic"]["upper"])
        cmp(f"run{r}.bc.lower", s["bc"][0], o["bc_z0_diagnostic"]["lower"])
        cmp(f"run{r}.bc.upper", s["bc"][1], o["bc_z0_diagnostic"]["upper"])
        cmp(f"run{r}.z0", s["z0"], o["z0"])
        cmp(f"run{r}.n_less_than_theta", s["n_less"], o["z0_counts"]["n_less"])
        cmp(f"run{r}.A7", s["A7"], o["A7"])
        cmp(f"run{r}.A8", s["A8"], o["A8_abs_z0"])
        cmp(f"run{r}.silverman_h", s["h"], o["kernel_mcse"]["h"])
        for e in ("lower", "upper"):
            for m in ("0.5", "1.0", "2.0"):
                cmp(f"run{r}.kernel_SE.{e}.x{m}", s["kernel_SE"][e][m], o["kernel_mcse"][e]["SE"][m])
        cmp(f"run{r}.mean", s["mean"], o["distribution_dC_corrected"]["mean"])
        cmp(f"run{r}.median", s["median"], o["distribution_dC_corrected"]["median"])
        cmp(f"run{r}.sd", s["sd"], o["distribution_dC_corrected"]["sd"])
        cmp(f"run{r}.share_below_zero", s["share_below_zero"], o["distribution_dC_corrected"]["share_below_zero"])
        cmp(f"run{r}.mean_dC_app", comp[r]["mean_dC_app"], o["mean_dC_app"])
        cmp(f"run{r}.mean_dOpt", comp[r]["mean_dOpt"], o["mean_dOpt"])
        for m in ("M1", "M2"):
            cmp(f"run{r}.mean_C_corrected.{m}", comp[r]["mean_C_corrected"][m], o["mean_C_corrected"][m])
    for e in ("lower", "upper"):
        cmp(f"between_run.{e}.sd", between[e]["sd_ddof1"], official["between_run"][e]["sd_ddof1"])
        cmp(f"between_run.{e}.range", between[e]["range"], official["between_run"][e]["range"])
    cmp("pooled.percentile.lower", pooled["percentile"][0], official["pooled_descriptive"]["percentile"][0])
    cmp("pooled.percentile.upper", pooled["percentile"][1], official["pooled_descriptive"]["percentile"][1])
    cmp("ep1.rho_max", rho_max, official["sizing_decision_at_level"]["ep1"]["rho_max"])
    max_diff = max(d["abs_diff"] for d in diffs)

    checks = {
        "five_runs_present": sorted(integrity) == list(RUN_IDS),
        "750_records_each": all(integrity[r]["n_records"] == B_OUT for r in RUN_IDS),
        "100_inner_pairs_each": all(integrity[r]["inner_per_outer"] == [B_IN]
                                    and integrity[r]["paired_complete_per_outer"] == [B_IN] for r in RUN_IDS),
        "total_outer": sum(integrity[r]["n_records"] for r in RUN_IDS),
        "total_inner_pairs": sum(integrity[r]["n_records"] * B_IN for r in RUN_IDS),
        "headers_bound_to_config_and_manifest": all(integrity[r]["header_ok"] for r in RUN_IDS),
        "record_hashes_all_valid": all(integrity[r]["record_hash_mismatches"] == 0 for r in RUN_IDS),
        "checkpoint_sha_equals_end_record": all(integrity[r]["end_record_sha_match"] for r in RUN_IDS),
        "prefix500_byte_identical": all(integrity[r]["prefix500_equals_inv1_end"] for r in RUN_IDS),
        "topup_started_from_500_checkpoint": all(integrity[r]["topup_started_from_inv1"] for r in RUN_IDS),
        "rebuilt_dC_matches_stored_max_abs": max(integrity[r]["max_abs_diff_rebuilt_vs_stored"] for r in RUN_IDS),
        "theta_identity_dCapp_minus_dOpt_abs_err": abs((dapp_D - dopt_D) - theta),
        "verification_verdict": ver["verdict"], "verification_reporting_permitted": ver["reporting_permitted"],
        "verification_counts": ver["status_counts"],
        "official_vs_independent_max_abs_diff": max_diff,
        "ep1_pass_rho_max_le_0.05": rho_max <= 0.05,
    }
    tol = 1e-12
    reproduced = (all(checks[k] is True for k in checks if k.endswith(("_present", "_each", "_valid", "_record",
                                                                        "_identical", "_checkpoint", "_manifest")))
                  and checks["total_outer"] == 3750 and checks["total_inner_pairs"] == 375000
                  and checks["rebuilt_dC_matches_stored_max_abs"] <= tol and max_diff <= tol
                  and ver["verdict"] == "PASS")
    classification = "INDEPENDENTLY REPRODUCED" if reproduced else (
        "REPRODUCED WITH MINOR NUMERICAL DIFFERENCES" if max_diff <= 1e-9 else "NOT REPRODUCED — MATERIAL DISCREPANCY")
    res = {"generated_utc": datetime.now(timezone.utc).isoformat(), "classification": classification,
           "tolerance_abs": tol, "theta_hat": theta,
           "theta_components": {"dC_app_D": dapp_D, "dOpt_D_2500": dopt_D},
           "checks": checks, "integrity": {str(r): integrity[r] for r in RUN_IDS},
           "per_run": {str(r): {**stats[r], **comp[r]} for r in RUN_IDS},
           "between_run": between, "pooled_descriptive": pooled, "ep1_rho_max": rho_max,
           "reconciliation": diffs,
           "sources": {p.relative_to(ROOT).as_posix(): sha(p.read_bytes())
                       for p in (OFFICIAL, VERIFY, EXECUTION, CONFIG, Path(__file__).resolve())}}
    out_dir.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(res, indent=1, sort_keys=True) + "\n").encode()
    fd = os.open(out_dir / "independent_audit_results.json", os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    os.write(fd, data); os.fsync(fd); os.close(fd)
    print(json.dumps({"classification": classification, "max_abs_diff": max_diff,
                      "rebuilt_vs_stored": checks["rebuilt_dC_matches_stored_max_abs"],
                      "results_sha256": sha(data)}, indent=1))
    return 0 if classification != "NOT REPRODUCED — MATERIAL DISCREPANCY" else 3


if __name__ == "__main__":
    sys.exit(main())
