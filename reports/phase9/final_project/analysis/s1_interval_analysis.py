"""Stage B/C/D -- plan v2.1 section 9.1 steps 5-6 (except the H5 spot-check, already done):

  per run:   percentile (primary), basic and z0 / BC(z0) (diagnostics); kernel endpoint MC-SE at
             Silverman x {0.5, 1, 2}; non-estimable accounting (6.2, 6.3, 6.5); bounding (6.4);
             MC-1 running trace (W5); inner-noise inflation (1.2, W3); A7 (W6); A8 (W7)
  across:    between-run endpoint SD (ddof = 1, df = 4); EP-1; W4; listed endpoints (2.1 #5)

Inputs are opened read-only: the five checkpoints and the frozen phase9_execution_results.json
(for theta_hat). Outputs go ONLY to reports/phase9/final_project/{results,tables}/.

Every quantity derived from dC_corrected depends on M2's lambda tuning and therefore carries the
label PROVISIONAL (CV-seed protocol deviation, not yet adjudicated).

Usage: .venv/bin/python -B -u reports/phase9/final_project/analysis/s1_interval_analysis.py
"""
from __future__ import annotations

import csv, json, math, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stats_core as sc                                   # noqa: E402
from common import (FP, N_RUNS, B_OUT, B_IN, EXPECTED_CKPT, CKPT_DIR, load_run,  # noqa: E402
                    sha256, theta_hat)

PROV = ("PROVISIONAL -- calculated using the as-run implementation; "
        "CV-seed protocol deviation not yet adjudicated.")
BW_MULT = (0.5, 1.0, 2.0)
EP1_FRACTION = 0.05            # plan section 2: MC-SE(endpoint) <= 0.05 x (upper - lower)
A7_FRACTION, A8_LIMIT, H1_LIMIT = 0.25, 0.25, 0.05
W4_FACTOR, W3_LIMIT = 2.0, 0.05


def summarise(v):
    v = np.asarray(v, dtype=np.float64)
    return {"n": int(v.size), "mean": float(v.mean()), "sd": float(v.std(ddof=1)),
            "median": sc.quantile(v, .5), "min": float(v.min()), "max": float(v.max()),
            "q025": sc.quantile(v, .025), "q975": sc.quantile(v, .975),
            "skewness": float(sps.skew(v, bias=False)), "excess_kurtosis": float(sps.kurtosis(v, bias=False)),
            "n_negative": int((v < 0).sum()), "frac_negative": float((v < 0).mean())}


def per_run(r, theta):
    recs = load_run(r)
    assert [x["k"] for x in recs] == list(range(B_OUT))
    # 6.5 -- unreliable outer replicates (apparent failure / paired < 95 % / bounding width > MC-SE)
    bounds = [sc.bounding_record(x, B_IN) for x in recs]
    unrel_app = [x["k"] for x in recs if x["status"] == "unreliable_apparent_failed"]
    unrel_pair = [x["k"] for x in recs if x.get("B_in_paired_complete", 0) < 0.95 * B_IN]
    usable = [x for x in recs if x["status"] == "ok"]
    v = np.array([x["dC_corrected"] for x in usable], dtype=np.float64)

    pct = sc.percentile_ci(v); bas = sc.basic_ci(v, theta); z = sc.z0(v, theta); bc = sc.bc_z0_ci(v, theta)
    width = pct[1] - pct[0]

    h1 = sc.silverman_bandwidth(v); hs = sc.scipy_silverman_bandwidth(v)
    kern = {}
    for name, p, q in (("lower", 0.025, pct[0]), ("upper", 0.975, pct[1])):
        kern[name] = {f"bw_x{m}": sc.kernel_quantile_se(v, p, q, h1 * m) for m in BW_MULT}
        kern[name]["crosscheck_scipy_silverman_x1"] = sc.kernel_quantile_se(v, p, q, hs)
        kern[name]["p"] = p; kern[name]["q"] = q
    kern_bc = {}
    if not bc["boundary"]:
        for name, p, q in (("lower", bc["alpha1"], bc["lower"]), ("upper", bc["alpha2"], bc["upper"])):
            kern_bc[name] = sc.kernel_quantile_se(v, p, q, h1)

    # 6.4 bounding: interval under min- and max-imputation
    bmin = np.array([b["dC_min"] for b, x in zip(bounds, recs) if x["status"] == "ok"])
    bmax = np.array([b["dC_max"] for b, x in zip(bounds, recs) if x["status"] == "ok"])
    pmin, pmax = sc.percentile_ci(bmin), sc.percentile_ci(bmax)
    move = max(abs(pmin[0] - pct[0]), abs(pmax[0] - pct[0]), abs(pmin[1] - pct[1]), abs(pmax[1] - pct[1]))
    kse_lo, kse_hi = kern["lower"]["bw_x1.0"]["SE"], kern["upper"]["bw_x1.0"]["SE"]
    unrel_bound = [x["k"] for b, x in zip(bounds, recs) if b["width"] > min(kse_lo, kse_hi)]

    trace = sc.running_trace([x["dC_corrected"] for x in recs])
    fin = [trace[B] for B in sc.MC1_FINAL]
    w5_spread = {"lower": max(t[0] for t in fin) - min(t[0] for t in fin),
                 "upper": max(t[1] for t in fin) - min(t[1] for t in fin)}
    fine = {B: sc.percentile_ci(v[:B]) for B in range(50, B_OUT + 1, 10)}

    noise = sc.inner_noise(usable)
    n_unrel = len(set(unrel_app) | set(unrel_pair) | set(unrel_bound))
    failures = sum(len(x["failures"]) for x in recs)
    a7 = {"lower": abs(pct[0] - bas[0]) / width, "upper": abs(pct[1] - bas[1]) / width}

    return {
        "run_id": r, "checkpoint_sha256": sha256(CKPT_DIR / f"run{r}_outer.jsonl"),
        "label": PROV,
        "n_records": len(recs), "n_usable": int(v.size),
        "accounting_6_2": {"B_in_attempted": sorted({x["B_in_attempted"] for x in recs}),
                           "B_in_paired_complete_min": min(x["B_in_paired_complete"] for x in recs),
                           "B_in_paired_complete_max": max(x["B_in_paired_complete"] for x in recs),
                           "inner_failures_total": failures,
                           "failures_by_model": {m: sum(f.get("model") == m for x in recs for f in x["failures"])
                                                 for m in ("M1", "M2")}},
        "representativeness_6_3": ("VACUOUS -- 0 failed inner replicates; no failed-vs-successful "
                                   "comparison is possible" if failures == 0 else "REQUIRED -- not implemented"),
        "unreliable_6_5": {"apparent_failed": unrel_app, "paired_below_95pct": unrel_pair,
                           "bounding_width_gt_MCSE": unrel_bound, "n_unreliable": n_unrel,
                           "fraction": n_unrel / len(recs), "H1_breach": n_unrel / len(recs) > H1_LIMIT},
        "bounding_6_4": {"n_outer_with_missing": sum(b["n_missing"] > 0 for b in bounds),
                         "interval_min_imputation": pmin, "interval_max_imputation": pmax,
                         "max_endpoint_movement": move,
                         "endpoint_MCSE_reference": [kse_lo, kse_hi],
                         "H2_breach": move > min(kse_lo, kse_hi),
                         "status": "VACUOUS (no missing paired values)" if all(b["n_missing"] == 0 for b in bounds)
                         else "COMPUTED"},
        "distribution": {"dC_corrected": summarise(v),
                         "dC_app": summarise([x["dC_app"] for x in usable]),
                         "dOpt": summarise([x["dOpt"] for x in usable]),
                         "opt_M1_mean": summarise([x["opt_M1_mean"] for x in usable]),
                         "opt_M2_mean": summarise([x["opt_M2_mean"] for x in usable]),
                         "C_app_M1": summarise([x["C_app"]["M1"] for x in usable]),
                         "C_app_M2": summarise([x["C_app"]["M2"] for x in usable]),
                         "n_distinct_outer": summarise([x["n_distinct_outer"] for x in usable]),
                         "lambda_apparent_M2": summarise([x["lambda_apparent"]["M2"] for x in usable]),
                         "lambda_apparent_M1_distinct": sorted({x["lambda_apparent"]["M1"] for x in usable})},
        "percentile_primary": {"lower": pct[0], "upper": pct[1], "width": width,
                               "quantile_levels": [0.025, 0.975], "quantile_method": "type 7 (numpy 'linear')"},
        "basic_diagnostic": {"lower": bas[0], "upper": bas[1], "width": bas[1] - bas[0], "theta_hat": theta},
        "z0_diagnostic": z,
        "bc_z0_diagnostic": bc,
        "kernel_MCSE_2_2": {"silverman_h": h1, "scipy_silverman_h_crosscheck": hs,
                            "bandwidth_multipliers": list(BW_MULT), "percentile_endpoints": kern,
                            "bc_endpoints_bw_x1": kern_bc},
        "MC1_running_trace": {"grid": {str(B): list(t) for B, t in trace.items()},
                              "final_40pct_grid": list(sc.MC1_FINAL), "final_40pct_spread": w5_spread,
                              "fine_trace_step10": {str(B): list(t) for B, t in fine.items()}},
        "inner_noise_1_2": noise,
        "warnings": {
            "W1_inner_failures_gt0": failures > 0,
            "W2_failure_concentration": "not applicable (0 failures)" if failures == 0 else "evaluate",
            "W3_SD_inflation": noise["SD_inflation"], "W3_breach": noise["SD_inflation"] > W3_LIMIT,
            "W5_spread": w5_spread, "W5_bound": EP1_FRACTION * width,
            "W5_breach": max(w5_spread.values()) > EP1_FRACTION * width,
            "W6_A7_fraction_of_width": a7, "W6_breach": max(a7.values()) > A7_FRACTION,
            "W7_A8_abs_z0": abs(z["z0"]), "W7_breach": abs(z["z0"]) > A8_LIMIT,
        },
    }


def main() -> int:
    t = theta_hat(); theta = t["value"]
    runs = [per_run(r, theta) for r in range(N_RUNS)]

    # ---- across runs: listed endpoints (2.1 #5) and between-run SD (2.1 #3; ddof = 1, df = 4)
    def col(get):
        return [get(x) for x in runs]
    endpoints = {
        "percentile_lower": col(lambda x: x["percentile_primary"]["lower"]),
        "percentile_upper": col(lambda x: x["percentile_primary"]["upper"]),
        "basic_lower": col(lambda x: x["basic_diagnostic"]["lower"]),
        "basic_upper": col(lambda x: x["basic_diagnostic"]["upper"]),
        "bc_z0_lower": col(lambda x: x["bc_z0_diagnostic"]["lower"]),
        "bc_z0_upper": col(lambda x: x["bc_z0_diagnostic"]["upper"]),
        "z0": col(lambda x: x["z0_diagnostic"]["z0"]),
        "percentile_width": col(lambda x: x["percentile_primary"]["width"]),
        "median_dC_corrected": col(lambda x: x["distribution"]["dC_corrected"]["median"]),
        "mean_dC_corrected": col(lambda x: x["distribution"]["dC_corrected"]["mean"]),
    }
    between = {k: {"values": v, "mean": float(np.mean(v)), "sd_ddof1": float(np.std(v, ddof=1)),
                   "range": float(max(v) - min(v)), "df": len(v) - 1} for k, v in endpoints.items()}

    # ---- EP-1 (plan section 2): MC-SE(endpoint) <= 0.05 x (upper - lower), BOTH endpoints.
    # The plan does not single out one MC-SE estimator, so EVERY estimator it prescribes is evaluated:
    # kernel at Silverman x {0.5, 1, 2} (2.2) and the between-run SD (2.1 #3).
    grid = []
    for x in runs:
        w = x["percentile_primary"]["width"]
        for e in ("lower", "upper"):
            ests = {f"kernel_bw_x{m}": x["kernel_MCSE_2_2"]["percentile_endpoints"][e][f"bw_x{m}"]["SE"]
                    for m in BW_MULT}
            ests["between_run_SD"] = between[f"percentile_{e}"]["sd_ddof1"]
            for name, se in ests.items():
                grid.append({"run": x["run_id"], "endpoint": e, "estimator": name, "MCSE": se,
                             "width": w, "ratio": se / w, "bound": EP1_FRACTION * w,
                             "pass": se <= EP1_FRACTION * w})
    n_pass = sum(g["pass"] for g in grid)
    worst = max(grid, key=lambda g: g["ratio"])
    ep1_all = n_pass == len(grid)
    by_est = {}
    for g in grid:
        s = by_est.setdefault(g["estimator"], {"n": 0, "n_pass": 0, "max_ratio": 0.0})
        s["n"] += 1; s["n_pass"] += g["pass"]; s["max_ratio"] = max(s["max_ratio"], g["ratio"])
    for s in by_est.values():
        s["all_pass"] = s["n_pass"] == s["n"]
    failing = [g for g in grid if not g["pass"]]
    # plan 2.4 ladder arithmetic, shown for the failing cases only (NOT executed)
    ladder = (750, 1000, 1500, 2000, 3000)
    req = 500 * (worst["ratio"] / EP1_FRACTION) ** 2
    rung = next((L for L in ladder if L >= req), None) if not ep1_all else None
    if ep1_all:
        status = "PROVISIONAL -- numerical PASS under every prescribed MC-SE estimator"
    elif by_est["kernel_bw_x1.0"]["all_pass"] and by_est["between_run_SD"]["all_pass"]:
        status = ("UNRESOLVED (decision required) -- PASS under the plan's named Silverman bandwidth and "
                  "under the between-run SD; FAIL in " + str(len(failing)) + "/" + str(len(grid)) +
                  " evaluations at a bandwidth-sensitivity multiplier. The plan does not state which "
                  "estimator governs EP-1.")
    else:
        status = "FAIL under the plan's named Silverman bandwidth and/or the between-run SD"
    ep1 = {"criterion": "MC-SE(endpoint) <= 0.05 x (upper - lower) for both endpoints (plan v2.1 section 2)",
           "interval": "percentile (primary), per run, B_out = 500",
           "n_evaluations": len(grid), "n_pass": n_pass, "all_pass": ep1_all,
           "by_estimator": by_est, "failing_cases": failing,
           "max_ratio": worst["ratio"], "max_ratio_case": worst,
           "numerical_outcome": "PASS" if ep1_all else ("FAIL" if n_pass == 0 else "MIXED"),
           "reported_status": status,
           "plan_2_4_ladder_if_strict": {"required_B_out": req, "ladder_rung": rung,
                                         "executed": False,
                                         "why_not": ("top-up would write new outer replicates to the run "
                                                     "checkpoints and is a resource/scientific decision; "
                                                     "not authorised in this work order")},
           "grid": grid, "label": PROV}

    # ---- W4: kernel (bw x1, per run) vs between-run SD, disagreement > 2x
    w4 = []
    for e in ("lower", "upper"):
        b = between[f"percentile_{e}"]["sd_ddof1"]
        for x in runs:
            k = x["kernel_MCSE_2_2"]["percentile_endpoints"][e]["bw_x1.0"]["SE"]
            w4.append({"run": x["run_id"], "endpoint": e, "kernel_bw_x1": k, "between_run_SD": b,
                       "ratio_kernel_over_between": k / b, "breach": max(k / b, b / k) > W4_FACTOR})

    # ---- H6 (11.1): paired-complete count must be equal for M1 and M2 in every outer record
    h6_count = sum(not (len(x_["inner_C"]["M1"]) == len(x_["inner_C"]["M2"]) == x_["B_in_paired_complete"])
                   for r in range(N_RUNS) for x_ in load_run(r))

    # ---- observed decomposition of the outer-vs-theta shift (descriptive, from stored values)
    allrecs = [x_ for r in range(N_RUNS) for x_ in load_run(r)]
    m_app = float(np.mean([x_["dC_app"] for x_ in allrecs])); m_opt = float(np.mean([x_["dOpt"] for x_ in allrecs]))
    bias = {"note": "DESCRIPTIVE -- pooled means over 2500 outer replicates vs the original-data values",
            "dC_app_D": t["dC_app_D"], "mean_dC_app_k": m_app, "shift_dC_app": m_app - t["dC_app_D"],
            "dOpt_D": t["dOpt_D"], "mean_dOpt_k": m_opt, "shift_dOpt": m_opt - t["dOpt_D"],
            "theta_hat": theta, "mean_dC_corrected_k": m_app - m_opt,
            "shift_dC_corrected": (m_app - m_opt) - theta}
    # z0 under theta_hat +/- 2 MC-SE (theta_hat is itself a 5-seed Monte-Carlo estimate)
    se_t = t["MC_SE_of_mean"]
    z0sens = {"theta_MC_SE": se_t,
              "runs": {x["run_id"]: {f"{m:+d}SE": sc.z0([y["dC_corrected"] for y in load_run(x["run_id"])],
                                                        theta + m * se_t)["z0"] for m in (-2, 0, 2)} for x in runs},
              "note": "checks whether theta_hat's own Monte-Carlo error could account for |z0| > 0.25"}

    # ---- descriptive pooled summary (NOT a plan-designated interval)
    allv = np.concatenate([[r_["dC_corrected"] for r_ in load_run(r)] for r in range(N_RUNS)])
    pooled = {"note": ("DESCRIPTIVE ONLY -- the plan does not designate a pooled interval; the 2500 "
                       "draws come from 5 independent runs of the same design"),
              "n": int(allv.size), "percentile": sc.percentile_ci(allv), "basic": sc.basic_ci(allv, theta),
              "z0": sc.z0(allv, theta), "distribution": summarise(allv)}

    out = {"analysis": "Phase 9 nested-bootstrap interval for optimism-corrected Delta-C, plan v2.1 section 9.1",
           "generated_utc": datetime.now(timezone.utc).isoformat(), "label": PROV,
           "conditional_statement": ("Conditional on the observed 136 recurrence events and 316 censored "
                                     "observations in the 452-patient A-19 population."),
           "theta_hat": t, "alpha": sc.ALPHA, "B_out": B_OUT, "B_in": B_IN,
           "checkpoint_hashes_16": {r: sha256(CKPT_DIR / f"run{r}_outer.jsonl")[:16] for r in range(N_RUNS)},
           "checkpoint_hashes_expected": EXPECTED_CKPT,
           "runs": runs, "between_run": between, "EP1": ep1, "W4": w4, "pooled_descriptive": pooled,
           "bias_decomposition_observed": bias,
           "z0_sensitivity_to_theta_MC_error": z0sens,
           "hard_failures": {
               "H1": any(x["unreliable_6_5"]["H1_breach"] for x in runs),
               "H2": any(x["bounding_6_4"]["H2_breach"] for x in runs),
               "H3": ("NOT TRIGGERED -- H3 is 'EP-1 unmet after exhausting the B_out ladder'; the ladder "
                      "has not been run. EP-1 status: " + ep1["reported_status"]),
               "H4": "evaluated separately (pre/post hash verification)",
               "H5": "PASS (reports/phase9/h5_spotcheck/h5_comparison.json)",
               "H6": h6_count > 0, "H6_records_with_unequal_pair_counts": h6_count}}
    (FP / "results" / "s1_interval_results.json").write_text(json.dumps(out, indent=1, default=float))

    # ---- tables (CSV, generated from the same objects)
    T = FP / "tables"
    with open(T / "T1_intervals_per_run.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "n", "mean", "median", "sd", "pct_lower", "pct_upper", "pct_width",
                    "basic_lower", "basic_upper", "z0", "p_less", "bc_alpha1", "bc_alpha2", "bc_lower", "bc_upper"])
        for x in runs:
            d = x["distribution"]["dC_corrected"]; p = x["percentile_primary"]; b = x["basic_diagnostic"]
            z = x["z0_diagnostic"]; c = x["bc_z0_diagnostic"]
            w.writerow([x["run_id"], x["n_usable"], d["mean"], d["median"], d["sd"], p["lower"], p["upper"],
                        p["width"], b["lower"], b["upper"], z["z0"], z["p_less"], c["alpha1"], c["alpha2"],
                        c["lower"], c["upper"]])
    with open(T / "T2_kernel_MCSE.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "endpoint", "p", "q", "silverman_h", "SE_bw_x0.5", "SE_bw_x1.0", "SE_bw_x2.0",
                    "fhat_bw_x0.5", "fhat_bw_x1.0", "fhat_bw_x2.0", "SE_scipy_silverman_crosscheck"])
        for x in runs:
            k = x["kernel_MCSE_2_2"]
            for e in ("lower", "upper"):
                pe = k["percentile_endpoints"][e]
                w.writerow([x["run_id"], e, pe["p"], pe["q"], k["silverman_h"],
                            *[pe[f"bw_x{m}"]["SE"] for m in BW_MULT], *[pe[f"bw_x{m}"]["f_hat"] for m in BW_MULT],
                            pe["crosscheck_scipy_silverman_x1"]["SE"]])
    with open(T / "T3_between_run.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["quantity", "run0", "run1", "run2", "run3", "run4", "mean", "sd_ddof1", "range"])
        for k, v in between.items():
            w.writerow([k, *v["values"], v["mean"], v["sd_ddof1"], v["range"]])
    with open(T / "T4_EP1_grid.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["run", "endpoint", "estimator", "MCSE", "width", "ratio", "bound", "pass"])
        for g in grid:
            w.writerow([g["run"], g["endpoint"], g["estimator"], g["MCSE"], g["width"], g["ratio"], g["bound"], g["pass"]])
    with open(T / "T5_warnings_hardfailures.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "inner_failures", "n_unreliable_6_5", "H1_breach", "bounding_status", "H2_breach",
                    "W3_SD_inflation", "W3_breach", "W5_spread_lower", "W5_spread_upper", "W5_bound", "W5_breach",
                    "W6_A7_lower", "W6_A7_upper", "W6_breach", "W7_abs_z0", "W7_breach"])
        for x in runs:
            a = x["warnings"]
            w.writerow([x["run_id"], x["accounting_6_2"]["inner_failures_total"], x["unreliable_6_5"]["n_unreliable"],
                        x["unreliable_6_5"]["H1_breach"], x["bounding_6_4"]["status"], x["bounding_6_4"]["H2_breach"],
                        a["W3_SD_inflation"], a["W3_breach"], a["W5_spread"]["lower"], a["W5_spread"]["upper"],
                        a["W5_bound"], a["W5_breach"], a["W6_A7_fraction_of_width"]["lower"],
                        a["W6_A7_fraction_of_width"]["upper"], a["W6_breach"], a["W7_A8_abs_z0"], a["W7_breach"]])
    with open(T / "T6_W4_kernel_vs_between.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["run", "endpoint", "kernel_bw_x1", "between_run_SD", "ratio", "breach"])
        for g in w4:
            w.writerow([g["run"], g["endpoint"], g["kernel_bw_x1"], g["between_run_SD"], g["ratio_kernel_over_between"], g["breach"]])

    print(f"theta_hat = {theta!r}")
    for x in runs:
        p = x["percentile_primary"]; b = x["basic_diagnostic"]; z = x["z0_diagnostic"]; c = x["bc_z0_diagnostic"]
        k = x["kernel_MCSE_2_2"]["percentile_endpoints"]
        print(f"run{x['run_id']}: pct [{p['lower']:.5f}, {p['upper']:.5f}] w={p['width']:.5f} | basic "
              f"[{b['lower']:.5f}, {b['upper']:.5f}] | z0={z['z0']:+.4f} (p={z['p_less']:.3f}, ties={z['n_tie']}) | "
              f"BC [{c['lower']:.5f}, {c['upper']:.5f}] | kSE lo {k['lower']['bw_x1.0']['SE']:.5f} "
              f"hi {k['upper']['bw_x1.0']['SE']:.5f}")
    for k in ("percentile_lower", "percentile_upper"):
        print(f"between-run {k}: mean {between[k]['mean']:.5f} sd {between[k]['sd_ddof1']:.6f}")
    print(f"EP-1: {ep1['n_pass']}/{ep1['n_evaluations']} pass, max ratio {ep1['max_ratio']:.4f} -> "
          f"{ep1['numerical_outcome']}")
    print("EP-1 by estimator:", {k: (v['n_pass'], v['n'], round(v['max_ratio'], 4)) for k, v in by_est.items()})
    print("EP-1 status:", ep1["reported_status"])
    print("ladder if strict:", ep1["plan_2_4_ladder_if_strict"]["required_B_out"], "->",
          ep1["plan_2_4_ladder_if_strict"]["ladder_rung"])
    print("bias decomposition:", {k: v for k, v in bias.items() if k != "note"})
    print("z0 sensitivity:", z0sens["runs"])
    print("hard failures:", out["hard_failures"])
    print("warnings:", {x["run_id"]: {k: v for k, v in x["warnings"].items() if k.endswith("breach")} for x in runs})
    print("W4 breaches:", sum(g["breach"] for g in w4))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
