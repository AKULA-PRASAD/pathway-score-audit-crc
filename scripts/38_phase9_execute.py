"""Phase 9 EXECUTION -- M1, M2 primary, M2 sensitivity. Authorised 2026-11-03.

Runs exactly the three authorised analyses under protocol v3.1. M3 is NOT fitted
and no transcriptome-wide selection is performed. No fixed-time / absolute-risk
quantity is computed (U1).
"""

from __future__ import annotations

import hashlib, json, sys, time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import paths                                              # noqa: E402
from src.phase9_execution import (AGE_COL, ConvergenceFailure, build_clinical,  # noqa: E402
                                  c_index, calibration_slope, fit_cox,
                                  linear_predictor, pipeline_fit, rcs_knots)
from src.score_standardisation import NonEstimableReplicate         # noqa: E402
from src.reproducibility import environment_report                  # noqa: E402

ANALYSIS = paths.PROCESSED_DIR / "analysis"
PANELS = paths.ROOT / "config" / "pathways"
OUT = paths.ROOT / "reports" / "phase9" / "execution"
EV, TM = "rfs_event_full", "rfs_time_full"
M1_VARS = [AGE_COL, "sex", "tnm.stage", "tumor.location", "chemotherapy.adjuvant"]
SEEDS = [20261102, 20261103, 20261104, 20261105, 20261106]
LAM_GRID = np.geomspace(1e-3, 10.0, 12)
B_DEFAULT = 500
NONEST_THRESHOLD = 0.05


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def boot_indices(e, rng):
    """Bootstrap resample, STRATIFIED BY EVENT STATUS, patient-level."""
    idx = np.arange(len(e)); out = []
    for lab in (1.0, 0.0):
        s = idx[np.asarray(e) == lab]
        out.append(rng.choice(s, size=len(s), replace=True))
    return np.sort(np.concatenate(out))


def risk_set_comparisons(t, e, strata=None):
    """N = sum over events of (|risk set| - 1). Model-free (D8D-4-R1)."""
    t = np.asarray(t, float); e = np.asarray(e, float)
    s = np.zeros(len(t)) if strata is None else np.asarray(strata)
    n = 0; per = {}
    for i in np.where(e == 1)[0]:
        m = int(((t >= t[i]) & (s == s[i])).sum()) - 1
        n += m
        per.setdefault(str(s[i]), {"events": 0, "comparisons": 0,
                                   "events_with_comparator": 0})
        per[str(s[i])]["events"] += 1
        per[str(s[i])]["comparisons"] += m
        per[str(s[i])]["events_with_comparator"] += int(m > 0)
    return n, per


def main() -> int:
    t_start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    B = int(sys.argv[1]) if len(sys.argv) > 1 else B_DEFAULT
    print("=" * 78); print(f"PHASE 9 EXECUTION -- B={B}, seeds={SEEDS}"); print("=" * 78)

    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    cov = pd.read_csv(ANALYSIS / "covariates.csv", keep_default_na=False,
                      dtype=str).set_index("sample_id")
    ids = [s for s in s3 if all(cov.loc[s, v] != "" for v in M1_VARS)]
    ep = pd.read_csv(ANALYSIS / "endpoint.csv").set_index("sample_id").reindex(ids)
    t = ep[TM].to_numpy(float); e = ep[EV].to_numpy(float)
    X = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")[ids].astype(np.float64)
    part = pd.read_csv(ANALYSIS / "partition_and_batch.csv").set_index("sample_id")
    batch = part.loc[ids, "combat_batch"].to_numpy()
    panels = {"M2_primary": json.loads((PANELS / "m2_panel_4set.json").read_text()),
              "M2_sensitivity": json.loads((PANELS / "m2_panel_5set.json").read_text())}
    MODELS = {"M1": None, "M2_primary": panels["M2_primary"],
              "M2_sensitivity": panels["M2_sensitivity"]}
    print(f"  n={len(ids)}  events={int(e.sum())}  censored={int((e==0).sum())}")

    # ---------------- apparent fits ----------------
    apparent, fitted = {}, {}
    for name, panel in MODELS.items():
        rng = np.random.default_rng(SEEDS[0])
        m, d_tr, _, dg = pipeline_fit(X, cov, panel, ids, ids, t, e, LAM_GRID, rng)
        lp = linear_predictor(m, d_tr)
        apparent[name] = {
            "C": c_index(t, e, lp), "calibration_slope": calibration_slope(t, e, lp),
            "n_parameters_df": int(d_tr.shape[1]), "lambda": dg["lambda"],
            "age_knots": dg["age_knots"], "age_linear_fallback": dg["age_linear_fallback"],
            "coefficients": {k: float(v) for k, v in m.params_.items()},
            "se": {k: float(v) for k, v in m.standard_errors_.items()},
            "hazard_ratios": {k: float(np.exp(v)) for k, v in m.params_.items()},
            "hr_ci_lower": {k: float(np.exp(v)) for k, v in
                            m.confidence_intervals_.iloc[:, 0].items()},
            "hr_ci_upper": {k: float(np.exp(v)) for k, v in
                            m.confidence_intervals_.iloc[:, 1].items()},
            "log_likelihood": float(m.log_likelihood_),
            "score_diag": dg.get("score_standardisation"),
            "coverage": dg.get("score_coverage"),
        }
        fitted[name] = (m, d_tr)
        print(f"  apparent {name:16s} C={apparent[name]['C']:.4f} "
              f"slope={apparent[name]['calibration_slope']:.3f} "
              f"df={d_tr.shape[1]} lambda={dg['lambda']:.4g}")

    # ---------------- bootstrap optimism, per seed ----------------
    per_seed, nonest = {m: [] for m in MODELS}, []
    traces = {m: {} for m in MODELS}
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        opt = {m: [] for m in MODELS}; opt_sl = {m: [] for m in MODELS}
        orig_C = {m: [] for m in MODELS}; dC_orig = []
        for b in range(B):
            bi = boot_indices(e, rng)
            bids = [ids[i] for i in bi]
            tb, eb = t[bi], e[bi]
            rep = {}
            for name, panel in MODELS.items():
                try:
                    m, d_tr, d_or, _ = pipeline_fit(X, cov, panel, bids, ids,
                                                    tb, eb, LAM_GRID, rng)
                    lp_b = linear_predictor(m, d_tr)
                    lp_o = linear_predictor(m, d_or)
                    cb, co = c_index(tb, eb, lp_b), c_index(t, e, lp_o)
                    opt[name].append(cb - co); orig_C[name].append(co)
                    opt_sl[name].append(calibration_slope(tb, eb, lp_b)
                                        - calibration_slope(t, e, lp_o))
                    rep[name] = co
                except (ConvergenceFailure, NonEstimableReplicate) as exc:
                    nonest.append({"model": name, "seed": seed, "replicate": b,
                                   "reason": getattr(exc, "reason", str(exc))})
            if {"M1", "M2_primary"} <= rep.keys():
                dC_orig.append(rep["M2_primary"] - rep["M1"])
        for name in MODELS:
            o = np.array(opt[name], float)
            sl = np.array(opt_sl[name], float)
            corr = apparent[name]["C"] - np.nanmean(o)
            per_seed[name].append({
                "seed": seed, "B_effective": int(len(o)),
                "mean_optimism_C": float(np.nanmean(o)),
                "corrected_C": float(corr),
                "mean_optimism_slope": float(np.nanmean(sl)),
                "corrected_slope": float(apparent[name]["calibration_slope"]
                                         - np.nanmean(sl)),
                "orig_C_p2.5": float(np.nanpercentile(orig_C[name], 2.5)),
                "orig_C_p97.5": float(np.nanpercentile(orig_C[name], 97.5)),
            })
            # MC-1: running trace of the corrected estimate
            run = [float(apparent[name]["C"] - np.nanmean(o[:k]))
                   for k in (50, 100, 200, 300, 400, len(o)) if k <= len(o)]
            traces[name][str(seed)] = run
        per_seed.setdefault("_dC", []).append({
            "seed": seed, "dC_orig_mean": float(np.nanmean(dC_orig)),
            "dC_orig_p2.5": float(np.nanpercentile(dC_orig, 2.5)),
            "dC_orig_p97.5": float(np.nanpercentile(dC_orig, 97.5))})
        print(f"  seed {seed}: " + "  ".join(
            f"{m} C_corr={per_seed[m][-1]['corrected_C']:.4f}" for m in MODELS))

    # ---------------- MC-2/3/4 ----------------
    mc = {}
    for name in MODELS:
        cs = np.array([r["corrected_C"] for r in per_seed[name]])
        samp_sd = float(np.nanstd(
            [r["orig_C_p97.5"] - r["orig_C_p2.5"] for r in per_seed[name]]))
        boot_sd = float((np.mean([r["orig_C_p97.5"] for r in per_seed[name]])
                         - np.mean([r["orig_C_p2.5"] for r in per_seed[name]])) / 3.92)
        mc_se = float(cs.std(ddof=1))
        mc[name] = {"corrected_C_mean": float(cs.mean()),
                    "monte_carlo_SE": mc_se, "sampling_SE_approx": boot_sd,
                    "MC3_ratio": float(mc_se / boot_sd) if boot_sd else float("inf"),
                    "MC3_pass": bool(boot_sd and mc_se <= 0.10 * boot_sd),
                    "per_seed": per_seed[name], "running_trace_MC1": traces[name]}
        print(f"  MC {name:16s} C_corr={cs.mean():.4f} MC-SE={mc_se:.5f} "
              f"samplingSE~{boot_sd:.5f} ratio={mc[name]['MC3_ratio']:.3f} "
              f"MC-3 {'PASS' if mc[name]['MC3_pass'] else 'FAIL'}")

    n_total = B * len(SEEDS)
    nonest_rate = {m: sum(1 for x in nonest if x["model"] == m) / n_total for m in MODELS}

    # ---------------- batch sensitivity (D8D-4-R1) ----------------
    N_un, _ = risk_set_comparisons(t, e)
    N_st, per_stratum = risk_set_comparisons(t, e, batch)
    for b_ in np.unique(batch):
        per_stratum.setdefault(str(b_), {"events": 0, "comparisons": 0,
                                         "events_with_comparator": 0})
        per_stratum[str(b_)]["n"] = int((batch == b_).sum())
    rng = np.random.default_rng(SEEDS[0])
    strat_res, lobo = {}, []
    try:
        m, d_tr, _, dg = pipeline_fit(X, cov, panels["M2_primary"], ids, ids,
                                      t, e, LAM_GRID, rng, strata_train=batch)
        H = m.variance_matrix_.to_numpy(float)
        ev_ = np.linalg.eigvalsh(H)
        strat_res = {"converged": True,
                     "coefficients": {k: float(v) for k, v in m.params_.items()},
                     "se": {k: float(v) for k, v in m.standard_errors_.items()},
                     "covariance_min_eigenvalue": float(ev_.min()),
                     "covariance_positive_definite": bool(ev_.min() > 0),
                     "max_abs_beta": float(np.abs(m.params_).max()),
                     "max_se": float(m.standard_errors_.max()),
                     "separation_sentinel_tripped": False}
    except (ConvergenceFailure, NonEstimableReplicate) as exc:
        strat_res = {"converged": False, "reason": getattr(exc, "reason", str(exc))}
    for b_ in sorted(np.unique(batch)):
        keep = batch != b_
        kids = [i for i, k in zip(ids, keep) if k]
        try:
            mm, dd, _, _ = pipeline_fit(X, cov, panels["M2_primary"], kids, kids,
                                        t[keep], e[keep], LAM_GRID,
                                        np.random.default_rng(SEEDS[0]))
            lobo.append({"omitted_batch": str(b_), "n_omitted": int((~keep).sum()),
                         "n_fit": len(kids),
                         "events_fit": int(e[keep].sum()),
                         "C_apparent": c_index(t[keep], e[keep], linear_predictor(mm, dd)),
                         "coefficients": {k: float(v) for k, v in mm.params_.items()}})
        except Exception as exc:
            lobo.append({"omitted_batch": str(b_), "failed": str(exc)[:120]})

    out = {
        "phase": "9", "stage": "EXECUTION", "protocol_version": "3.1",
        "executed_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": round(time.time() - t_start, 1),
        "authorised_analyses": list(MODELS), "M3_fitted": False,
        "endpoint": {"event_col": EV, "time_col": TM,
                     "basis": "Register A-07 cause-specific recurrence, full follow-up"},
        "population": {"n": len(ids), "events": int(e.sum()),
                       "censored": int((e == 0).sum())},
        "B": B, "seeds": SEEDS, "lambda_grid": [float(x) for x in LAM_GRID],
        "apparent": apparent, "monte_carlo": mc,
        "delta_C": per_seed["_dC"],
        "non_estimable_replicates": {
            "records": nonest, "n_total_replicates_per_model": n_total,
            "rate": nonest_rate,
            "threshold": NONEST_THRESHOLD,
            "any_model_unreliable": any(v > NONEST_THRESHOLD for v in nonest_rate.values()),
            "redrawn": False},
        "batch_sensitivity": {
            "strata": "combat_batch, nine individual strata, NO pooling",
            "stratified_fit": strat_res,
            "leave_one_batch_out": lobo,
            "risk_set_comparisons": {
                "N_unstratified": int(N_un), "N_stratified": int(N_st),
                "retention": float(N_st / N_un) if N_un else None,
                "per_stratum": per_stratum,
                "interpretation": "descriptive only; counts COMPARISONS, not Fisher "
                                  "information; no threshold applied"}},
        "prohibited_quantities_not_computed": [
            "t-year / fixed-horizon absolute risk", "CIF at named years",
            "time-dependent AUC", "Brier at a fixed horizon",
            "fixed-horizon cloglog calibration"],
        "input_hashes": {str(Path(p).relative_to(paths.ROOT)): sha(p) for p in (
            ANALYSIS / "expression_median_S0.parquet", ANALYSIS / "covariates.csv",
            ANALYSIS / "endpoint.csv", ANALYSIS / "cohort_membership.csv",
            ANALYSIS / "partition_and_batch.csv",
            PANELS / "m2_panel_4set.json", PANELS / "m2_panel_5set.json",
            PANELS / "M2_FREEZE_RECORD.json",
            paths.ROOT / "src" / "pathway_scoring.py",
            paths.ROOT / "src" / "score_standardisation.py",
            paths.ROOT / "src" / "phase9_execution.py")},
        "environment": environment_report(),
    }
    (OUT / "phase9_execution_results.json").write_text(
        json.dumps(out, indent=2, default=str) + "\n")
    print(f"\n  runtime {out['runtime_seconds']}s -> "
          f"{(OUT/'phase9_execution_results.json').relative_to(paths.ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
