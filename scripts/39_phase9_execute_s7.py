"""Phase 9 estimand S-7 -- representative-IQR representation sensitivity.

Identical to the primary execution in every respect EXCEPT the probe-to-gene
representation: `expression_repiqr_disc443_S0.parquet` (frozen Phase 2B, A-05)
replaces the median matrix. Same population, endpoint, clinical specification,
age spline, panels, M2-A2 standardisation, P6B-A1 penalty scope, nested tuning,
seeds and B.

Reuses `src/phase9_execution.py` and the resampling helpers of
`scripts/38_phase9_execute.py` unchanged, so the only difference is the matrix.

Writes to reports/phase9/execution_s7/ -- the primary outputs are never touched.
The batch-stratified / leave-one-batch-out investigation is NOT repeated: the
frozen S-7 definition (Phase 6A A-08) prescribes an aggregation-choice re-fit
only.
"""

from __future__ import annotations

import hashlib, importlib.util, json, sys, time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import paths                                              # noqa: E402
from src.phase9_execution import (AGE_COL, ConvergenceFailure, c_index,  # noqa: E402
                                  calibration_slope, linear_predictor,
                                  pipeline_fit)
from src.score_standardisation import NonEstimableReplicate         # noqa: E402
from src.reproducibility import environment_report                  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "p9primary", ROOT / "scripts" / "38_phase9_execute.py")
_p = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_p)
boot_indices = _p.boot_indices                     # same resampling, same object

ANALYSIS = paths.PROCESSED_DIR / "analysis"
PANELS = paths.ROOT / "config" / "pathways"
OUT = paths.ROOT / "reports" / "phase9" / "execution_s7"
MATRIX = ANALYSIS / "expression_repiqr_disc443_S0.parquet"
EV, TM = "rfs_event_full", "rfs_time_full"
M1_VARS = [AGE_COL, "sex", "tnm.stage", "tumor.location", "chemotherapy.adjuvant"]
SEEDS = [20261102, 20261103, 20261104, 20261105, 20261106]
LAM_GRID = np.geomspace(1e-3, 10.0, 12)


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main() -> int:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    print("=" * 78); print(f"PHASE 9 S-7 (representative-IQR) -- B={B}"); print("=" * 78)

    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    cov = pd.read_csv(ANALYSIS / "covariates.csv", keep_default_na=False,
                      dtype=str).set_index("sample_id")
    ids = [s for s in s3 if all(cov.loc[s, v] != "" for v in M1_VARS)]
    ep = pd.read_csv(ANALYSIS / "endpoint.csv").set_index("sample_id").reindex(ids)
    t = ep[TM].to_numpy(float); e = ep[EV].to_numpy(float)
    X = pd.read_parquet(MATRIX)[ids].astype(np.float64)
    panels = {"M2_primary": json.loads((PANELS / "m2_panel_4set.json").read_text()),
              "M2_sensitivity": json.loads((PANELS / "m2_panel_5set.json").read_text())}
    MODELS = {"M1": None, "M2_primary": panels["M2_primary"],
              "M2_sensitivity": panels["M2_sensitivity"]}
    print(f"  matrix: {MATRIX.name} ({sha(MATRIX)[:16]})")
    print(f"  n={len(ids)} events={int(e.sum())} censored={int((e==0).sum())}")

    apparent = {}
    for name, panel in MODELS.items():
        m, d_tr, _, dg = pipeline_fit(X, cov, panel, ids, ids, t, e, LAM_GRID,
                                      np.random.default_rng(SEEDS[0]))
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
            "coverage": dg.get("score_coverage")}
        print(f"  apparent {name:16s} C={apparent[name]['C']:.4f} "
              f"slope={apparent[name]['calibration_slope']:.3f} "
              f"df={d_tr.shape[1]} lambda={dg['lambda']:.4g}")

    per_seed = {m: [] for m in MODELS}; per_seed["_dC"] = []
    traces = {m: {} for m in MODELS}; nonest = []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        opt = {m: [] for m in MODELS}; opt_sl = {m: [] for m in MODELS}
        origC = {m: [] for m in MODELS}; dC = []
        for b in range(B):
            bi = boot_indices(e, rng); bids = [ids[i] for i in bi]
            tb, eb = t[bi], e[bi]; rep = {}
            for name, panel in MODELS.items():
                try:
                    m, d_tr, d_or, _ = pipeline_fit(X, cov, panel, bids, ids,
                                                    tb, eb, LAM_GRID, rng)
                    lb, lo = linear_predictor(m, d_tr), linear_predictor(m, d_or)
                    cb, co = c_index(tb, eb, lb), c_index(t, e, lo)
                    opt[name].append(cb - co); origC[name].append(co)
                    opt_sl[name].append(calibration_slope(tb, eb, lb)
                                        - calibration_slope(t, e, lo))
                    rep[name] = co
                except (ConvergenceFailure, NonEstimableReplicate) as exc:
                    nonest.append({"model": name, "seed": seed, "replicate": b,
                                   "reason": getattr(exc, "reason", str(exc))})
            if {"M1", "M2_primary"} <= rep.keys():
                dC.append(rep["M2_primary"] - rep["M1"])
        for name in MODELS:
            o = np.array(opt[name], float); sl = np.array(opt_sl[name], float)
            per_seed[name].append({
                "seed": seed, "B_effective": int(len(o)),
                "mean_optimism_C": float(np.nanmean(o)),
                "corrected_C": float(apparent[name]["C"] - np.nanmean(o)),
                "mean_optimism_slope": float(np.nanmean(sl)),
                "corrected_slope": float(apparent[name]["calibration_slope"]
                                         - np.nanmean(sl)),
                "orig_C_p2.5": float(np.nanpercentile(origC[name], 2.5)),
                "orig_C_p97.5": float(np.nanpercentile(origC[name], 97.5))})
            traces[name][str(seed)] = [
                float(apparent[name]["C"] - np.nanmean(o[:k]))
                for k in (50, 100, 200, 300, 400, len(o)) if k <= len(o)]
        per_seed["_dC"].append({
            "seed": seed, "dC_orig_mean": float(np.nanmean(dC)),
            "dC_orig_p2.5": float(np.nanpercentile(dC, 2.5)),
            "dC_orig_p97.5": float(np.nanpercentile(dC, 97.5))})
        print(f"  seed {seed}: " + "  ".join(
            f"{m} C_corr={per_seed[m][-1]['corrected_C']:.4f}" for m in MODELS))

    mc = {}
    for name in MODELS:
        cs = np.array([r["corrected_C"] for r in per_seed[name]])
        boot_sd = float((np.mean([r["orig_C_p97.5"] for r in per_seed[name]])
                         - np.mean([r["orig_C_p2.5"] for r in per_seed[name]])) / 3.92)
        mc_se = float(cs.std(ddof=1))
        mc[name] = {"corrected_C_mean": float(cs.mean()), "monte_carlo_SE": mc_se,
                    "sampling_SE_approx": boot_sd,
                    "MC3_ratio": float(mc_se / boot_sd) if boot_sd else float("inf"),
                    "MC3_pass": bool(boot_sd and mc_se <= 0.10 * boot_sd),
                    "per_seed": per_seed[name], "running_trace_MC1": traces[name]}
        print(f"  MC {name:16s} C_corr={cs.mean():.4f} MC-SE={mc_se:.5f} "
              f"samplingSE~{boot_sd:.5f} ratio={mc[name]['MC3_ratio']:.3f} "
              f"MC-3 {'PASS' if mc[name]['MC3_pass'] else 'FAIL'}")

    n_total = B * len(SEEDS)
    out = {
        "phase": "9", "stage": "EXECUTION", "estimand": "S-7",
        "representation": "representative-IQR (discovery_443), frozen Phase 2B A-05",
        "role": "SENSITIVITY-ONLY -- median aggregation remains PRIMARY (D8D-3)",
        "protocol_version": "3.1",
        "executed_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": round(time.time() - t0, 1),
        "authorised_analyses": list(MODELS), "M3_fitted": False,
        "batch_sensitivity_repeated": False,
        "batch_sensitivity_note": "NOT repeated -- the frozen S-7 definition "
                                  "(Phase 6A A-08) prescribes an aggregation-choice "
                                  "re-fit only, not the D8D-4 batch investigation",
        "endpoint": {"event_col": EV, "time_col": TM},
        "population": {"n": len(ids), "events": int(e.sum()),
                       "censored": int((e == 0).sum())},
        "B": B, "seeds": SEEDS, "lambda_grid": [float(x) for x in LAM_GRID],
        "apparent": apparent, "monte_carlo": mc, "delta_C": per_seed["_dC"],
        "non_estimable_replicates": {
            "records": nonest, "n_total_replicates_per_model": n_total,
            "rate": {m: sum(1 for x in nonest if x["model"] == m) / n_total
                     for m in MODELS},
            "threshold": 0.05, "redrawn": False},
        "prohibited_quantities_not_computed": [
            "t-year / fixed-horizon absolute risk", "CIF at named years",
            "time-dependent AUC", "Brier at a fixed horizon",
            "fixed-horizon cloglog calibration"],
        "input_hashes": {str(Path(p).relative_to(paths.ROOT)): sha(p) for p in (
            MATRIX, ANALYSIS / "covariates.csv", ANALYSIS / "endpoint.csv",
            ANALYSIS / "cohort_membership.csv",
            paths.ROOT / "data/processed/frozen/frozen_gene_to_probe_discovery443.csv",
            paths.ROOT / "data/processed/frozen/selection_population_discovery443.txt",
            PANELS / "m2_panel_4set.json", PANELS / "m2_panel_5set.json",
            paths.ROOT / "src" / "pathway_scoring.py",
            paths.ROOT / "src" / "score_standardisation.py",
            paths.ROOT / "src" / "phase9_execution.py")},
        "environment": environment_report()}
    (OUT / "phase9_s7_results.json").write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"\n  runtime {out['runtime_seconds']}s -> "
          f"{(OUT/'phase9_s7_results.json').relative_to(paths.ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
