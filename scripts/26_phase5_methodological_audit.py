"""Phase 5 audit: (a) divergences from Riley 2022; (b) C-index CI estimator
audit; (c) transcriptome-wide selection risk at 136 events.

ENTIRELY SYNTHETIC. No project data read, no model fitted to real data,
Phase 3 bundle untouched.

Usage:
    .venv/bin/python scripts/26_phase5_methodological_audit.py [--reps N]
"""

import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                      # noqa: E402
from scipy.stats import norm                               # noqa: E402
from src.precision_sim import harrell_c, simulate          # noqa: E402

OUT = paths.ROOT / "reports" / "phase5"
SEED = 20261102


def section(t): print("\n" + "=" * 78); print(t); print("=" * 78)


def univariate_score_tests(time, event, Z):
    """Cox score statistics at beta=0 for every column of Z, vectorised.

    At beta=0 all weights are 1, so risk-set sums are reverse cumulative sums
    over time-sorted rows. Returns |z-statistic| per predictor.
    """
    order = np.argsort(time)
    d = event[order].astype(bool)
    X = Z[order]                                   # n x p
    n = X.shape[0]
    s0 = np.arange(n, 0, -1).astype(float)         # risk-set sizes
    s1 = np.cumsum(X[::-1], axis=0)[::-1]          # n x p
    s2 = np.cumsum((X ** 2)[::-1], axis=0)[::-1]
    mu = s1 / s0[:, None]
    U = (X[d] - mu[d]).sum(axis=0)
    info = (s2[d] / s0[d, None] - mu[d] ** 2).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.abs(U / np.sqrt(info)), U, info


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reps", type=int, default=40)
    ap.add_argument("--genes", type=int, default=21815)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    findings = {}

    # ---------------------------------------------------------------
    section("A. C-INDEX CONFIDENCE-INTERVAL ESTIMATOR AUDIT")
    # ---------------------------------------------------------------
    print("  Estimator used in script 25 : Harrell's C (concordant/comparable")
    print("                                pairs, ties credited 0.5)")
    print("  Censoring assumption        : pairs comparable only when the")
    print("                                EARLIER time is an observed event")
    print("                                -> Harrell's C depends on the")
    print("                                   censoring distribution (known")
    print("                                   limitation; Uno's IPCW C is the")
    print("                                   censoring-robust alternative --")
    print("                                   NOT implemented)")
    print("  CI method in script 25      : 1.96 x empirical SD of C across")
    print("                                replications")
    print("  What that quantity IS       : a Monte-Carlo estimate of the")
    print("                                SAMPLING SD of the estimator")
    print("  What Riley reports instead  : the MEAN WITHIN-DATASET 95% CI width")
    print("  => related but NOT the same quantity. Documented as an")
    print("     approximation, not as Riley's measure.\n")

    print("  -- empirical check: does 1.96*SD track a within-dataset CI? --")
    rows = []
    for ev in (48, 136, 365):
        n = int(round(ev / 0.30))
        cs = []
        for _ in range(300):
            t, e, lp = simulate(n, 0.555, 0.30, rng)
            cs.append(harrell_c(t, e, lp, rng))
        cs = np.array(cs)
        # bootstrap a single dataset for a within-dataset CI
        t, e, lp = simulate(n, 0.555, 0.30, rng)
        boot = []
        for _ in range(300):
            idx = rng.integers(0, n, n)
            c = harrell_c(t[idx], e[idx], lp[idx], rng)
            if np.isfinite(c):
                boot.append(c)
        boot = np.array(boot)
        bw = float(np.percentile(boot, 97.5) - np.percentile(boot, 2.5))
        mc = float(2 * 1.96 * np.std(cs, ddof=1))
        rows.append({"events": ev, "monte_carlo_width": round(mc, 4),
                     "bootstrap_within_dataset_width": round(bw, 4),
                     "ratio": round(bw / mc, 3)})
        print(f"    {ev:4d} events: MC width {mc:.4f} | bootstrap within-dataset "
              f"width {bw:.4f} | ratio {bw/mc:.3f}")
    findings["c_index_audit"] = {
        "estimator": "Harrell's C", "censoring_robust": False,
        "alternative_not_implemented": "Uno's IPCW C",
        "ci_method_script25": "1.96 x empirical SD across replications "
                              "(Monte-Carlo sampling SD)",
        "riley_reports": "mean within-dataset 95% CI width",
        "equivalent": False, "replicates_script25": 400,
        "comparison": rows,
    }
    pd.DataFrame(rows).to_csv(OUT / "c_index_ci_audit.csv", index=False)

    # ---------------------------------------------------------------
    section("B. TRANSCRIPTOME-WIDE SELECTION RISK AT 136 EVENTS (null genes)")
    # ---------------------------------------------------------------
    print(f"  Design: {args.genes} predictors, ALL NULL (no association with")
    print("  outcome by construction). True C-index = 0.500 for any model.")
    print("  n = 453, ~136 events -- the S3 cohort's actual size.\n")
    n, ef = 453, 0.30
    res = []
    for k in (5, 10, 25, 50, 100):
        apparent = []
        for _ in range(args.reps):
            t, e, lp0 = simulate(n, 0.0, ef, rng)       # sigma=0 -> no signal
            Z = rng.standard_normal((n, args.genes))
            zstat, U, info = univariate_score_tests(t, e, Z)
            top = np.argsort(zstat)[::-1][:k]
            sign = np.sign(U[top])
            score = (Z[:, top] * sign).sum(axis=1)
            apparent.append(harrell_c(t, e, score, rng))
        a = np.array(apparent)
        res.append({"genes_selected": k,
                    "apparent_C_mean": round(float(np.mean(a)), 4),
                    "apparent_C_sd": round(float(np.std(a, ddof=1)), 4),
                    "true_C": 0.5,
                    "optimism": round(float(np.mean(a) - 0.5), 4)})
        print(f"    top {k:3d} genes selected: apparent C = {np.mean(a):.4f} "
              f"(SD {np.std(a, ddof=1):.4f})  ->  OPTIMISM = "
              f"{np.mean(a)-0.5:+.4f} on a TRUE C of 0.500")
    print("\n  Every gene is pure noise. Any apparent discrimination above")
    print("  0.500 is selection-induced optimism and nothing else.")
    pd.DataFrame(res).to_csv(OUT / "selection_optimism_null.csv", index=False)
    findings["selection_risk"] = {
        "n": n, "events_approx": 136, "genes": args.genes,
        "all_predictors_null": True, "reps": args.reps, "results": res,
    }

    # ---------------------------------------------------------------
    section("C. MAXIMUM UNIVARIATE SCORE STATISTIC UNDER THE NULL")
    # ---------------------------------------------------------------
    t, e, _ = simulate(n, 0.0, ef, rng)
    Z = rng.standard_normal((n, args.genes))
    zstat, _, _ = univariate_score_tests(t, e, Z)
    qs = np.percentile(zstat, [50, 95, 99, 99.9])
    print(f"  {args.genes} null genes, one simulated cohort (~136 events):")
    print(f"    median |z|   : {qs[0]:.3f}")
    print(f"    95th pct |z| : {qs[1]:.3f}")
    print(f"    99.9th pct   : {qs[3]:.3f}")
    print(f"    MAXIMUM |z|  : {zstat.max():.3f}   (nominal two-sided p = {2*norm.sf(zstat.max()):.2e})")
    print(f"    genes with |z| > 1.96 : {(zstat>1.96).sum()} "
          f"({(zstat>1.96).mean()*100:.1f}%, expected ~5%)")
    print(f"    genes with |z| > 3.00 : {(zstat>3.0).sum()}")
    print("\n  => with ~21,800 null genes a handful ALWAYS reach |z|>4 by chance.")
    print("     Univariate screening cannot distinguish them from real signal")
    print("     at this event count.")
    findings["null_max_z"] = {
        "median": float(qs[0]), "pct95": float(qs[1]), "pct999": float(qs[3]),
        "max": float(zstat.max()),
        "n_above_1.96": int((zstat > 1.96).sum()),
        "n_above_3": int((zstat > 3.0).sum()),
    }

    (OUT / "methodological_audit.json").write_text(
        json.dumps(findings, indent=2, default=str) + "\n")
    section("DONE")
    print(f"  -> {OUT.relative_to(paths.ROOT)}/")
    print("  All data synthetic. No real data read; no model fitted to real data.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
