"""Phase 5: precision-based feasibility assessment for external validation.

SIMULATION ONLY. Entirely synthetic data. No project data is read, no
predictive model is fitted to real data, no expression data is touched, and
the Phase 3 bundle is not modified.

Framework: Riley et al., Stat Med 2022;41(7):1280-1295 -- target a confidence
interval width for the CALIBRATION SLOPE, the measure that typically needs the
largest sample size. Their stated targets:
    SE 0.051 -> CI width 0.2  (slope 0.9-1.1)   [primary]
    SE 0.102 -> CI width 0.4                     [relaxed]

Usage:
    .venv/bin/python scripts/25_validation_precision_simulation.py [--reps N]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                        # noqa: E402
from src.precision_sim import (cox_univariate, harrell_c,     # noqa: E402
                               n_for_events, simulate)

OUT = paths.ROOT / "reports" / "phase5"
SEED = 20261102

TARGET_SE_PRIMARY = 0.051   # Riley 2022: CI width 0.2
TARGET_SE_RELAXED = 0.102   # Riley 2022: CI width 0.4

# Event counts actually available (Phases 1-4). NOT pooled.
AVAILABLE = {
    "GSE39582 S3 (full follow-up)": 136,
    "GSE39582 S3 (60-unit, conditional)": 130,
    "GSE14333 (external, exploratory)": 48,
    "GSE33113 (external, stage II only, composite)": 18,
}
REFERENCE = [18, 48, 100, 130, 136, 200, 365, 500, 1000]


def calibrate_sigma(target_c: float, event_fraction: float,
                    rng: np.random.Generator) -> float:
    """Find the linear-predictor SD giving a target Harrell C. Bisection."""
    lo, hi = 0.05, 4.0
    for _ in range(40):
        mid = (lo + hi) / 2
        t, e, lp = simulate(40000, mid, event_fraction, rng)
        c = harrell_c(t, e, lp, rng)
        if c < target_c:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def run(events: int, sigma: float, event_fraction: float, reps: int,
        rng: np.random.Generator) -> dict:
    n = n_for_events(events, event_fraction)
    slopes, ses, cs, nev = [], [], [], []
    for _ in range(reps):
        t, e, lp = simulate(n, sigma, event_fraction, rng)
        if e.sum() < 5:
            continue
        b, se = cox_univariate(t, e, lp)
        if not np.isfinite(b) or not np.isfinite(se):
            continue
        slopes.append(b); ses.append(se); nev.append(e.sum())
        cs.append(harrell_c(t, e, lp, rng))
    slopes = np.array(slopes); ses = np.array(ses); cs = np.array(cs)
    return {
        "target_events": events, "n": n,
        "mean_events": float(np.mean(nev)),
        "reps_used": int(len(slopes)),
        "slope_mean": float(np.mean(slopes)),
        "slope_empirical_sd": float(np.std(slopes, ddof=1)),
        "slope_mean_model_se": float(np.mean(ses)),
        "slope_ci_width": float(2 * 1.96 * np.mean(ses)),
        "slope_ci_lo": float(np.mean(slopes) - 1.96 * np.mean(ses)),
        "slope_ci_hi": float(np.mean(slopes) + 1.96 * np.mean(ses)),
        "meets_primary_SE_0051": bool(np.mean(ses) <= TARGET_SE_PRIMARY),
        "meets_relaxed_SE_0102": bool(np.mean(ses) <= TARGET_SE_RELAXED),
        "c_mean": float(np.nanmean(cs)),
        "c_sd": float(np.nanstd(cs, ddof=1)),
        "c_ci_width": float(2 * 1.96 * np.nanstd(cs, ddof=1)),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--event-fraction", type=float, default=0.30,
                    help="assumed proportion with the event by end of follow-up")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    print("=" * 78)
    print("PHASE 5 -- EXTERNAL VALIDATION PRECISION ASSESSMENT (simulation)")
    print("=" * 78)
    print("  Framework: Riley et al. Stat Med 2022;41:1280-1295")
    print(f"  Targets  : SE <= {TARGET_SE_PRIMARY} (CI width 0.2) primary;"
          f" <= {TARGET_SE_RELAXED} (width 0.4) relaxed")
    print(f"  reps     : {args.reps}   seed: {SEED}")
    print(f"  ALL DATA SYNTHETIC. No project data read; no model fitted to real data.")

    print("\n-- ASSUMPTIONS (explicit; every one is a choice) --")
    ef = args.event_fraction
    print(f"  event fraction by end of follow-up : {ef}")
    print("  linear predictor                   : Normal(0, sigma^2)")
    print("  baseline hazard                    : exponential")
    print("  censoring                          : administrative only, at the")
    print("                                       quantile giving the event fraction")
    print("  calibration                        : PERFECT (true slope = 1)")
    print("                                       -> an OPTIMISTIC assumption")

    print("\n-- MODEL STRENGTH CALIBRATION --")
    strengths = {}
    for target_c in (0.60, 0.65, 0.70):
        s = calibrate_sigma(target_c, ef, rng)
        t, e, lp = simulate(40000, s, ef, rng)
        got = harrell_c(t, e, lp, rng)
        strengths[target_c] = s
        print(f"  target C {target_c:.2f} -> sigma_lp = {s:.4f}  (achieved C {got:.4f})")

    results = []
    for target_c, sigma in strengths.items():
        print(f"\n{'='*78}\n  MODEL STRENGTH: C = {target_c:.2f} (sigma_lp = {sigma:.3f})\n{'='*78}")
        print(f"  {'events':>7} {'n':>7} {'slope SE':>9} {'slope 95% CI':>18} "
              f"{'width':>7} {'SE<=.051':>9} {'SE<=.102':>9} {'C 95% CI width':>15}")
        for ev in REFERENCE:
            r = run(ev, sigma, ef, args.reps, rng)
            r["target_c"] = target_c
            r["sigma_lp"] = sigma
            r["event_fraction"] = ef
            results.append(r)
            print(f"  {r['target_events']:>7} {r['n']:>7} {r['slope_mean_model_se']:>9.4f} "
                  f"{'[' + format(r['slope_ci_lo'], '.3f') + ', ' + format(r['slope_ci_hi'], '.3f') + ']':>18} "
                  f"{r['slope_ci_width']:>7.3f} "
                  f"{'YES' if r['meets_primary_SE_0051'] else 'no':>9} "
                  f"{'YES' if r['meets_relaxed_SE_0102'] else 'no':>9} "
                  f"{r['c_ci_width']:>15.4f}")

    df = pd.DataFrame(results)
    df.to_csv(OUT / "precision_simulation.csv", index=False)

    print("\n" + "=" * 78)
    print("  AVAILABLE COHORTS -- not pooled, assessed separately")
    print("=" * 78)
    summary = {}
    for label, ev in AVAILABLE.items():
        sub = df[df.target_events == min(REFERENCE, key=lambda x: abs(x - ev))]
        rows = []
        for _, r in sub.iterrows():
            rows.append({"target_c": r.target_c,
                         "slope_SE": round(r.slope_mean_model_se, 4),
                         "slope_CI_width": round(r.slope_ci_width, 3),
                         "C_CI_width": round(r.c_ci_width, 4),
                         "meets_primary": bool(r.meets_primary_SE_0051),
                         "meets_relaxed": bool(r.meets_relaxed_SE_0102)})
        summary[label] = {"events": ev, "by_model_strength": rows}
        print(f"\n  {label}  ({ev} events)")
        for x in rows:
            print(f"    C={x['target_c']:.2f}: slope SE {x['slope_SE']:.4f}, "
                  f"CI width {x['slope_CI_width']:.3f}, C CI width {x['C_CI_width']:.4f}  "
                  f"-> primary {'MET' if x['meets_primary'] else 'NOT MET'}, "
                  f"relaxed {'MET' if x['meets_relaxed'] else 'NOT MET'}")

    payload = {
        "framework": "Riley et al. Stat Med 2022;41(7):1280-1295",
        "targets": {"primary_SE": TARGET_SE_PRIMARY,
                    "primary_CI_width": 0.2,
                    "relaxed_SE": TARGET_SE_RELAXED,
                    "relaxed_CI_width": 0.4},
        "assumptions": {
            "linear_predictor": "Normal(0, sigma^2)",
            "baseline_hazard": "exponential",
            "censoring": "administrative only",
            "event_fraction": ef,
            "true_calibration_slope": 1.0,
            "note": "perfect calibration assumed -> OPTIMISTIC; a miscalibrated "
                    "model needs MORE events, not fewer",
        },
        "reps": args.reps, "seed": SEED,
        "model_strengths": {str(k): v for k, v in strengths.items()},
        "available_cohorts": summary,
        "results": results,
        "data_source": "ENTIRELY SYNTHETIC",
    }
    (OUT / "precision_simulation.json").write_text(
        json.dumps(payload, indent=2, default=str) + "\n")
    print(f"\n  -> {OUT.relative_to(paths.ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
