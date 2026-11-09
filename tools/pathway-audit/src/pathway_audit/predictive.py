"""Optional: does a score add optimism-corrected discrimination to a clinical Cox model?

Single-level paired Harrell bootstrap: for each resample, both models are refitted (scores re-standardised in the
resample) and evaluated on the resample and on the original data; corrected ΔC = apparent ΔC − mean(optimism difference).
This is a screening estimate. It gives no interval; a confirmatory study should use a nested bootstrap for an interval
on the corrected difference (as in the GSE39582 primary analysis).
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _fit_c(df, cols, t, e, eval_df, eval_t, eval_e):
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    d = df[cols].copy()
    d["_T"], d["_E"] = t, e
    m = CoxPHFitter(penalizer=0.0).fit(d, "_T", "_E")
    lp = eval_df[cols].to_numpy(float) @ m.params_[cols].to_numpy(float)
    return concordance_index(eval_t, -lp, eval_e)


def delta_c(clinical: pd.DataFrame, scores: pd.DataFrame, time, event, B: int = 200, seed: int = 0) -> dict:
    clinical = clinical.reset_index(drop=True).astype(float)
    scores = scores.reset_index(drop=True).astype(float)
    t, e = np.asarray(time, float), np.asarray(event, float)
    ccols, scols = list(clinical.columns), list(scores.columns)

    def design(idx_fit, idx_eval):
        s = scores.iloc[idx_fit]
        mu, sd = s.mean(), s.std(ddof=1)
        fit = pd.concat([clinical.iloc[idx_fit].reset_index(drop=True), ((s - mu) / sd).reset_index(drop=True)], axis=1)
        ev = pd.concat([clinical.iloc[idx_eval].reset_index(drop=True), ((scores.iloc[idx_eval] - mu) / sd).reset_index(drop=True)], axis=1)
        return fit, ev

    all_idx = np.arange(len(t))
    f, ev = design(all_idx, all_idx)
    app = {k: _fit_c(f, c, t, e, ev, t, e) for k, c in (("clinical", ccols), ("clinical+scores", ccols + scols))}
    rng = np.random.default_rng(seed)
    ev_idx, ce_idx = all_idx[e == 1], all_idx[e == 0]
    opt = []
    for _ in range(B):
        bi = np.r_[rng.choice(ev_idx, len(ev_idx)), rng.choice(ce_idx, len(ce_idx))]
        fb, eb = design(bi, bi)
        _, eo = design(bi, all_idx)
        d = 0.0
        for k, c in (("clinical", ccols), ("clinical+scores", ccols + scols)):
            o = _fit_c(fb, c, t[bi], e[bi], eb, t[bi], e[bi]) - _fit_c(fb, c, t[bi], e[bi], eo, t, e)
            d += o if k == "clinical+scores" else -o
        opt.append(d)
    dapp = app["clinical+scores"] - app["clinical"]
    return {"apparent_C": app, "apparent_delta_C": dapp, "optimism_delta": float(np.mean(opt)),
            "corrected_delta_C": float(dapp - np.mean(opt)), "B": B, "n": int(len(t)), "events": int(e.sum()),
            "note": "screening estimate: single-level bootstrap, no interval"}
