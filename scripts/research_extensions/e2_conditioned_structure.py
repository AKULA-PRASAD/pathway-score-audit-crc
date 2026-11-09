"""E2. Composition- and subtype-conditioned structure (registered extension E2). Outcome-free; estimation only."""
import json

import numpy as np
import pandas as pd
from scipy import stats

import rcommon as C
from rcommon import X

SC = ["EMT", "E2F", "Wnt", "p53"]


def rank_resid(y, Z):
    yr = stats.rankdata(y)
    Zr = np.column_stack([np.ones(len(y))] + [stats.rankdata(z) for z in Z.T])
    b, *_ = np.linalg.lstsq(Zr, yr, rcond=None)
    return yr - Zr @ b


def partial_rho(a, b, Z):
    return np.corrcoef(rank_resid(a, Z), rank_resid(b, Z))[0, 1]


def ols_coef(y, x, Z):
    A = np.column_stack([np.ones(len(y)), x] + ([Z] if Z is not None and Z.size else []))
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    return b[1]


def main():
    C.dirs()
    ctx = X.context()
    z, _ = X.cohort_scores(ctx)
    comp = pd.read_csv(X.OUT / "x5_composition" / "composition_estimates.csv", index_col=0).loc[ctx.ids]
    cov = ctx.cov.loc[ctx.ids]
    cit = X.soft_characteristics().loc[ctx.ids, "cit.molecularsubtype"].to_numpy()
    rng = np.random.default_rng(C.SEED)
    Zc = comp[["StromalScore", "ImmuneScore"]].to_numpy()
    # 1. partial correlations
    pr = []
    for i, a in enumerate(SC):
        for b in SC[i + 1:]:
            raw = C.boot_ci(lambda x, y: stats.spearmanr(x, y)[0], [z[a].to_numpy(), z[b].to_numpy()], rng)
            par = C.boot_ci(lambda x, y, w: partial_rho(x, y, w), [z[a].to_numpy(), z[b].to_numpy(), Zc], rng)
            pr.append(dict(pair=f"{a}–{b}", rho=raw[0], rho_lo=raw[1], rho_hi=raw[2], partial_rho=par[0], partial_lo=par[1], partial_hi=par[2]))
    pr = pd.DataFrame(pr)
    pr.to_csv(C.RES / "e2_partial_correlations.csv", index=False)
    # 2. genotype effects adjusted for composition (+ dMMR where specified)
    mmr = cov["mmr.status"].astype(str).to_numpy()
    rows = []
    for ann, pos, neg, add_mmr in [("tp53.mutation", "M", "WT", True), ("braf.mutation", "M", "WT", True), ("cimp.status", "+", "-", True),
                                   ("kras.mutation", "M", "WT", False), ("cin.status", "+", "-", False)]:
        v = cov[ann].astype(str).to_numpy()
        sel = np.isin(v, [pos, neg]) & (np.isin(mmr, ["dMMR", "pMMR"]) if add_mmr else True)
        x = (v[sel] == pos).astype(float)
        adj = Zc[sel]
        if add_mmr:
            adj = np.column_stack([adj, (mmr[sel] == "dMMR").astype(float)])
        for s in SC:
            y = z[s].to_numpy()[sel]
            un = C.boot_ci(lambda yy, xx: ols_coef(yy, xx, None), [y, x], rng)
            ad = C.boot_ci(lambda yy, xx, ww: ols_coef(yy, xx, ww), [y, x, adj], rng)
            rows.append(dict(annotation=ann, score=s, n=int(sel.sum()), n_pos=int(x.sum()), adjusted_for="stromal+immune" + ("+dMMR" if add_mmr else ""),
                             diff_unadj=un[0], unadj_lo=un[1], unadj_hi=un[2], diff_adj=ad[0], adj_lo=ad[1], adj_hi=ad[2],
                             robust_to_composition=bool((un[1] > 0 or un[2] < 0) and (ad[1] > 0 or ad[2] < 0) and np.sign(un[0]) == np.sign(ad[0]))))
    gt = pd.DataFrame(rows)
    gt.to_csv(C.RES / "e2_genotype_adjusted.csv", index=False)
    # 3. within-subtype EMT–stroma
    ws = []
    for k in sorted(np.unique(cit)):
        m = cit == k
        for c in ("MCP_Fibroblasts", "StromalScore"):
            r = C.boot_ci(lambda x, y: stats.spearmanr(x, y)[0], [z["EMT"].to_numpy()[m], comp[c].to_numpy()[m]], rng)
            ws.append(dict(cit=k, n=int(m.sum()), composition=c, rho=r[0], lo=r[1], hi=r[2]))
    ws = pd.DataFrame(ws)
    ws.to_csv(C.RES / "e2_within_subtype_emt_stroma.csv", index=False)
    fib = ws[ws.composition == "MCP_Fibroblasts"]
    n_ok = int((fib.rho >= 0.7).sum())
    X.write_json(C.REP / "e2_summary.json", {"extension": "E2", "label": "post-hoc extension, outcome-free; estimation only", "n": len(ctx.ids),
                                             "bootstrap": {"B": 2000, "seed": C.SEED},
                                             "within_subtype_rule": {"criterion": "rho(EMT, MCP fibroblasts) >= 0.7 in >= 5 of 6 CIT subtypes", "n_subtypes_meeting": n_ok,
                                                                     "met": n_ok >= 5}})
    print(pr.round(3).to_string())
    print(gt.round(3).to_string())
    print(ws.round(3).to_string())
    print("subtypes meeting:", n_ok)


if __name__ == "__main__":
    main()
