"""X7. Exploratory composition-adjusted Cox models. Registered plan §X7.

All models are UNPENALISED Cox fits on the A-19 cohort. R0 is an unpenalised analogue of M2 used only as the reference
for coefficient changes; the frozen ridge M2 estimates are shown alongside and remain the primary result.
Exploratory discrimination: single-level Harrell bootstrap optimism correction (B = 200), every data-dependent step
(knots, gene-level and score standardisation, composition standardisation) refitted inside each resample.
Not comparable with the frozen nested-bootstrap ΔC estimand; no interval is claimed.
"""
import json
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

import xcommon as X
from src.pathway_scoring import score_panel
from src.phase9_execution import build_clinical, rcs_knots, AGE_COL

warnings.filterwarnings("ignore")
B = 200
SCORES = ["EMT", "E2F", "Wnt", "p53"]
SETS = {"ESTIMATE": ["StromalScore", "ImmuneScore"], "MCP": ["MCP_Fibroblasts", "MCP_Endothelial cells", "MCP_Cytotoxic lymphocytes"]}
MODELS = {"C0": ([], []), "C1": ([], SETS["ESTIMATE"]), "R0": (SCORES, []), "R1": (SCORES, SETS["ESTIMATE"]), "R2": (SCORES, SETS["MCP"])}


def features(ctx, comp, train_ids, apply_ids):
    """Design for apply_ids with every parameter estimated on train_ids (may contain duplicates)."""
    age_tr = pd.to_numeric(ctx.cov.loc[train_ids, AGE_COL]).to_numpy(float)
    knots, linear = rcs_knots(age_tr)
    clin_ap = build_clinical(ctx.cov, apply_ids, knots, linear).reset_index(drop=True)
    raw, _ = score_panel(ctx.X, ctx.panel, train_ids=list(train_ids))
    raw = raw.rename(columns=X.LABEL)[SCORES]
    mu, sd = raw.loc[train_ids].mean(), raw.loc[train_ids].std(ddof=1)
    sc_ap = ((raw.loc[apply_ids] - mu) / sd).reset_index(drop=True)
    cmu, csd = comp.loc[train_ids].mean(), comp.loc[train_ids].std(ddof=1)
    cp_ap = ((comp.loc[apply_ids] - cmu) / csd).reset_index(drop=True)
    return pd.concat([clin_ap, sc_ap, cp_ap], axis=1)


def fit(design, cols, t, e):
    df = design[cols].copy()
    df["_T"], df["_E"] = t, e
    return CoxPHFitter().fit(df, "_T", "_E")


def cidx(m, design, cols, t, e):
    lp = design[cols].to_numpy(float) @ m.params_[cols].to_numpy(float)
    return concordance_index(t, -lp, e)


def vif(df):
    Xm = df.to_numpy(float)
    out = {}
    for j, c in enumerate(df.columns):
        y = Xm[:, j]
        A = np.column_stack([np.ones(len(y)), np.delete(Xm, j, axis=1)])
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        r2 = 1 - ((y - A @ beta) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        out[c] = 1 / (1 - r2)
    return out


def main():
    ctx = X.context()
    out = X.outdir("x7_composition_adjusted_cox")
    comp = pd.read_csv(X.OUT / "x5_composition" / "composition_estimates.csv", index_col=0).loc[ctx.ids][SETS["ESTIMATE"] + SETS["MCP"]]
    ids = list(ctx.ids)
    t, e = ctx.t, ctx.e
    D = features(ctx, comp, ids, ids)
    clin_cols = [c for c in D.columns if c not in SCORES and c not in comp.columns]
    cols = {k: clin_cols + s + c for k, (s, c) in MODELS.items()}
    # ---- full-cohort fits
    rows, apparent = [], {}
    for k, cl in cols.items():
        m = fit(D, cl, t, e)
        apparent[k] = cidx(m, D, cl, t, e)
        s = m.summary
        for term in s.index:
            if term in SCORES or term in comp.columns:
                rows.append(dict(model=k, term=term, HR=s.loc[term, "exp(coef)"], lo=s.loc[term, "exp(coef) lower 95%"],
                                 hi=s.loc[term, "exp(coef) upper 95%"], beta=s.loc[term, "coef"], se=s.loc[term, "se(coef)"], p=s.loc[term, "p"]))
    est = pd.DataFrame(rows)
    ref = est[est.model == "R0"].set_index("term")
    est["delta_beta_vs_R0"] = [r.beta - ref.loc[r.term, "beta"] if (r.term in SCORES and r.model in ("R1", "R2")) else np.nan for r in est.itertuples()]
    exr = json.loads((X.ROOT / "reports/phase9/execution/phase9_execution_results.json").read_text())["apparent"]["M2_primary"]
    frozen = {lab: {"HR": exr["hazard_ratios"][pid], "lo": exr["hr_ci_lower"][pid], "hi": exr["hr_ci_upper"][pid]} for pid, lab in X.PATHWAYS}
    est.to_csv(out / "adjusted_estimates.csv", index=False)
    vifs = {"R1": vif(D[SCORES + SETS["ESTIMATE"]]), "R2": vif(D[SCORES + SETS["MCP"]])}
    # ---- deviation 8 (after independent review): likelihood-ratio tests for added blocks, and an MMR-adjusted sensitivity
    from scipy import stats as _st
    ll = {k: fit(D, cl, t, e).log_likelihood_ for k, cl in cols.items()}
    lrt = {}
    for big, small, df in (("R1", "R0", 2), ("R2", "R0", 3), ("C1", "C0", 2), ("R0", "C0", 4)):
        stat_ = 2 * (ll[big] - ll[small])
        lrt[f"{big} vs {small}"] = {"chi2": stat_, "df": df, "p": float(_st.chi2.sf(stat_, df))}
    mmr = ctx.cov.loc[ids, "mmr.status"].astype(str).to_numpy()
    km = np.isin(mmr, ["dMMR", "pMMR"])
    Dm = D.loc[km].reset_index(drop=True).copy()
    Dm["dMMR"] = (mmr[km] == "dMMR").astype(float)
    mmr_rows = []
    for name, cl in (("R0 (MMR-known subset)", cols["R0"]), ("R0 + dMMR", cols["R0"] + ["dMMR"])):
        m = fit(Dm, cl, t[km], e[km])
        for term in SCORES + (["dMMR"] if "dMMR" in cl else []):
            r = m.summary.loc[term]
            mmr_rows.append(dict(model=name, term=term, HR=r["exp(coef)"], lo=r["exp(coef) lower 95%"], hi=r["exp(coef) upper 95%"], p=r["p"]))
    mmr_sens = pd.DataFrame(mmr_rows)
    mmr_sens.to_csv(out / "mmr_adjusted_sensitivity.csv", index=False)
    # ---- bootstrap optimism (single level)
    rng = np.random.default_rng(X.SEED)
    opt = {k: [] for k in cols}
    fails = {k: 0 for k in cols}
    idx = np.arange(len(ids))
    ev = e.astype(bool)
    for b in range(B):
        bi = np.r_[rng.choice(idx[ev], ev.sum(), replace=True), rng.choice(idx[~ev], (~ev).sum(), replace=True)]   # event-stratified
        bids = [ids[i] for i in bi]
        Db = features(ctx, comp, bids, bids)
        Do = features(ctx, comp, bids, ids)
        for k, cl in cols.items():
            try:
                m = fit(Db, cl, t[bi], e[bi])
                opt[k].append(cidx(m, Db, cl, t[bi], e[bi]) - cidx(m, Do, cl, t, e))
            except Exception:
                fails[k] += 1
    disc = pd.DataFrame({k: {"apparent_C": apparent[k], "optimism": float(np.mean(opt[k])), "corrected_C": apparent[k] - float(np.mean(opt[k])),
                             "n_boot_ok": len(opt[k]), "n_fail": fails[k]} for k in cols}).T
    disc.to_csv(out / "exploratory_discrimination.csv")
    X.write_json(out / "x7_summary.json", {"plan_section": "X7", "label": "exploratory outcome analysis; unpenalised; not comparable with frozen nested estimand",
                                          "n": len(ids), "events": int(e.sum()), "models": {k: v for k, v in cols.items()},
                                          "frozen_ridge_M2_HR": frozen, "vif": vifs, "lrt_blocks": lrt,
                                          "mmr_adjusted_sensitivity": {"n": int(km.sum()), "events": int(e[km].sum()), "rows": mmr_sens.to_dict("records"),
                                                                       "caveat": "ComBat protected the combat_mod covariate (aligned with MMR), which can inflate MMR-group differences in the deposited expression (Nygaard et al., Biostatistics 2016)"}, "bootstrap": {"B": B, "seed": X.SEED, "scheme": "event-stratified, single level"},
                                          "discrimination": disc.to_dict("index")})
    print(est.round(3).to_string())
    print(disc.round(4).to_string())
    print({k: {kk: round(vv, 2) for kk, vv in v.items()} for k, v in vifs.items()})
    print(lrt)
    print(mmr_sens.round(3).to_string())


if __name__ == "__main__":
    main()
