"""X1. Molecular anchoring of the four frozen pathway scores (outcome-free). Registered plan §X1."""
import numpy as np
import pandas as pd
from scipy import stats

import xcommon as X

ANN = {"tp53.mutation": ("M", "WT"), "kras.mutation": ("M", "WT"), "braf.mutation": ("M", "WT"),
       "cimp.status": ("+", "-"), "cin.status": ("+", "-"), "cit.molecularsubtype": None}


def main():
    ctx = X.context()
    z, _ = X.cohort_scores(ctx)
    cov = ctx.cov.loc[ctx.ids].copy()
    soft = X.soft_characteristics()
    cov["cit.molecularsubtype"] = soft.loc[ctx.ids, "cit.molecularsubtype"]
    _v = cov["tp53.mutation"].astype(str) != ""
    assert (soft.loc[ctx.ids, "tp53.mutation"][_v].astype(str) == cov["tp53.mutation"][_v].astype(str)).all()
    out = X.outdir("x1_molecular_anchoring")
    rows, miss = [], {}
    for ann, levels in ANN.items():
        v = cov[ann].astype(str).str.strip()
        valid = v.isin(levels) if levels else (~v.isin(["", "N/A", "NA", "nan"]))
        miss[ann] = {"n_valid": int(valid.sum()), "n_missing_or_NA": int((~valid).sum()),
                     "levels": v[valid].value_counts().to_dict()}
        for lab in z.columns:
            s = z.loc[valid, lab]
            g = v[valid]
            if levels:
                a, b = s[g == levels[0]], s[g == levels[1]]
                U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                rbc = 2 * U / (len(a) * len(b)) - 1          # rank-biserial: P(a>b) - P(a<b)
                rows.append(dict(annotation=ann, pathway=lab, test="Mann-Whitney U", group1=levels[0], n1=len(a),
                                 group2=levels[1], n2=len(b), median1=a.median(), median2=b.median(),
                                 median_diff_sd=a.median() - b.median(), effect=rbc, effect_name="rank-biserial",
                                 statistic=U, p=p))
            else:
                groups = [s[g == k] for k in sorted(g.unique())]
                H, p = stats.kruskal(*groups)
                n, k = len(s), len(groups)
                eps2 = (H - k + 1) / (n - k)
                rows.append(dict(annotation=ann, pathway=lab, test="Kruskal-Wallis", group1="; ".join(f"{k_}:{(g == k_).sum()}" for k_ in sorted(g.unique())),
                                 n1=n, group2="", n2=np.nan, median1=np.nan, median2=np.nan, median_diff_sd=np.nan,
                                 effect=eps2, effect_name="epsilon-squared", statistic=H, p=p))
    res = pd.DataFrame(rows)
    res["q_BH"] = X.bh(res["p"].to_numpy())
    # deviation 7 (after independent review): CIT subtypes were derived from these arrays (circular); BH also reported over the 20 non-CIT tests
    nc = res.annotation != "cit.molecularsubtype"
    res["q_BH_20_excl_CIT"] = np.nan
    res.loc[nc, "q_BH_20_excl_CIT"] = X.bh(res.loc[nc, "p"].to_numpy())
    res.to_csv(out / "anchoring_tests.csv", index=False)
    # subtype medians for interpretation
    cit = cov["cit.molecularsubtype"].astype(str)
    okc = ~cit.isin(["", "N/A", "NA"])
    med = z.loc[okc].groupby(cit[okc]).median().round(3)
    med["n"] = cit[okc].value_counts()
    med.to_csv(out / "cit_subtype_score_medians.csv")

    # pre-specified secondary: TP53 within pMMR, batch-stratified permutation
    part = pd.read_csv(X.ANALYSIS / "partition_and_batch.csv").set_index("sample_id").loc[ctx.ids]
    tp = cov["tp53.mutation"].astype(str)
    mmr = cov["mmr.status"].astype(str)
    sel = tp.isin(["M", "WT"]) & (mmr == "pMMR")
    s, gm, bt = z.loc[sel, "p53"].to_numpy(), (tp[sel] == "M").to_numpy(), part.loc[sel, "combat_batch"].to_numpy()

    def stat(lbl):
        return s[lbl].mean() - s[~lbl].mean()
    obs = stat(gm)
    rng = np.random.default_rng(X.SEED)
    perm = np.empty(10000)
    for i in range(10000):
        lab = gm.copy()
        for b in np.unique(bt):
            idx = np.where(bt == b)[0]
            lab[idx] = rng.permutation(lab[idx])
        perm[i] = stat(lab)
    p_perm = (1 + np.sum(np.abs(perm) >= abs(obs))) / (1 + len(perm))
    U, p_mw = stats.mannwhitneyu(s[gm], s[~gm], alternative="two-sided")
    sec = {"population": "pMMR tumours with TP53 status", "n_mutant": int(gm.sum()), "n_wildtype": int((~gm).sum()),
           "mean_diff_sd_mutant_minus_wt": obs, "batch_stratified_permutation_p": p_perm, "n_permutations": 10000,
           "mann_whitney_p": p_mw, "rank_biserial": 2 * U / (gm.sum() * (~gm).sum()) - 1}
    # TP53 x MMR cross-tab (explains confounding)
    xt = pd.crosstab(tp[tp.isin(["M", "WT"])], mmr[tp.isin(["M", "WT"])])
    X.write_json(out / "x1_summary.json", {"plan_section": "X1", "n_cohort": len(ctx.ids), "missingness": miss,
                                          "n_primary_tests": len(res), "multiplicity": "BH-FDR across 24 tests",
                                          "secondary_tp53_within_pMMR": sec, "tp53_by_mmr": xt.to_dict(),
                                          "label": "exploratory, outcome-free"})
    print(res[["annotation", "pathway", "n1", "n2", "median_diff_sd", "effect", "p", "q_BH"]].round(4).to_string())
    print(sec)
    print(xt)
    print(med)


if __name__ == "__main__":
    main()
