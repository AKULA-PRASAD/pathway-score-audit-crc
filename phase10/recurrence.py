"""Track A: gene-expression survival ML for cause-specific recurrence (exploratory, post hoc).

Cohort, endpoint and clinical coding come from the frozen Phase 9 loader. Strict nesting: every
data-dependent step (spline knots, pathway-score gene moments and score standardisation, variance
filter, gene scaling) is re-fitted inside each inner-training fold during tuning and on the outer
training set for the final fit. Metric: Harrell's C, within-fold pooled over identical comparable
pairs. No fixed-time metric is computed: the time unit of rfs.delay is undocumented.

Usage:  python -m phase10.recurrence [--benchmark] [--repeats N] [--null N]
"""
from __future__ import annotations

import argparse
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from phase10 import common as C
from src.pathway_scoring import score_panel
from src.phase9_corrected.data import LAM_GRID, load_context
from src.phase9_execution import (AGE_COL, ConvergenceFailure, build_clinical, c_index, fit_cox,
                                  linear_predictor, rcs_knots)
from src.score_standardisation import standardise_scores_training_fold

OUT = C.REPORTS / "recurrence"
PHASE9 = {"estimand": "bootstrap optimism-corrected Harrell C on the full cohort (Phase 9, frozen)",
          "M1": 0.6016131687304218, "M2_four_set": 0.6059295115344525, "delta": 0.004316342804,
          "source": "reports/phase9/corrected_rerun_report/corrected_interval_results.json"}


# ---------------------------------------------------------------- designs (frozen functions)
def clinical_design(ctx, tr, ap):
    knots, lin = rcs_knots(pd.to_numeric(ctx.cov.loc[tr, AGE_COL]).to_numpy(float))
    return build_clinical(ctx.cov, tr, knots, lin), build_clinical(ctx.cov, ap, knots, lin)


def m2_design(ctx, tr, ap):
    """Same construction as frozen pipeline_fit: knots, scores and standardisation from `tr`."""
    d_tr, d_ap = clinical_design(ctx, tr, ap)
    sc_all, _ = score_panel(ctx.X, ctx.panel, train_ids=list(tr))
    std_all, _ = standardise_scores_training_fold(sc_all, list(tr))
    d_tr = pd.concat([d_tr, std_all.loc[tr].set_axis(d_tr.index)], axis=1)
    d_ap = pd.concat([d_ap, std_all.loc[ap].set_axis(d_ap.index)], axis=1)
    return d_tr, d_ap, 6 if d_tr.shape[1] >= 6 else d_tr.shape[1] - 4


def n_clinical(d):
    return sum(not c.startswith("HALLMARK") for c in d.columns)


def m2_fit(ctx, tr, t_tr, e_tr, lam):
    d_tr, _, _ = m2_design(ctx, tr, tr)
    nc = n_clinical(d_tr)
    pen = np.r_[np.zeros(nc), np.full(d_tr.shape[1] - nc, lam)]
    return fit_cox(d_tr, t_tr, e_tr, pen)


def permute_pairs(t, e, perm_seed):
    """One permutation applied to (time, event) together: each patient's pair moves as a unit."""
    pi = np.random.default_rng(perm_seed).permutation(len(t))
    return np.asarray(t)[pi], np.asarray(e)[pi]


def _idx(ids, sub):
    pos = {s: i for i, s in enumerate(ids)}
    return np.array([pos[s] for s in sub])


# ---------------------------------------------------------------- models for one outer fold
def inner_splits(e_tr, base, r, f, k):
    skf = StratifiedKFold(k, shuffle=True, random_state=C.seed(base, 2, r, f, 99))
    return list(skf.split(np.zeros(len(e_tr)), e_tr))


def tune_m2(ctx, tr, t_tr, e_tr, splits):
    scores = {lam: [] for lam in LAM_GRID}
    for itr, iva in splits:
        if e_tr[itr].sum() < 5 or e_tr[iva].sum() < 2:          # frozen skip rule
            continue
        a, b = [tr[i] for i in itr], [tr[i] for i in iva]
        d_a, d_b, _ = m2_design(ctx, a, b)
        nc = n_clinical(d_a)
        for lam in LAM_GRID:
            pen = np.r_[np.zeros(nc), np.full(d_a.shape[1] - nc, lam)]
            try:
                m = fit_cox(d_a, t_tr[itr], e_tr[itr], pen)
                scores[lam].append(c_index(t_tr[iva], e_tr[iva], linear_predictor(m, d_b)))
            except ConvergenceFailure:
                scores[lam].append(np.nan)
    best, best_s = LAM_GRID[0], -np.inf
    for lam in LAM_GRID:                                          # strict '>' keeps the smaller lambda
        s = np.nanmean(scores[lam]) if np.any(np.isfinite(scores[lam])) else -np.inf
        if s > best_s:
            best, best_s = lam, s
    return float(best), float(best_s)


def _surv(t, e):
    from sksurv.util import Surv
    return Surv.from_arrays(event=np.asarray(e).astype(bool), time=np.asarray(t, float))


def m3_design(ctx, Xs, a, b, k):
    """Clinical 6 df (unpenalised) + top-k variance genes; all scaling fitted on `a`."""
    d_a, d_b = clinical_design(ctx, a, b)
    vf = C.VarianceFilter(k).fit(Xs.loc[a])
    g_a, g_b = vf.transform(Xs.loc[a]), vf.transform(Xs.loc[b])
    sc = StandardScaler().fit(np.hstack([d_a.to_numpy(), g_a.to_numpy()]))
    A = sc.transform(np.hstack([d_a.to_numpy(), g_a.to_numpy()]))
    B = sc.transform(np.hstack([d_b.to_numpy(), g_b.to_numpy()]))
    return A, B, d_a.shape[1], vf.columns_


def coxnet(l1, alphas, n_clin, n_gene):
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    return CoxnetSurvivalAnalysis(l1_ratio=l1, alphas=list(alphas), penalty_factor=np.r_[np.zeros(n_clin), np.ones(n_gene)],
                                  normalize=False, fit_baseline_model=False, max_iter=100000)


def tune_m3(ctx, Xs, tr, t_tr, e_tr, splits, l1s, alphas, k):
    res = {(l1, a): [] for l1 in l1s for a in alphas}
    failures = []
    for j, (itr, iva) in enumerate(splits):
        a_ids, b_ids = [tr[i] for i in itr], [tr[i] for i in iva]
        A, B, nc, cols = m3_design(ctx, Xs, a_ids, b_ids, k)
        for l1 in l1s:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    m = coxnet(l1, alphas, nc, len(cols)).fit(A, _surv(t_tr[itr], e_tr[itr]))
                for al in alphas:
                    res[(l1, al)].append(c_index(t_tr[iva], e_tr[iva], m.predict(B, alpha=al)))
            except Exception as exc:
                failures.append(f"inner{j} l1={l1}: {type(exc).__name__}: {exc}")
                for al in alphas:
                    res[(l1, al)].append(np.nan)
    mean = {key: (np.nanmean(v) if np.any(np.isfinite(v)) else -np.inf) for key, v in res.items()}
    top = max(mean.values())
    # ties -> larger alpha (sparser), then larger l1_ratio
    best = max((key for key, v in mean.items() if v == top), key=lambda kk: (kk[1], kk[0]))
    return best, float(top), failures


def outer_fold(ctx, Xs, cfg, folds_rf, r, f, perm_seed=None, models=("M1", "M1_breslow", "M2_strict", "M3_enet_cox")):
    """Fit every requested model on one outer training set and score the untouched test fold."""
    ids = ctx.ids
    te = folds_rf.index[folds_rf == f].tolist()
    tr = folds_rf.index[folds_rf != f].tolist()
    t, e = ctx.t, ctx.e
    t_tr, e_tr = t[_idx(ids, tr)].copy(), e[_idx(ids, tr)].copy()
    t_te, e_te = t[_idx(ids, te)], e[_idx(ids, te)]
    if perm_seed is not None:                   # N1: permute (time, event) PAIRS within training only
        t_tr, e_tr = permute_pairs(t_tr, e_tr, perm_seed)
    base, k_in, k = cfg["seed_base"], cfg["inner_folds"], cfg["variance_top_k"]
    splits = inner_splits(e_tr, base, r, f, k_in)
    rows, extra = [], {}
    for name in models:
        t0 = time.perf_counter()
        status, risk, info = "ok", None, {}
        try:
            if name == "M1":
                d_tr, d_te = clinical_design(ctx, tr, te)
                m = fit_cox(d_tr, t_tr, e_tr, np.zeros(d_tr.shape[1]))
                risk = linear_predictor(m, d_te)
            elif name == "M1_breslow":
                from sksurv.linear_model import CoxPHSurvivalAnalysis
                d_tr, d_te = clinical_design(ctx, tr, te)
                m = CoxPHSurvivalAnalysis(alpha=0.0, ties="breslow").fit(d_tr.to_numpy(), _surv(t_tr, e_tr))
                risk = m.predict(d_te.to_numpy())
            elif name == "M2_strict":
                lam, inner = tune_m2(ctx, tr, t_tr, e_tr, splits)
                d_tr, d_te, _ = m2_design(ctx, tr, te)
                nc = n_clinical(d_tr)
                m = fit_cox(d_tr, t_tr, e_tr, np.r_[np.zeros(nc), np.full(d_tr.shape[1] - nc, lam)])
                risk = linear_predictor(m, d_te)
                info = {"lambda": lam, "inner_c": inner}
            elif name == "M3_enet_cox":
                l1s = cfg["recurrence"]["enet_cox"]["l1_ratio"]
                alphas = C.grid(cfg["recurrence"]["enet_cox"]["alphas"])
                (l1, al), inner, fails = tune_m3(ctx, Xs, tr, t_tr, e_tr, splits, l1s, alphas, k)
                A, B, nc, cols = m3_design(ctx, Xs, tr, te, k)
                path = [a for a in alphas if a >= al]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    m = coxnet(l1, path, nc, len(cols)).fit(A, _surv(t_tr, e_tr))
                risk = m.predict(B, alpha=al)
                coef = m.coef_[:, -1]
                sel = {cols[i - nc]: float(coef[i]) for i in np.flatnonzero(coef[nc:]) + nc}
                info = {"l1_ratio": l1, "alpha": al, "inner_c": inner, "alpha_at_grid_edge": al in (alphas[0], alphas[-1]),
                        "n_genes_selected": len(sel), "selected": sel, "inner_failures": fails}
        except Exception as exc:
            status = f"failed: {type(exc).__name__}: {exc}"
        secs = time.perf_counter() - t0
        if status == "ok":
            num, den = C.harrell_pair_counts(t_te, e_te, risk)
            ref = c_index(t_te, e_te, risk)
        else:
            num = den = ref = np.nan
        rows.append({"model": name, "repeat": r, "fold": f, "n_test": len(te), "events_test": int(e_te.sum()),
                     "c_num": num, "c_den": den, "fold_c": num / den if den else np.nan, "fold_c_lifelines": ref,
                     "status": status, "seconds": round(secs, 1), "null_seed": perm_seed})
        extra[name] = info
    return rows, {"repeat": r, "fold": f, "null_seed": perm_seed, **extra}


def summarise(fold_df):
    per_rep = (fold_df[fold_df.status == "ok"].groupby(["model", "repeat"])
               .apply(lambda g: g.c_num.sum() / g.c_den.sum(), include_groups=False).rename("c_within_fold_pooled").reset_index())
    out, deltas = {}, {}
    for m, g in per_rep.groupby("model"):
        out[m] = C.summarise(g.c_within_fold_pooled)
    for m, ref in [("M1_breslow", "M1"), ("M2_strict", "M1"), ("M3_enet_cox", "M1"), ("M3_enet_cox", "M1_breslow")]:
        ds = []
        for r in sorted(fold_df.repeat.unique()):
            a = fold_df[(fold_df.model == m) & (fold_df.repeat == r)].set_index("fold")
            b = fold_df[(fold_df.model == ref) & (fold_df.repeat == r)].set_index("fold")
            if a.empty or b.empty:
                continue
            both = a.index[(a.status == "ok") & (b.loc[a.index].status == "ok")]
            ds.append((a.loc[both].c_num.sum() - b.loc[both].c_num.sum()) / a.loc[both].c_den.sum())
        deltas[f"{m}_minus_{ref}"] = C.summarise(ds)
    return per_rep, out, deltas


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", action="store_true")
    ap.add_argument("--repeats", type=int, default=None)
    ap.add_argument("--null", type=int, default=0, help="N1 replicates (training-fold (time,event) permutation, M3 only)")
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args(argv)
    cfg = C.load_config()
    ctx = load_context(C.ROOT)
    Xs = ctx.X.T
    reps = a.repeats or cfg["repeats"]
    folds = C.make_folds(ctx.ids, ctx.e, "recurrence", reps, cfg["outer_folds"], cfg["seed_base"])
    print(f"Track A cohort: n={len(ctx.ids)} events={int(ctx.e.sum())} censored={int((1 - ctx.e).sum())}")
    fr = {r: folds[folds.repeat == r].set_index("sample_id")["fold"] for r in range(reps)}
    if a.benchmark:
        t0 = time.perf_counter()
        rows, info = outer_fold(ctx, Xs, cfg, fr[0], 0, 0)
        for x in rows:
            print(f"  {x['model']:12s} {x['status'][:40]:40s} {x['seconds']:6.1f}s")
        print("  M3:", {k: v for k, v in info["M3_enet_cox"].items() if k != "selected"})
        print(f"benchmark one outer fold: {time.perf_counter() - t0:.1f}s")
        return
    OUT.mkdir(parents=True, exist_ok=True)
    folds.to_csv(OUT / "folds.csv", index=False)
    t0 = time.perf_counter()
    jobs = [(r, f, None) for r in range(reps) for f in range(cfg["outer_folds"])]
    jobs += [(i % reps, f, C.seed(cfg["seed_base"], 2, 5000 + i)) for i in range(a.null) for f in range(cfg["outer_folds"])]
    results = Parallel(n_jobs=a.jobs, verbose=5)(
        delayed(outer_fold)(ctx, Xs, cfg, fr[r], r, f, s, ("M1", "M1_breslow", "M2_strict", "M3_enet_cox") if s is None else ("M3_enet_cox",))
        for r, f, s in jobs)
    rows = [x for res in results for x in res[0]]
    infos = [res[1] for res in results]
    fold_df = pd.DataFrame(rows)
    real = fold_df[fold_df.null_seed.isna()]
    null = fold_df[fold_df.null_seed.notna()]
    per_rep, summ, deltas = summarise(real)
    n_fits = int(((real.model == "M3_enet_cox") & (real.status == "ok")).sum())
    sel = {}
    for inf in infos:
        if inf["null_seed"] is None and inf.get("M3_enet_cox", {}).get("selected") is not None:
            for g, c in inf["M3_enet_cox"]["selected"].items():
                sel.setdefault(g, []).append(c)
    sel_tab = pd.DataFrame([{"gene": g, "selected_in": len(v), "of_fits": n_fits, "frequency": len(v) / n_fits,
                             "positive_sign_fraction": float(np.mean(np.asarray(v) > 0)), "median_coef": float(np.median(v))}
                            for g, v in sel.items()])
    if len(sel_tab):
        sel_tab = sel_tab.sort_values(["selected_in", "gene"], ascending=[False, True])
    null_rep = []
    for s, g in null.groupby("null_seed"):
        ok = g[g.status == "ok"]
        null_rep.append({"null_seed": int(s), "c_within_fold_pooled": ok.c_num.sum() / ok.c_den.sum(), "folds_ok": len(ok)})
    null_rep = pd.DataFrame(null_rep)
    m3 = [inf["M3_enet_cox"] for inf in infos if inf["null_seed"] is None and inf.get("M3_enet_cox")]
    m2 = [inf["M2_strict"] for inf in infos if inf["null_seed"] is None and inf.get("M2_strict")]
    fold_df.to_csv(OUT / "per_fold.csv", index=False)
    per_rep.to_csv(OUT / "per_repeat.csv", index=False)
    sel_tab.to_csv(OUT / "m3_selection_stability.csv", index=False)
    null_rep.to_csv(OUT / "n1_null_replicates.csv", index=False)
    C.write_json(OUT / "fit_details.json", infos)
    nsum = C.summarise(null_rep.c_within_fold_pooled) if len(null_rep) else None
    alarm = None
    if nsum and nsum["n_repeats"] > 1:
        alarm = bool(nsum["mean"] > 0.52 and nsum["mean"] > 0.5 + 2 * nsum["mc_se"])
    C.write_json(OUT / "summary.json", {
        "track": "A: survival ML for cause-specific recurrence, GSE39582 (exploratory, post hoc)",
        "cohort": {"n": len(ctx.ids), "events": int(ctx.e.sum()), "censored": int((1 - ctx.e).sum())},
        "endpoint_caveats": ["rfs.delay time unit undocumented: no fixed-time metrics",
                             "deaths without recurrence (46) censored: cause-specific hazard, not cumulative incidence"],
        "design": {"outer": f"{reps} repeats x {cfg['outer_folds']}-fold, stratified by event", "inner": f"{cfg['inner_folds']}-fold",
                   "folds_sha256": C.frame_sha(folds), "strict_nesting": True, "variance_top_k": cfg["variance_top_k"]},
        "metric_note": "Harrell C (lifelines rules), within-fold pooled; paired deltas use identical comparable pairs. SD/MC SE = partition variability only.",
        "models": summ, "paired_delta_c": deltas,
        "m2_lambda_counts": pd.Series([x["lambda"] for x in m2]).value_counts().to_dict(),
        "m3_tuning": {"l1_ratio_counts": pd.Series([x["l1_ratio"] for x in m3]).value_counts().to_dict(),
                      "alpha_at_grid_edge": int(sum(x["alpha_at_grid_edge"] for x in m3)),
                      "genes_selected_per_fit": C.summarise([x["n_genes_selected"] for x in m3]),
                      "inner_failures": int(sum(len(x["inner_failures"]) for x in m3))},
        "n1_null": {"replicates": int(len(null_rep)), "summary": nsum, "alarm_rule": "mean > 0.52 and > 0.5 + 2*MC SE", "alarm": alarm},
        "failed_outer_fits": real[real.status != "ok"][["model", "repeat", "fold", "status"]].to_dict("records"),
        "phase9_reference_not_comparable": PHASE9,
        "runtime_seconds": round(time.perf_counter() - t0, 1), "environment": C.environment_record(), "git_head": C.git_head(),
    })
    print(f"done in {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
