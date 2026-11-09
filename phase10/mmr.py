"""Track B: MMR-associated expression classification in GSE39582 (exploratory).

GSE39582's deposited ComBat correction protected `combat_mod`, which is concordant with MMR
status (X4). Within-cohort results here are therefore likely biased upward and are labelled so.

Usage:  python -m phase10.mmr [--benchmark] [--repeats N]
"""
from __future__ import annotations

import argparse
import time
import warnings

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from phase10 import common as C

OUT = C.REPORTS / "mmr"


# ---------------------------------------------------------------- data
def load_gse(shared_genes):
    """409 = the 452 Phase 9 cohort ∩ non-missing mmr.status. Expression restricted to shared genes."""
    from src.phase9_corrected.data import load_context
    ctx = load_context(C.ROOT)
    cov = ctx.cov.loc[ctx.ids]
    lab = cov["mmr.status"]
    ids = [i for i in ctx.ids if lab[i] in ("dMMR", "pMMR")]
    y = pd.Series((lab[ids] == "dMMR").astype(int).to_numpy(), index=ids, name="dMMR")
    Xg = ctx.X.loc[list(shared_genes), ids].T.astype(float)
    Xg.index.name = "sample_id"
    c = cov.loc[ids]
    Xc = pd.DataFrame({
        "age": pd.to_numeric(c["age.at.diagnosis_(year)"]),
        "male": (c["sex"] == "Male").astype(float),
        "stage3": (c["tnm.stage"] == "3.0").astype(float),
        "proximal": (c["tumor.location"] == "proximal").astype(float),
    }, index=pd.Index(ids, name="sample_id"))
    return ids, y, Xg, Xc


def shared_gene_list():
    """Gene symbols present in both GSE39582 and TCGA HiSeqV2 (identifiers only, no labels)."""
    import gzip
    with gzip.open(C.TCGA_DIR / "HiSeqV2.gz", "rt") as fh:
        fh.readline()
        tcga = {line.split("\t", 1)[0] for line in fh}
    gse = pd.read_parquet(C.GSE_EXPR, columns=[]).index
    return sorted(set(map(str, gse)) & tcga)


# ---------------------------------------------------------------- models
def build_models(cfg, genes, cv, rs):
    m = cfg["mmr"]
    k = cfg["variance_top_k"]
    enet = GridSearchCV(
        Pipeline([("sel", C.ColumnSelect(genes)), ("var", C.VarianceFilter(k)), ("sc", StandardScaler()),
                  ("clf", LogisticRegression(solver="saga", max_iter=m["enet"]["max_iter"], tol=m["enet"]["tol"],
                                             class_weight=m["class_weight"], random_state=rs))]),
        {"clf__l1_ratio": m["enet"]["l1_ratio"], "clf__C": C.grid(m["enet"]["C"])},
        scoring="roc_auc", cv=cv, n_jobs=-1, refit=True)
    rf = GridSearchCV(
        Pipeline([("sel", C.ColumnSelect(genes)), ("var", C.VarianceFilter(k)), ("sc", StandardScaler()),
                  ("clf", RandomForestClassifier(n_estimators=m["rf"]["n_estimators"], class_weight=m["class_weight"],
                                                 n_jobs=1, random_state=rs))]),
        {"clf__max_features": m["rf"]["max_features"], "clf__min_samples_leaf": m["rf"]["min_samples_leaf"]},
        scoring="roc_auc", cv=cv, n_jobs=-1, refit=True)
    mlh1 = Pipeline([("sel", C.ColumnSelect([m["single_gene"]])), ("sc", StandardScaler()),
                     ("clf", LogisticRegression(C=np.inf, max_iter=1000))])
    clin = Pipeline([("sc", StandardScaler()), ("clf", LogisticRegression(C=np.inf, max_iter=1000))])
    # (estimator, which input frame)
    return {"enet_logistic": (enet, "genes"), "random_forest": (rf, "genes"),
            "mlh1_only": (mlh1, "genes"), "clinical_only": (clin, "clinical")}


def _pipe(est):
    return est.best_estimator_ if hasattr(est, "best_estimator_") else est


def describe_fit(name, est):
    """Selected features and tuning outcome of one fitted model (descriptive)."""
    p = _pipe(est)
    info = {"best_params": {k: (v if not isinstance(v, np.floating) else float(v))
                            for k, v in getattr(est, "best_params_", {}).items()}}
    clf = p.named_steps["clf"]
    cols = p.named_steps["var"].columns_ if "var" in p.named_steps else list(
        p.named_steps["sel"].columns if "sel" in p.named_steps else [])
    if name == "enet_logistic":
        coef = clf.coef_.ravel()
        nz = np.flatnonzero(coef)
        info["selected"] = {cols[i]: float(coef[i]) for i in nz}
        info["n_iter"] = int(np.max(clf.n_iter_))
        info["hit_max_iter"] = bool(np.max(clf.n_iter_) >= clf.max_iter)
    elif name == "random_forest":
        imp = clf.feature_importances_
        top = np.argsort(-imp)[:20]
        info["top20_impurity"] = {cols[i]: float(imp[i]) for i in top}
    elif name == "mlh1_only":
        info["coef_per_sd"] = float(clf.coef_.ravel()[0])
    return info


# ---------------------------------------------------------------- CV
def run_cv(cfg, ids, y, frames, folds, models_factory, repeats=None, only=None, log=print):
    """Outer CV over the saved fold table. Returns per-fold rows, OOF predictions, fit details."""
    base, k_in = cfg["seed_base"], cfg["inner_folds"]
    reps = sorted(folds["repeat"].unique())[: repeats or None]
    fold_rows, oof, details = [], [], []
    for r in reps:
        fr = folds[folds.repeat == r].set_index("sample_id")["fold"]
        for f in sorted(fr.unique()):
            if only and (r, f) != only:
                continue
            te = fr.index[fr == f].tolist()
            tr = fr.index[fr != f].tolist()
            cv = C.inner_cv(base, "mmr", r, f, k_in)
            for name, (est, which) in models_factory(cv, C.seed(base, 1, r, f, 7)).items():
                X = frames[which]
                t0 = time.perf_counter()
                status = "ok"
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        est.fit(X.loc[tr], y.loc[tr])
                    p = est.predict_proba(X.loc[te])[:, 1]
                except Exception as exc:                       # recorded, never hidden
                    status, p = f"failed: {type(exc).__name__}: {exc}", np.full(len(te), np.nan)
                secs = time.perf_counter() - t0
                num, den = C.auc_pair_counts(y.loc[te], p) if status == "ok" else (np.nan, np.nan)
                fold_rows.append({"model": name, "repeat": r, "fold": f, "n_test": len(te),
                                  "n_pos_test": int(y.loc[te].sum()), "auc_num": num, "auc_den": den,
                                  "fold_auroc": num / den if den else np.nan, "status": status,
                                  "seconds": round(secs, 2)})
                oof += [{"model": name, "repeat": r, "sample_id": s, "y": int(y[s]), "p": float(pi)}
                        for s, pi in zip(te, p)]
                if status == "ok":
                    details.append({"model": name, "repeat": r, "fold": f, **describe_fit(name, est)})
                log(f"  r{r} f{f} {name:14s} {status:4s} {secs:7.1f}s AUROC_fold={num / den if den else float('nan'):.3f}")
    return pd.DataFrame(fold_rows), pd.DataFrame(oof), details


def summarise_cv(folds_df, oof, clip):
    """Repeat-level metrics. AUROC: within-fold pooled pairs. AUPRC/Brier/calibration: one OOF
    prediction per patient per repeat (never stacked across repeats)."""
    per_rep = []
    for (m, r), g in folds_df.groupby(["model", "repeat"]):
        ok = g[g.status == "ok"]
        o = oof[(oof.model == m) & (oof.repeat == r)].dropna(subset=["p"])
        row = {"model": m, "repeat": int(r), "n_failed_folds": int((g.status != "ok").sum()),
               "auroc_within_fold_pooled": ok.auc_num.sum() / ok.auc_den.sum()}
        if len(o):
            prev = o.y.mean()
            brier = brier_score_loss(o.y, o.p)
            a, b = C.calibration(o.y, o.p, clip)
            row.update({"auprc": average_precision_score(o.y, o.p), "prevalence": prev, "brier": brier,
                        "scaled_brier": 1 - brier / (prev * (1 - prev)),
                        "calibration_intercept": a, "calibration_slope": b})
        per_rep.append(row)
    per_rep = pd.DataFrame(per_rep)
    # paired differences on identical pairs (same folds, same denominators)
    paired = []
    for ref in ("mlh1_only", "clinical_only"):
        for m in per_rep.model.unique():
            if m == ref:
                continue
            for r in per_rep.repeat.unique():
                a = folds_df[(folds_df.model == m) & (folds_df.repeat == r)].set_index("fold")
                b = folds_df[(folds_df.model == ref) & (folds_df.repeat == r)].set_index("fold")
                both = a.index[(a.status == "ok") & (b.loc[a.index].status == "ok")]
                d = (a.loc[both].auc_num.sum() - b.loc[both].auc_num.sum()) / a.loc[both].auc_den.sum()
                paired.append({"model": m, "reference": ref, "repeat": int(r), "delta_auroc": d})
    paired = pd.DataFrame(paired)
    summary = {}
    for m, g in per_rep.groupby("model"):
        summary[m] = {k: C.summarise(g[k]) for k in g.columns if k not in ("model", "repeat")}
    deltas = {f"{m}_minus_{ref}": C.summarise(g.delta_auroc) for (m, ref), g in paired.groupby(["model", "reference"])}
    return per_rep, paired, summary, deltas


def selection_stability(details, n_fits):
    rows = {}
    for d in details:
        if d["model"] != "enet_logistic":
            continue
        for g, c in d["selected"].items():
            rows.setdefault(g, []).append(c)
    tab = pd.DataFrame([{"gene": g, "selected_in": len(v), "of_fits": n_fits, "frequency": len(v) / n_fits,
                         "positive_sign_fraction": float(np.mean(np.asarray(v) > 0)),
                         "median_coef": float(np.median(v))} for g, v in rows.items()])
    return tab.sort_values(["selected_in", "gene"], ascending=[False, True]) if len(tab) else tab


def rf_top20_stability(details, n_fits):
    rows = {}
    for d in details:
        if d["model"] == "random_forest":
            for g in d["top20_impurity"]:
                rows[g] = rows.get(g, 0) + 1
    tab = pd.DataFrame([{"gene": g, "in_top20": n, "of_fits": n_fits, "frequency": n / n_fits} for g, n in rows.items()])
    return tab.sort_values(["in_top20", "gene"], ascending=[False, True]) if len(tab) else tab


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", action="store_true", help="time one outer fold only; writes nothing")
    ap.add_argument("--repeats", type=int, default=None)
    a = ap.parse_args(argv)
    cfg = C.load_config()
    genes = shared_gene_list()
    ids, y, Xg, Xc = load_gse(genes)
    frames = {"genes": Xg, "clinical": Xc}
    reps = a.repeats or cfg["repeats"]
    folds = C.make_folds(ids, y.to_numpy(), "mmr", reps, cfg["outer_folds"], cfg["seed_base"])
    print(f"Track B cohort: n={len(ids)} dMMR={int(y.sum())} pMMR={int((1 - y).sum())}; shared genes={len(genes)}")

    def factory(cv, rs):
        return build_models(cfg, genes, cv, rs)

    if a.benchmark:
        t0 = time.perf_counter()
        run_cv(cfg, ids, y, frames, folds, factory, only=(0, 0))
        print(f"benchmark one outer fold: {time.perf_counter() - t0:.1f}s for 4 models")
        return

    OUT.mkdir(parents=True, exist_ok=True)
    (C.REPORTS / "shared_genes.txt").write_text("\n".join(genes) + "\n")
    folds.to_csv(OUT / "folds.csv", index=False)
    t0 = time.perf_counter()
    fold_df, oof, details = run_cv(cfg, ids, y, frames, folds, factory)
    per_rep, paired, summary, deltas = summarise_cv(fold_df, oof, cfg["mmr"]["prob_clip"])
    n_fits = int(((fold_df.model == "enet_logistic") & (fold_df.status == "ok")).sum())
    sel = selection_stability(details, n_fits)
    rft = rf_top20_stability(details, int(((fold_df.model == "random_forest") & (fold_df.status == "ok")).sum()))
    fold_df.to_csv(OUT / "per_fold.csv", index=False)
    oof.to_csv(OUT / "oof_predictions.csv", index=False)
    per_rep.to_csv(OUT / "per_repeat.csv", index=False)
    paired.to_csv(OUT / "paired_delta_auroc.csv", index=False)
    sel.to_csv(OUT / "enet_selection_stability.csv", index=False)
    rft.to_csv(OUT / "rf_top20_stability.csv", index=False)
    C.write_json(OUT / "fit_details.json", details)
    C.write_json(OUT / "summary.json", {
        "track": "B: MMR-associated expression classification, GSE39582 (exploratory)",
        "caveat": "Deposited ComBat protected combat_mod (concordant with MMR; X4). Within-cohort metrics are likely biased upward.",
        "label": "mmr.status as deposited (dMMR=1, pMMR=0); assay: see methods",
        "cohort": {"n": len(ids), "dMMR": int(y.sum()), "pMMR": int((1 - y).sum()), "shared_genes": len(genes)},
        "design": {"outer": f"{reps} repeats x {cfg['outer_folds']}-fold stratified", "inner": f"{cfg['inner_folds']}-fold stratified",
                   "folds_sha256": C.frame_sha(folds), "variance_top_k": cfg["variance_top_k"]},
        "metrics_note": "AUROC = within-fold pooled positive-negative pairs. AUPRC, Brier, calibration computed within each repeat (one prediction per patient). SD and MC SE describe partition variability, not population uncertainty.",
        "models": summary, "paired_delta_auroc": deltas,
        "failed_folds": fold_df[fold_df.status != "ok"][["model", "repeat", "fold", "status"]].to_dict("records"),
        "runtime_seconds": round(time.perf_counter() - t0, 1), "environment": C.environment_record(),
        "git_head": C.git_head(),
    })
    print(f"done in {time.perf_counter() - t0:.0f}s")


def fit_final(cfg, genes, ids, y, frames):
    """Final models on all labelled GSE39582 patients (inner-CV tuning) for the TCGA freeze."""
    cv = C.inner_cv(cfg["seed_base"], "mmr", 999, 999, cfg["inner_folds"])
    out = {}
    for name, (est, which) in build_models(cfg, genes, cv, C.seed(cfg["seed_base"], 1, 999)).items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            est.fit(frames[which].loc[ids], y.loc[ids])
        out[name] = (est, which)
    return out


if __name__ == "__main__":
    main()
