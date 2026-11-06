"""Phase 9 execution engine -- M1, M2 primary, M2 sensitivity.

Implements the frozen Phase 9 v3.1 protocol. Every data-dependent quantity is
estimated inside the training fold and applied unchanged to held-out data.

  * age: restricted cubic spline, 3 knots at the 10/50/90 percentiles of the
    TRAINING fold (numpy.percentile type 7), recomputed per fold; degenerate
    cases fall back to linear age and are RECORDED (age spline specification).
  * pathway scores: src.pathway_scoring.score_panel, gene-level sigma_ddof=1
    from the specification (M2-A1), training fold only.
  * score standardisation: src.score_standardisation (M2-A2), training fold
    only, zero-variance dropped, non-finite -> NonEstimableReplicate.
  * penalty: clinical block UNPENALISED, molecular block penalised (P6B-A1),
    lambda by inner CV inside the training fold.

No fixed-time / absolute-risk quantity is computed anywhere in this module.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

from src.pathway_scoring import score_panel
from src.score_standardisation import (NonEstimableReplicate,
                                       standardise_scores_training_fold)

__all__ = ["rcs_basis", "build_clinical", "fit_cox", "tune_lambda",
           "pipeline_fit", "c_index", "calibration_slope", "ConvergenceFailure"]

CLINICAL_BINARIES = {"sex": ("Male", "Female"), "tnm.stage": ("3.0", "2.0"),
                     "tumor.location": ("proximal", "distal"),
                     "chemotherapy.adjuvant": ("Y", "N")}
AGE_COL = "age.at.diagnosis_(year)"
REFERENCE_AGE = 69            # reporting only; fixed by the spline specification
KNOT_PCTS = (10, 50, 90)
SEP_ABS_BETA, SEP_ABS_SE = 10.0, 10.0
MAX_ITER, TOL = 100, 1e-9


class ConvergenceFailure(RuntimeError):
    def __init__(self, reason, detail=None):
        super().__init__(reason)
        self.reason, self.detail = reason, detail or {}


def rcs_knots(age_train: np.ndarray):
    """3 knots from TRAINING ages; returns (knots, fell_back_to_linear)."""
    k = np.percentile(np.asarray(age_train, float), KNOT_PCTS)   # type 7
    uniq = np.unique(np.round(k, 6))
    return (k, False) if len(uniq) == 3 else (uniq, True)


def rcs_basis(age, knots):
    """Harrell restricted cubic spline basis, 3 knots -> 2 columns."""
    x = np.asarray(age, float)
    t1, t2, t3 = knots
    cu = lambda z: np.maximum(z, 0.0) ** 3
    x2 = (cu(x - t1)
          - cu(x - t2) * (t3 - t1) / (t3 - t2)
          + cu(x - t3) * (t2 - t1) / (t3 - t2)) / (t3 - t1) ** 2
    return np.column_stack([x, x2])


def build_clinical(cov: pd.DataFrame, ids, knots, linear_age: bool):
    """Clinical design for `ids` using TRAINING-derived knots."""
    sub = cov.loc[ids]
    age = pd.to_numeric(sub[AGE_COL]).to_numpy(float)
    cols, names = [], []
    if linear_age:
        cols.append(age[:, None]); names.append("age_linear")
    else:
        cols.append(rcs_basis(age, knots)); names += ["age_rcs1", "age_rcs2"]
    for c, (pos, _ref) in CLINICAL_BINARIES.items():
        cols.append((sub[c].to_numpy() == pos).astype(float)[:, None])
        names.append(f"{c}={pos}")
    return pd.DataFrame(np.hstack(cols), index=pd.Index(list(ids), name="sample_id"),
                        columns=names)


def fit_cox(design: pd.DataFrame, t, e, penalizer: np.ndarray, strata=None):
    """Cox fit with a PER-COVARIATE penalty vector (P6B-A1)."""
    df = design.reset_index(drop=True).copy()
    df["_T"], df["_E"] = np.asarray(t, float), np.asarray(e, float)
    kw = {}
    if strata is not None:
        df["_S"] = list(strata); kw["strata"] = ["_S"]
    cph = CoxPHFitter(penalizer=np.asarray(penalizer, float), l1_ratio=0.0)
    try:
        cph.fit(df, duration_col="_T", event_col="_E", show_progress=False, **kw)
    except Exception as exc:
        raise ConvergenceFailure(f"fit failed: {type(exc).__name__}: {exc}") from exc
    beta = cph.params_.to_numpy(float)
    se = cph.standard_errors_.to_numpy(float)
    if not np.all(np.isfinite(beta)):
        raise ConvergenceFailure("non-finite coefficient")
    if np.nanmax(np.abs(beta)) > SEP_ABS_BETA or np.nanmax(se[np.isfinite(se)] if np.isfinite(se).any() else [0]) > SEP_ABS_SE:
        raise ConvergenceFailure("separation sentinel tripped",
                                 {"max_abs_beta": float(np.nanmax(np.abs(beta))),
                                  "max_se": float(np.nanmax(se))})
    return cph


def linear_predictor(cph: CoxPHFitter, design: pd.DataFrame) -> np.ndarray:
    cols = list(cph.params_.index)
    return design[cols].to_numpy(float) @ cph.params_.to_numpy(float)


def c_index(t, e, lp) -> float:
    """Harrell's C. Higher hazard -> shorter survival, hence -lp."""
    return float(concordance_index(np.asarray(t, float), -np.asarray(lp, float),
                                   np.asarray(e, float)))


def calibration_slope(t, e, lp) -> float:
    """Slope of the LOG-HAZARD scale: Cox fit on the linear predictor alone."""
    lp = np.asarray(lp, float)
    if np.std(lp) < 1e-12:
        return float("nan")
    d = pd.DataFrame({"lp": lp, "_T": np.asarray(t, float), "_E": np.asarray(e, float)})
    cph = CoxPHFitter()
    try:
        cph.fit(d, duration_col="_T", event_col="_E")
    except Exception:
        return float("nan")
    return float(cph.params_["lp"])


def _event_stratified_folds(e, k, rng):
    e = np.asarray(e, float); idx = np.arange(len(e))
    folds = [[] for _ in range(k)]
    for lab in (1.0, 0.0):
        sel = rng.permutation(idx[e == lab])
        for j, i in enumerate(sel):
            folds[j % k].append(i)
    return [np.array(sorted(f)) for f in folds]


def tune_lambda(design, t, e, n_clin, grid, rng, inner_k=5):
    """Inner-CV ridge lambda for the MOLECULAR block only (P6B-A1)."""
    p = design.shape[1]
    folds = _event_stratified_folds(e, inner_k, rng)
    best, best_score = grid[0], -np.inf
    for lam in grid:
        pen = np.r_[np.zeros(n_clin), np.full(p - n_clin, lam)]
        scores = []
        for f in folds:
            tr = np.setdiff1d(np.arange(len(e)), f)
            if np.asarray(e)[tr].sum() < 5 or np.asarray(e)[f].sum() < 2:
                continue
            try:
                m = fit_cox(design.iloc[tr], np.asarray(t)[tr], np.asarray(e)[tr], pen)
                scores.append(c_index(np.asarray(t)[f], np.asarray(e)[f],
                                      linear_predictor(m, design.iloc[f])))
            except ConvergenceFailure:
                scores.append(np.nan)
        s = np.nanmean(scores) if np.any(np.isfinite(scores)) else -np.inf
        if s > best_score:
            best, best_score = lam, s
    return float(best), float(best_score)


def pipeline_fit(X, cov, panel, train_ids, apply_ids, t_train, e_train,
                 lam_grid, rng, strata_train=None):
    """The ENTIRE pipeline, fitted on `train_ids`, applied to `apply_ids`.

    Returns (model, design_train, design_apply, diagnostics). Every estimated
    quantity -- knots, gene-level mu/sigma, score standardisation, lambda --
    comes from `train_ids` only.
    """
    diag = {}
    age_tr = pd.to_numeric(cov.loc[train_ids, AGE_COL]).to_numpy(float)
    knots, linear = rcs_knots(age_tr)
    diag["age_knots"] = [float(k) for k in knots]
    diag["age_linear_fallback"] = bool(linear)

    d_tr = build_clinical(cov, train_ids, knots, linear)
    d_ap = build_clinical(cov, apply_ids, knots, linear)
    n_clin = d_tr.shape[1]
    diag["n_clinical_df"] = n_clin

    if panel is not None:
        sc_all, cov_diag = score_panel(X, panel, train_ids=list(train_ids))
        std_all, std_diag = standardise_scores_training_fold(sc_all, list(train_ids))
        diag["score_coverage"] = {k: {kk: vv for kk, vv in v.items()
                                      if kk != "zero_variance_genes"}
                                  for k, v in cov_diag.items()}
        diag["score_standardisation"] = {"n_dropped": std_diag["n_dropped"],
                                         "dropped": std_diag["dropped_pathways"]}
        d_tr = pd.concat([d_tr, std_all.loc[train_ids].set_axis(d_tr.index)], axis=1)
        d_ap = pd.concat([d_ap, std_all.loc[apply_ids].set_axis(d_ap.index)], axis=1)

    p = d_tr.shape[1]
    if p > n_clin:
        lam, inner_c = tune_lambda(d_tr, t_train, e_train, n_clin, lam_grid, rng)
        diag["lambda"] = lam
        diag["inner_cv_c"] = inner_c
    else:
        lam = 0.0
        diag["lambda"] = 0.0
    pen = np.r_[np.zeros(n_clin), np.full(p - n_clin, lam)]
    diag["penalizer_vector"] = [float(x) for x in pen]
    model = fit_cox(d_tr, t_train, e_train, pen, strata=strata_train)
    return model, d_tr, d_ap, diag
