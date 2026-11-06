"""Amendment M2-A2: training-fold-only standardisation of M2 pathway scores.

A SECOND standardisation, applied to the SCORES produced by
`src.pathway_scoring.score_panel`. It does not touch, replace or re-implement
the frozen gene-level standardisation (`sigma_ddof = 1`, amendment M2-A1) that
produces those scores -- `src/pathway_scoring.py` is unmodified by this
amendment and its source hash is unchanged.

The rule, verbatim from M2-A2:

  * centre and scale each pathway score to unit variance using the mean and
    SAMPLE standard deviation (ddof = 1) of the TRAINING samples of that fold;
  * apply those exact training-fold parameters, unchanged, to held-out
    observations;
  * never compute a scaling parameter from a held-out observation;
  * training-fold SD == 0 or < 1e-12 -> DROP that score for that fold and
    record the exclusion (never set to 0 or 1, never silently kept);
  * non-finite training-fold mean or SD -> raise NonEstimableReplicate so the
    resampling loop records the replicate and its reason. Never silently
    redraw, impute or substitute;
  * the unstandardised scores are preserved unchanged by the caller as the
    auditable source features -- this function returns a NEW frame and never
    mutates its input.

Coefficient interpretation after standardisation (recorded by the amendment):
a molecular pathway coefficient is a **log cause-specific hazard ratio per one
training-fold standard deviation of that pathway score**.

Clinical predictor scaling is OUT OF SCOPE here and is unchanged by M2-A2; it
is governed by amendment P6B-A1 (clinical block force-included and
unpenalised, hence scale-equivariant).

NOTE: this module fits no model and reads no outcome.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["NonEstimableReplicate", "ZERO_SD_TOL", "standardise_scores_training_fold"]

#: a training-fold SD at or below this is treated as zero variance.
ZERO_SD_TOL = 1e-12


class NonEstimableReplicate(RuntimeError):
    """A training-fold scaling parameter was non-finite.

    Raised, never swallowed: the resampling loop catches it, records the
    replicate index, seed and reason, excludes that replicate from the
    optimism average and counts it. It never triggers a silent re-draw.
    """

    def __init__(self, reason: str, detail: dict | None = None):
        super().__init__(reason)
        self.reason = reason
        self.detail = detail or {}


def standardise_scores_training_fold(scores: pd.DataFrame, train_ids):
    """Standardise pathway scores using TRAINING-FOLD parameters only.

    Parameters
    ----------
    scores
        samples x pathways, as returned by `score_panel` (unstandardised).
    train_ids
        the training samples of this fold. Mandatory and never defaulted, so a
        full-cohort fit cannot happen by omission -- the same discipline the
        frozen scorer applies to its own `train_ids`.

    Returns
    -------
    (standardised, diagnostics)
        `standardised` is a NEW frame; `scores` is never mutated. Dropped
        pathways are absent from its columns. `diagnostics` records, per
        pathway, the training-fold mean and SD actually used and whether the
        pathway was dropped.
    """
    if train_ids is None:
        raise ValueError(
            "train_ids is mandatory -- M2-A2: scaling parameters are estimated "
            "on TRAINING SAMPLES ONLY")
    train_ids = list(train_ids)
    if len(train_ids) < 2:
        raise ValueError("need >= 2 training samples for a ddof=1 SD")
    missing = [s for s in train_ids if s not in scores.index]
    if missing:
        raise ValueError(f"train_ids not present in the score frame: {missing[:5]}")

    train = scores.loc[train_ids]
    out = {}
    diag: dict[str, dict] = {}
    for col in scores.columns:
        v = train[col].to_numpy(dtype=float)
        mu = float(v.mean())
        sd = float(v.std(ddof=1))          # SAMPLE SD, ddof=1, per M2-A2

        if not (np.isfinite(mu) and np.isfinite(sd)):
            raise NonEstimableReplicate(
                f"non-finite training-fold scaling parameter for {col}",
                {"pathway": col, "mean": mu, "sd": sd, "n_train": len(train_ids)})

        if sd <= ZERO_SD_TOL:
            diag[col] = {"dropped": True,
                         "reason": "training-fold SD <= ZERO_SD_TOL",
                         "train_mean": mu, "train_sd": sd,
                         "zero_sd_tol": ZERO_SD_TOL, "n_train": len(train_ids)}
            continue

        out[col] = (scores[col] - mu) / sd
        diag[col] = {"dropped": False, "train_mean": mu, "train_sd": sd,
                     "n_train": len(train_ids)}

    std = pd.DataFrame(out, index=scores.index)
    return std, {
        "per_pathway": diag,
        "n_dropped": sum(1 for d in diag.values() if d["dropped"]),
        "dropped_pathways": [k for k, d in diag.items() if d["dropped"]],
        "ddof": 1,
        "estimated_on": "TRAINING FOLD SAMPLES ONLY",
        "coefficient_interpretation":
            "log cause-specific hazard ratio per 1 training-fold SD of the score",
        "amendment": "M2-A2",
    }
