"""Production M2 pathway scorer -- the single implementation.

Extracted verbatim from `scripts/32_phase8a_eda.py` (Phase 8C recommendation 4)
so that Phase 8A, Phase 8C and any later phase import ONE implementation
instead of loading a script by path. The function body below is byte-identical
to the pre-move version; `tests/fixtures/premove_scorer_reference.json` pins
both its source hash and its outputs, and `tests/test_pathway_scoring_move.py`
fails if either drifts.

Nothing about the scoring semantics changed in the move:

  * frozen panel membership and ordering come from the panel dict, untouched;
  * `sigma_ddof` is read from the specification (amendment M2-A1) and an
    unpinned panel is REFUSED rather than defaulted;
  * mu/sigma are estimated on `train_ids` only, which has no default;
  * held-out samples cannot influence training parameters;
  * zero-variance exclusion and the >=80% / >=20-gene coverage rule are as
    frozen;
  * the return shape is unchanged: (DataFrame of scores, dict of diagnostics).

This module does NOT decide score-level standardisation, panel roles, or
representation policy -- those remain open Phase 8C recommendations 1-3 and
are deliberately not implemented here.
"""

from __future__ import annotations

import pandas as pd

__all__ = ["score_panel"]


def score_panel(X: pd.DataFrame, panel: dict, train_ids, *, ddof: int | None = None):
    """Frozen M2 score: mean of training-fold-standardised member genes.

        score_P(i) = mean_{g in P & index, sigma_g>0} (x_gi - mu_g) / sigma_g

    mu_g / sigma_g are estimated on `train_ids` ONLY and applied unchanged to
    every sample. `train_ids` has NO default: a caller must state the
    estimation population explicitly, so a complete-cohort fit can never
    happen by omission.

    sigma uses the ddof PINNED BY THE SPECIFICATION (amendment M2-A1,
    `scoring.sigma_ddof` = 1). It is read from the panel, never defaulted: if
    the panel does not pin it, this raises rather than silently inheriting
    numpy's ddof=0 or pandas' ddof=1. The `ddof` argument exists ONLY so the
    recorded sensitivity comparison can deviate deliberately; production
    scoring must leave it None.
    """
    if train_ids is None:
        raise ValueError("train_ids is mandatory -- frozen rule: "
                         "mu/sigma are estimated on TRAINING SAMPLES ONLY")
    if ddof is None:
        if "sigma_ddof" not in panel.get("scoring", {}):
            raise ValueError(
                "panel does not pin scoring.sigma_ddof -- refusing to fall back "
                "on a library default (amendment M2-A1)")
        ddof = int(panel["scoring"]["sigma_ddof"])
    train_ids = list(train_ids)
    if len(train_ids) < 2:
        raise ValueError("need >= 2 training samples to estimate sigma")

    sc = {}
    diag = {}
    for name, gs in panel["gene_sets"].items():
        members = list(dict.fromkeys(gs["members_resolved"]))
        present = [g for g in members if g in X.index]
        n_source = int(gs["n_source"])

        Xtr = X.loc[present, train_ids]
        mu = Xtr.mean(axis=1)
        sd = Xtr.std(axis=1, ddof=ddof)
        keep = present if len(present) == 0 else list(sd.index[sd.values > 0])
        dropped_zero_var = [g for g in present if g not in set(keep)]

        # frozen minimum-coverage rule
        frac_resolved = len(present) / n_source if n_source else 0.0
        ok = (frac_resolved >= 0.80) and (len(keep) >= 20)

        if ok:
            Z = X.loc[keep].sub(mu[keep], axis=0).div(sd[keep], axis=0)
            sc[name] = Z.mean(axis=0)
        diag[name] = {
            "n_source": n_source,
            "n_resolved_in_panel": len(members),
            "n_present_in_matrix": len(present),
            "n_absent_from_matrix": len(members) - len(present),
            "frac_of_source_resolved": round(frac_resolved, 6),
            "n_zero_variance_dropped": len(dropped_zero_var),
            "zero_variance_genes": dropped_zero_var,
            "n_genes_used": len(keep) if ok else 0,
            "coverage_rule_met": bool(ok),
            "score_dropped_for_this_fold": not ok,
        }
    return pd.DataFrame(sc), diag
