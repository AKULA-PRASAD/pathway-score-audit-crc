"""Phase 8C: frozen M2 pathway-score analysis.

Audits the behaviour, distribution, relationships and technical sensitivity of
the ALREADY-FROZEN M2 scores. Not a modelling phase.

HARD CONSTRAINTS ENFORCED IN CODE:
  * Outcome blindness at the loading boundary, two layers (file + column).
  * All 453 S3 samples retained. No exclusion, relabelling, cohort revision.
  * Panel membership and roles read from the CURRENT frozen specifications
    (revision 2, amendment M2-A1) and never reconstructed or reordered.
  * sigma_ddof read from the specification, never a library default.
  * Training-fold-only standardisation; held-out samples cannot influence
    training parameters.
  * No normalization, batch correction, CEL acquisition, representation change.
  * Phase 8A/8B artifacts are READ-ONLY and hash-verified before and after.
  * No recurrence association, survival model, feature selection, performance
    evaluation.

Usage:
    /Users/kalki/coad_transportability/.venv/bin/python \
        /Users/kalki/coad_transportability/scripts/34_phase8c_pathway_scores.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import paths                                        # noqa: E402
from src.m2_amendment import resolve_hash                     # noqa: E402
from src.reproducibility import environment_report, set_seeds   # noqa: E402

ANALYSIS = paths.PROCESSED_DIR / "analysis"
FROZEN = paths.PROCESSED_DIR / "frozen"
PANELS = paths.ROOT / "config" / "pathways"
P7C = paths.ROOT / "reports" / "phase7c"
P8A = paths.ROOT / "reports" / "phase8a"
P8B = paths.ROOT / "reports" / "phase8b"
OUT = paths.ROOT / "reports" / "phase8c"
FIG = OUT / "figures"

SEED = 20261102            # the Phase 8A seed; folds are REUSED, not created
N_PERM = 10_000

OUTCOME_COLUMNS = ("rfs.event", "rfs.delay", "os.event", "os.delay_(months)")
BLOCKED_FILES = ("endpoint.csv",)


def _guard_outcome_files() -> None:
    import pandas as _pd
    for fn in ("read_csv", "read_parquet", "read_table"):
        orig = getattr(_pd, fn)

        def blocked(path, *a, _orig=orig, **kw):
            if any(b in str(path) for b in BLOCKED_FILES):
                raise PermissionError(
                    f"PHASE 8C IS OUTCOME-BLIND: refusing to read {path}")
            return _orig(path, *a, **kw)
        setattr(_pd, fn, blocked)


def load_metadata_outcome_blind() -> tuple[pd.DataFrame, list[str]]:
    m = pd.read_csv(paths.SAMPLE_METADATA, keep_default_na=False, dtype=str)
    dropped = [c for c in m.columns if c in OUTCOME_COLUMNS]
    m = m.drop(columns=dropped)
    assert not [c for c in m.columns if c in OUTCOME_COLUMNS]
    return m.set_index("sample_id"), dropped


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def canonical_scorer():
    """The production scorer -- never reimplemented here.

    Since Phase 8C recommendation 4 this is a plain module import from
    `src.pathway_scoring` rather than loading `scripts/32_phase8a_eda.py` by
    path. Same function object, same semantics; the indirection is kept so the
    audit code reads the same and callers need no change.
    """
    from src.pathway_scoring import score_panel
    return score_panel


def score_panel_independent(X: pd.DataFrame, panel: dict, train_ids):
    """Deliberately naive re-implementation of the frozen rule, for audit.

    Explicit python loops, no pandas broadcasting, ddof read from the spec.
    Agreement with the production scorer is evidence neither carries a bug
    that the other would cancel.
    """
    ddof = int(panel["scoring"]["sigma_ddof"])
    tr = list(train_ids)
    out = {}
    for name, gs in panel["gene_sets"].items():
        members = list(dict.fromkeys(gs["members_resolved"]))
        present = [g for g in members if g in X.index]
        kept, mu, sd = [], {}, {}
        for g in present:
            v = X.loc[g, tr].to_numpy(float)
            m = v.sum() / len(v)
            s = float(np.sqrt(((v - m) ** 2).sum() / (len(v) - ddof)))
            if s > 0:
                kept.append(g); mu[g] = m; sd[g] = s
        if len(present) / gs["n_source"] < 0.80 or len(kept) < 20:
            continue
        acc = np.zeros(X.shape[1])
        for g in kept:
            acc += (X.loc[g].to_numpy(float) - mu[g]) / sd[g]
        out[name] = pd.Series(acc / len(kept), index=X.columns)
    return pd.DataFrame(out)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def kruskal_H_from_ranks(ranks, onehot):
    r = np.atleast_2d(ranks); n = r.shape[1]
    ni = onehot.sum(axis=0); Ri = r @ onehot
    return 12.0 / (n * (n + 1)) * (Ri ** 2 / ni).sum(axis=1) - 3 * (n + 1)


def stratified_perm_indices(strata, n_perm, rng):
    n = len(strata); base = np.arange(n)
    out = np.empty((n_perm, n), dtype=np.int64)
    blocks = [np.where(strata == s)[0] for s in np.unique(strata)]
    for i in range(n_perm):
        idx = base.copy()
        for b in blocks:
            idx[b] = rng.permutation(b)
        out[i] = idx
    return out


def phase8a_folds(s3):
    """The folds Phase 8A established. Regenerated, NOT newly chosen."""
    rng = np.random.default_rng(SEED)
    order = np.array(sorted(s3))
    rng.shuffle(order)
    return np.array_split(order, 5)


def section(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)


FROZEN_FILES = [
    FROZEN / "frozen_gene_to_probeset_median.csv",
    FROZEN / "frozen_gene_to_probe_discovery443.csv",
    FROZEN / "frozen_representation.json",
    ANALYSIS / "expression_median_S0.parquet",
    ANALYSIS / "expression_repiqr_disc443_S0.parquet",
    ANALYSIS / "cohort_membership.csv",
    ANALYSIS / "partition_and_batch.csv",
    ANALYSIS / "covariates.csv",
    ANALYSIS / "endpoint.csv",
    PANELS / "m2_panel_4set.json",
    PANELS / "m2_panel_5set.json",
    PANELS / "M2_FREEZE_RECORD.json",
    paths.ROOT / "config" / "analysis" / "A19_FREEZE_RECORD.json",
    paths.ROOT / "config" / "analysis" / "QC_DECISIONS_RECORD.json",
    paths.ROOT / "config" / "analysis" / "QC_SPECIFICATION.json",
    paths.SAMPLE_METADATA,
]
READONLY_PRIOR = sorted(P8A.glob("pathway_score*")) + sorted(P8A.glob("pca_*")) \
    + sorted(P8B.glob("*.csv")) + [P8B / "phase8b_summary.json"]


def main() -> int:
    _guard_outcome_files()
    set_seeds(SEED)
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import stats

    print("=" * 78)
    print("PHASE 8C -- FROZEN M2 PATHWAY-SCORE ANALYSIS")
    print("=" * 78)

    before = {p.name: sha(p) for p in FROZEN_FILES}
    prior_before = {str(p.relative_to(paths.ROOT)): sha(p) for p in READONLY_PRIOR}

    meta, dropped = load_metadata_outcome_blind()
    print(f"  outcome blindness: files {list(BLOCKED_FILES)}; columns dropped {dropped}")

    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    assert len(s3) == 453
    X = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")[s3].astype(np.float64)
    Xr = pd.read_parquet(ANALYSIS / "expression_repiqr_disc443_S0.parquet")[s3].astype(np.float64)
    part = pd.read_csv(ANALYSIS / "partition_and_batch.csv").set_index("sample_id").reindex(s3)
    cov = pd.read_csv(ANALYSIS / "covariates.csv", keep_default_na=False, dtype=str).set_index("sample_id").reindex(s3)
    qcm = pd.read_csv(P7C / "sample_qc_metrics.csv").set_index("sample_id").reindex(s3)
    meta3 = meta.reindex(s3)

    panels = {"4set_primary": json.loads((PANELS / "m2_panel_4set.json").read_text()),
              "5set_sensitivity": json.loads((PANELS / "m2_panel_5set.json").read_text())}
    record = json.loads((PANELS / "M2_FREEZE_RECORD.json").read_text())
    score_panel = canonical_scorer()

    summary: dict = {
        "phase": "8C", "version": "1.0", "cohort": "S3", "n_samples": len(s3),
        "seed": SEED, "samples_excluded": 0, "outcome_columns_opened": False,
        "outcome_columns_dropped_at_read": list(dropped),
        "outcome_files_blocked": list(BLOCKED_FILES),
        "modelling_performed": False, "recurrence_association_tested": False,
        "feature_selection_performed": False,
        "performance_evaluation_performed": False,
        "normalization_applied": "none", "batch_correction_applied": "none",
        "folds_newly_created": False,
        "folds_source": "the Phase 8A outcome-blind 5-way split, regenerated "
                        "from seed 20261102; no new folds were chosen",
    }

    # =====================================================================
    section("A. REPRODUCE AND AUDIT THE FROZEN SCORES")
    # =====================================================================
    audit = {}
    for pname, panel in panels.items():
        spec_role = panel["role"]
        rec_key = "primary" if pname.startswith("4") else "sensitivity"
        rec_sets = record["specifications"][rec_key]["gene_sets"]
        assert sorted(panel["gene_sets"]) == sorted(rec_sets), pname
        assert panel["specification_revision"] == 2
        ddof = panel["scoring"]["sigma_ddof"]
        assert ddof == 1

        sc_can, diag = score_panel(X, panel, train_ids=s3)
        sc_ind = score_panel_independent(X, panel, s3)
        d_impl = float(np.abs(sc_can[sc_ind.columns].values - sc_ind.values).max())

        p8a_file = P8A / f"pathway_scores_{pname}_EXPLORATORY_full_cohort.csv"
        p8a = pd.read_csv(p8a_file, index_col=0).reindex(s3)
        cols = [c for c in p8a.columns if c.startswith("HALLMARK")]
        d_8a = float(np.abs(sc_can[cols].values - p8a[cols].values).max())

        print(f"\n  --- {pname} (rev {panel['specification_revision']}, role {spec_role}, "
              f"ddof={ddof}) ---")
        for nm, d in diag.items():
            print(f"    {nm.replace('HALLMARK_',''):42s} source {d['n_source']:3d} "
                  f"resolved {d['n_resolved_in_panel']:3d} present {d['n_present_in_matrix']:3d} "
                  f"used {d['n_genes_used']:3d}  zero-var {d['n_zero_variance_dropped']} "
                  f"cov {d['frac_of_source_resolved']*100:.1f}%  "
                  f"{'MET' if d['coverage_rule_met'] else 'FAILED'}")
        print(f"    production vs independent re-implementation : max|diff| {d_impl:.3e}")
        print(f"    vs Phase 8A artifact (CSV rounded to 6 dp)  : max|diff| {d_8a:.3e}")
        assert d_impl < 1e-12, d_impl
        assert d_8a < 1e-6, d_8a

        # ddof is genuinely in force, not coincidental
        sc_d0, _ = score_panel(X, panel, train_ids=s3, ddof=0)
        ratio = float(np.nanmax(np.abs(sc_d0.values / sc_can.values)))
        # fold isolation: deletion AND poisoning
        folds = phase8a_folds(s3)
        iso = []
        for fi, test in enumerate(folds, 1):
            train = [x for x in sorted(s3) if x not in set(test)]
            a, _ = score_panel(X, panel, train_ids=train)
            b, _ = score_panel(X.drop(columns=list(test)), panel, train_ids=train)
            poisoned = X.copy(); poisoned.loc[:, list(test)] += 1000.0
            c, _ = score_panel(poisoned, panel, train_ids=train)
            iso.append({"fold": fi, "n_train": len(train), "n_test": len(test),
                        "bit_identical_when_heldout_deleted":
                            bool(np.array_equal(a.loc[train].values, b.loc[train].values)),
                        "bit_identical_when_heldout_poisoned":
                            bool(np.array_equal(a.loc[train].values, c.loc[train].values)),
                        "heldout_scores_do_move_when_poisoned":
                            bool(not np.allclose(a.loc[list(test)].values,
                                                 c.loc[list(test)].values)),
                        "max_abs_shift_vs_full_cohort_fit":
                            float(np.abs(a.loc[list(test)].values
                                         - sc_can.loc[list(test)].values).max())})
        assert all(f["bit_identical_when_heldout_deleted"] for f in iso)
        assert all(f["bit_identical_when_heldout_poisoned"] for f in iso)
        assert all(f["heldout_scores_do_move_when_poisoned"] for f in iso)
        print(f"    fold isolation (Phase 8A folds): 5/5 bit-identical under "
              f"deletion AND +1000 poisoning of held-out samples")

        # zero-variance path exercised on a synthetic constant gene
        Xz = X.copy()
        victim = list(dict.fromkeys(
            panel["gene_sets"][sorted(panel["gene_sets"])[0]]["members_resolved"]))[0]
        Xz.loc[victim] = 7.0
        _, dz = score_panel(Xz, panel, train_ids=s3)
        zdrop = dz[sorted(panel["gene_sets"])[0]]["n_zero_variance_dropped"]
        print(f"    zero-variance path: injecting a constant gene drops {zdrop} "
              f"(observed in real data: 0)")
        assert zdrop == 1

        sc_can.to_csv(OUT / f"scores_{pname}_reproduced_full_cohort.csv")
        audit[pname] = {
            "specification_revision": panel["specification_revision"],
            "role": spec_role, "sigma_ddof": ddof,
            "membership_matches_freeze_record": True,
            "coverage": {k: {kk: vv for kk, vv in v.items() if kk != "zero_variance_genes"}
                         for k, v in diag.items()},
            "max_abs_diff_production_vs_independent": d_impl,
            "max_abs_diff_vs_phase8a_artifact": d_8a,
            "phase8a_artifact_reproduced": True,
            "ddof0_over_ddof1_ratio": ratio,
            "expected_ratio_sqrt_n_over_n_minus_1": float(np.sqrt(len(s3)/(len(s3)-1))),
            "fold_isolation": iso,
            "zero_variance_path_exercised": {"injected": 1, "dropped": int(zdrop),
                                             "observed_in_real_data": 0},
            "per_sample_missingness": 0,
            "non_finite_in_matrix": int((~np.isfinite(X.values)).sum()),
        }
    summary["audit"] = audit
    sc4, _ = score_panel(X, panels["4set_primary"], train_ids=s3)
    sc5, _ = score_panel(X, panels["5set_sensitivity"], train_ids=s3)

    # the four shared scores must be IDENTICAL across panels: each score
    # depends only on its own members, so adding a 5th set cannot move them.
    shared = list(sc4.columns)
    d_panels = float(np.abs(sc4[shared].values - sc5[shared].values).max())
    print(f"\n  4-set vs 5-set, the four shared scores: max|diff| {d_panels:.3e} "
          f"(must be 0 -- a score depends only on its own members)")
    assert d_panels == 0.0
    summary["panel_independence_of_shared_scores"] = {"max_abs_diff": d_panels}

    # =====================================================================
    section("B. SCORE DISTRIBUTIONS AND BEHAVIOUR")
    # =====================================================================
    qs = [0, 1, 5, 25, 50, 75, 95, 99, 100]
    rows = []
    for c in sc5.columns:
        v = sc5[c].values
        k = audit["5set_sensitivity"]["coverage"][c]["n_genes_used"]
        var = float(v.var(ddof=1))
        # Var(mean of k unit-variance z-scores) = (1 + (k-1)*rbar)/k  ->  rbar
        rbar = (k * var - 1.0) / (k - 1)
        r = {"pathway": c, "n_genes_used": k, "mean": float(v.mean()),
             "sd": float(v.std(ddof=1)), "var": var,
             "skew": float(stats.skew(v)), "kurtosis": float(stats.kurtosis(v)),
             "iqr": float(np.percentile(v, 75) - np.percentile(v, 25)),
             "implied_mean_inter_gene_r": float(rbar),
             "sd_if_members_independent": float(1 / np.sqrt(k)),
             "n_missing": int(np.isnan(v).sum())}
        r.update({f"p{q}": float(np.percentile(v, q)) for q in qs})
        rows.append(r)
    dist = pd.DataFrame(rows).set_index("pathway")
    dist.round(6).to_csv(OUT / "score_distributions.csv")
    print(f"  {'pathway':34s} {'k':>4s} {'SD':>7s} {'IQR':>7s} {'skew':>6s} "
          f"{'kurt':>6s} {'SD_indep':>8s} {'implied r̄':>10s}")
    for c, r in dist.iterrows():
        print(f"  {c.replace('HALLMARK_',''):34s} {int(r.n_genes_used):4d} "
              f"{r['sd']:7.4f} {r['iqr']:7.4f} {r['skew']:6.2f} "
              f"{r['kurtosis']:6.2f} {r['sd_if_members_independent']:8.4f} "
              f"{r['implied_mean_inter_gene_r']:10.4f}")

    # leave-one-out influence -- reported, NEVER an exclusion criterion
    infl = []
    for c in sc5.columns:
        v = sc5[c].values
        sd_full = v.std(ddof=1)
        d = np.array([abs(np.delete(v, i).std(ddof=1) - sd_full) for i in range(len(v))])
        j = int(np.argmax(d))
        infl.append({"pathway": c, "max_abs_delta_sd": float(d[j]),
                     "relative_to_sd": float(d[j] / sd_full),
                     "sample_id": s3[j], "sample_score": float(v[j]),
                     "sample_z": float((v[j] - v.mean()) / sd_full)})
    infldf = pd.DataFrame(infl)
    infldf.round(6).to_csv(OUT / "score_influence_loo.csv", index=False)
    print(f"\n  leave-one-out influence on score SD (REPORTED, never an exclusion rule):")
    for _, r in infldf.iterrows():
        print(f"    {r.pathway.replace('HALLMARK_',''):34s} max ΔSD "
              f"{r.max_abs_delta_sd:.5f} ({r.relative_to_sd*100:.2f}% of SD) "
              f"by {r.sample_id} (z={r.sample_z:+.2f})")
    summary["distributions"] = {
        "table": dist.reset_index().to_dict("records"),
        "influence_loo": infldf.to_dict("records"),
        "transformation_vs_data":
            "score SD is a property of the TRANSFORMATION as much as the data: "
            "averaging k unit-variance z-scores gives SD 1/sqrt(k) if members "
            "were independent. The implied mean inter-gene correlation r̄ "
            "inverts Var = (1+(k-1)r̄)/k and separates set coherence from the "
            "shrinkage that averaging necessarily produces.",
        "variation_sufficiency":
            "all five scores have non-zero, non-degenerate dispersion across "
            "453 samples. Whether that dispersion is SUFFICIENT for a model is "
            "an outcome-dependent question and is NOT answered here.",
    }

    # =====================================================================
    section("C. RELATIONSHIPS AMONG PATHWAYS")
    # =====================================================================
    rel = {}
    for pname, sc in (("4set_primary", sc4), ("5set_sensitivity", sc5)):
        cm = sc.corr(); cs = sc.corr(method="spearman")
        cm.round(6).to_csv(OUT / f"score_correlation_pearson_{pname}.csv")
        cs.round(6).to_csv(OUT / f"score_correlation_spearman_{pname}.csv")
        ew = np.linalg.eigvalsh(cm.values)[::-1]
        pv = ew / ew.sum()
        off = cm.values[np.triu_indices(len(cm), k=1)]
        rec6c = panels[pname]["collinearity_all461_descriptive_only"]
        r = {"max_abs_r": float(np.abs(off).max()),
             "median_abs_r": float(np.median(np.abs(off))),
             "n_gt_0.5": int((np.abs(off) > 0.5).sum()),
             "n_gt_0.7": int((np.abs(off) > 0.7).sum()),
             "condition_number": float(ew[0] / ew[-1]),
             "min_eigenvalue": float(ew[-1]),
             "eigenvalues": [float(x) for x in ew],
             "pcs_80": int(np.searchsorted(np.cumsum(pv), 0.80) + 1),
             "kaiser": int((ew > 1).sum()),
             "phase6c_recorded_n461": {k: rec6c[k] for k in
                                       ("max_abs_r", "condition_number",
                                        "min_eigenvalue", "n_gt_0.5", "pcs_80",
                                        "kaiser") if k in rec6c}}
        rel[pname] = r
        print(f"\n  --- {pname} ---")
        print(f"    max|r| {r['max_abs_r']:.3f} (6C n=461: {rec6c['max_abs_r']:.3f})  "
              f"median|r| {r['median_abs_r']:.3f}")
        print(f"    condition {r['condition_number']:.3f} (6C {rec6c['condition_number']:.3f})  "
              f"min eig {r['min_eigenvalue']:.4f} (6C {rec6c['min_eigenvalue']:.4f})")
        print(f"    n(|r|>0.5) {r['n_gt_0.5']} (6C {rec6c['n_gt_0.5']})  "
              f"PCs for 80% {r['pcs_80']} (6C {rec6c['pcs_80']})  Kaiser {r['kaiser']}")

    cm5 = sc5.corr()
    pairs = [(a, b, float(cm5.loc[a, b])) for i, a in enumerate(cm5.columns)
             for b in cm5.columns[i+1:]]
    strongest = max(pairs, key=lambda t: abs(t[2]))
    print(f"\n  strongest 5-set pair: {strongest[0].replace('HALLMARK_','')} <-> "
          f"{strongest[1].replace('HALLMARK_','')}  r={strongest[2]:+.3f}")

    # PI3K incremental DESCRIPTIVE contribution (no outcome anywhere)
    pi = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
    A = np.c_[np.ones(len(sc5)), sc4.values]
    beta, *_ = np.linalg.lstsq(A, sc5[pi].values, rcond=None)
    resid = sc5[pi].values - A @ beta
    r2 = 1 - resid.var(ddof=1) / sc5[pi].values.var(ddof=1)
    print(f"  PI3K regressed on the four primary scores: R² = {r2:.4f} "
          f"-> {(1-r2)*100:.1f}% of its variance is not linearly captured by them")
    rel["pi3k_increment"] = {
        "r2_on_primary_four": float(r2),
        "residual_variance_fraction": float(1 - r2),
        "residual_sd": float(resid.std(ddof=1)),
        "condition_number_4set": rel["4set_primary"]["condition_number"],
        "condition_number_5set": rel["5set_sensitivity"]["condition_number"],
        "min_eigenvalue_4set": rel["4set_primary"]["min_eigenvalue"],
        "min_eigenvalue_5set": rel["5set_sensitivity"]["min_eigenvalue"],
        "interpretation_limit":
            "a DESCRIPTIVE redundancy measure computed among scores only. It "
            "says nothing about recurrence relevance, and is NOT grounds to "
            "promote, remove or reorder any pathway -- the M2-A1 immutability "
            "clause stands.",
    }
    summary["relationships"] = rel

    # =====================================================================
    section("D. TECHNICAL AND BIOLOGICAL STRUCTURE IN THE SCORES")
    # =====================================================================
    anyflag = qcm[["flag_S_S1", "flag_S_S2", "flag_S_S3"]].any(axis=1)
    annot = {
        "combat_batch": part["combat_batch"].astype(str),
        "dataset": part["dataset"].astype(str),
        "combat_mod": part["combat_mod"].astype(str),
        "mmr.status": cov["mmr.status"].replace({"": np.nan}),
        "qc_any_soft_flag": anyflag.map({True: "flagged", False: "unflagged"}),
        "rna.extraction.method": meta3["rna.extraction.method"].astype(str),
        "tumor.location": cov["tumor.location"].replace({"": np.nan}),
        "cit.molecularsubtype": meta3["cit.molecularsubtype"].astype(str),  # CIRCULAR
    }
    # Two different circularity problems, kept apart:
    #   CIRCULAR          -- the variable IS a clustering of this same matrix.
    #   MATRIX_DERIVED    -- the variable is a Tukey rule on distributional
    #                        properties (median, IQR, sum(1-rho)) of this same
    #                        matrix, so a sample with an unusual IQR will have
    #                        unusual z-scores BY CONSTRUCTION. A dispersion
    #                        association here is partly mechanical.
    CIRCULAR = {"cit.molecularsubtype"}
    MATRIX_DERIVED = {"qc_any_soft_flag"}
    rows = []
    for nm, s in annot.items():
        v = s.replace({"N/A": np.nan}).dropna()
        for c in sc5.columns:
            y = sc5.loc[v.index, c].values
            g = [y[v.values == lv] for lv in sorted(v.unique())]
            g = [x for x in g if len(x) >= 2]
            if len(g) < 2:
                continue
            kw = stats.kruskal(*g)
            lev = stats.levene(*g, center="median")
            n, k = sum(len(x) for x in g), len(g)
            rows.append({"annotation": nm, "pathway": c, "circular": nm in CIRCULAR,
                         "matrix_derived": nm in MATRIX_DERIVED,
                         "n": n, "n_levels": k,
                         "min_cell": int(min(len(x) for x in g)),
                         "location_H": kw.statistic, "location_p": kw.pvalue,
                         "location_epsilon_squared": (kw.statistic - k + 1) / (n - k),
                         "dispersion_W": lev.statistic, "dispersion_p": lev.pvalue})
    st = pd.DataFrame(rows)
    st["location_q_BH"] = bh(st.location_p.values)
    st["dispersion_q_BH"] = bh(st.dispersion_p.values)
    st.to_csv(OUT / "score_vs_annotation.csv", index=False)
    nc = st[~st.circular & ~st.matrix_derived]
    md = st[st.matrix_derived]
    print(f"  {len(st)} tests ({len(annot)} variables x 5 scores), location AND "
          f"dispersion, BH-adjusted")
    print(f"\n  non-circular hits with q<0.05 (location):")
    for _, r in nc[nc.location_q_BH < 0.05].nlargest(8, "location_epsilon_squared").iterrows():
        print(f"    {r.annotation:24s} {r.pathway.replace('HALLMARK_',''):34s} "
              f"eps2={r.location_epsilon_squared:+.4f} q={r.location_q_BH:.2e} "
              f"n={r['n']} mincell={r.min_cell}")
    disp_hits = nc[nc.dispersion_q_BH < 0.05]
    print(f"\n  non-circular hits with q<0.05 (dispersion): {len(disp_hits)}")
    for _, r in disp_hits.nlargest(5, "dispersion_W").iterrows():
        print(f"    {r.annotation:24s} {r.pathway.replace('HALLMARK_',''):34s} "
              f"W={r.dispersion_W:.3f} q={r.dispersion_q_BH:.2e}")
    md_disp = md[md.dispersion_q_BH < 0.05]
    print(f"\n  MATRIX-DERIVED variable (partly mechanical, reported separately): "
          f"{len(md_disp)} dispersion hits q<0.05")
    for _, r in md_disp.nlargest(5, "dispersion_W").iterrows():
        print(f"    {r.annotation:24s} {r.pathway.replace('HALLMARK_',''):34s} "
              f"W={r.dispersion_W:.3f} q={r.dispersion_q_BH:.2e}  [mechanical risk]")

    # Phase 8B discipline: anything nominally significant for batch or
    # extraction method gets a composition-conditioned permutation re-test.
    sv = part["combat_mod"].astype(str).values
    strat = []
    for nm in ("combat_batch", "rna.extraction.method", "qc_any_soft_flag"):
        gv = annot[nm].values
        oh = pd.get_dummies(pd.Series(gv)).values.astype(float)
        pidx = stratified_perm_indices(sv, N_PERM, np.random.default_rng(SEED))
        for c in sc5.columns:
            r = stats.rankdata(sc5[c].values)
            obs = float(kruskal_H_from_ranks(r, oh)[0])
            null = kruskal_H_from_ranks(r[pidx], oh)
            # dispersion variant: |deviation from median|
            dv = np.abs(sc5[c].values - np.median(sc5[c].values))
            rd = stats.rankdata(dv)
            obsd = float(kruskal_H_from_ranks(rd, oh)[0])
            nulld = kruskal_H_from_ranks(rd[pidx], oh)
            strat.append({
                "annotation": nm, "pathway": c, "strata": "combat_mod",
                "location_H": obs, "location_null_p95": float(np.percentile(null, 95)),
                "location_p_strat": float((1 + (null >= obs).sum()) / (1 + N_PERM)),
                "dispersion_H": obsd,
                "dispersion_null_p95": float(np.percentile(nulld, 95)),
                "dispersion_p_strat": float((1 + (nulld >= obsd).sum()) / (1 + N_PERM)),
                "n_permutations": N_PERM})
    stdf = pd.DataFrame(strat)
    stdf.to_csv(OUT / "score_vs_annotation_stratified.csv", index=False)
    surv = stdf[(stdf.location_p_strat < 0.05) | (stdf.dispersion_p_strat < 0.05)]
    print(f"\n  stratified by combat_mod ({N_PERM} permutations within strata):")
    print(f"    tests run: {len(stdf)}   surviving at p<0.05: {len(surv)}")
    for _, r in surv.iterrows():
        print(f"    {r.annotation:22s} {r.pathway.replace('HALLMARK_',''):34s} "
              f"loc p={r.location_p_strat:.4f}  disp p={r.dispersion_p_strat:.4f}")
    summary["technical_structure"] = {
        "n_tests": int(len(st)), "variables": sorted(annot),
        "circular_variables": sorted(CIRCULAR),
        "matrix_derived_variables": sorted(MATRIX_DERIVED),
        "circularity_note":
            "cit.molecularsubtype is the depositors' own EXPRESSION-DERIVED "
            "clustering of this same matrix. An M2 score is also a function of "
            "that matrix, so any association between them is near-tautological "
            "and is NOT independent biological validation.",
        "matrix_derived_note":
            "qc_any_soft_flag is a Tukey rule on the sample median, IQR and "
            "mean inter-sample correlation OF THIS SAME MATRIX (Phase 7B "
            "S-S1/S-S2/S-S3). A sample flagged for an unusual IQR has unusual "
            "z-scores by construction, so a DISPERSION association between the "
            "flags and the scores is partly mechanical, not an independent "
            "discovery of array failure. It is reported, and it is explicitly "
            "NOT grounds for exclusion -- Phase 7C D2 retains all 453.",
        "location_hits_q05_noncircular": int((nc.location_q_BH < 0.05).sum()),
        "dispersion_hits_q05_noncircular": int(len(disp_hits)),
        "matrix_derived_dispersion_hits_q05": int(len(md_disp)),
        "matrix_derived_detail": md_disp.to_dict("records"),
        "top_noncircular": nc.nlargest(8, "location_epsilon_squared").to_dict("records"),
        "stratified": stdf.to_dict("records"),
        "n_surviving_stratification": int(len(surv)),
        "limitation":
            "ComBat was applied to this deposit before release, so batch "
            "variance was removed by construction; a null batch result bounds "
            "what remains detectable and is not evidence that batch effects are "
            "absent. Sparse strata and multiplicity make single nominal "
            "p-values weak evidence.",
    }

    # =====================================================================
    section("E. SENSITIVITY AND STABILITY")
    # =====================================================================
    # E1 panel: already shown identical (max diff 0).
    # E2 representation: median vs representative-IQR, same gene index.
    sc4r, diag_r = score_panel(Xr, panels["4set_primary"], train_ids=s3)
    sc5r, _ = score_panel(Xr, panels["5set_sensitivity"], train_ids=s3)
    rep_rows = []
    for c in sc5.columns:
        a, b = sc5[c].values, sc5r[c].values
        rep_rows.append({"pathway": c,
                         "pearson_median_vs_repiqr": float(np.corrcoef(a, b)[0, 1]),
                         "spearman": float(stats.spearmanr(a, b).statistic),
                         "sd_median": float(a.std(ddof=1)),
                         "sd_repiqr": float(b.std(ddof=1)),
                         "mean_abs_difference": float(np.abs(a - b).mean())})
    repdf = pd.DataFrame(rep_rows); repdf.round(6).to_csv(OUT / "representation_sensitivity.csv", index=False)
    print("  median vs representative-IQR representation (same 21,815-gene index):")
    for _, r in repdf.iterrows():
        print(f"    {r.pathway.replace('HALLMARK_',''):34s} r={r.pearson_median_vs_repiqr:.4f} "
              f"rho={r.spearman:.4f}  SD {r.sd_median:.4f} -> {r.sd_repiqr:.4f}")
    cm5r = sc5r.corr()
    d_corr = float(np.abs(cm5.values - cm5r.values).max())
    ew_r = np.linalg.eigvalsh(cm5r.values)[::-1]
    print(f"    correlation structure: max|Δr| {d_corr:.4f}; condition "
          f"{rel['5set_sensitivity']['condition_number']:.3f} -> {ew_r[0]/ew_r[-1]:.3f}")

    # E3 training-fold composition, using the Phase 8A folds
    folds = phase8a_folds(s3)
    fold_rows = []
    for fi, test in enumerate(folds, 1):
        train = [x for x in sorted(s3) if x not in set(test)]
        sf, _ = score_panel(X, panels["5set_sensitivity"], train_ids=train)
        for c in sc5.columns:
            fold_rows.append({
                "fold": fi, "pathway": c, "n_train": len(train), "n_test": len(test),
                "heldout_mean": float(sf.loc[list(test), c].mean()),
                "heldout_sd": float(sf.loc[list(test), c].std(ddof=1)),
                "max_abs_shift_vs_full_cohort": float(
                    np.abs(sf.loc[list(test), c].values
                           - sc5.loc[list(test), c].values).max()),
                "corr_with_full_cohort_fit": float(np.corrcoef(
                    sf.loc[list(test), c].values, sc5.loc[list(test), c].values)[0, 1])})
    folddf = pd.DataFrame(fold_rows)
    folddf.round(6).to_csv(OUT / "fold_composition_sensitivity.csv", index=False)
    agg = folddf.groupby("pathway").agg(
        max_shift=("max_abs_shift_vs_full_cohort", "max"),
        min_corr=("corr_with_full_cohort_fit", "min"))
    print(f"\n  training-fold composition (the Phase 8A 5-way split, reused):")
    for c, r in agg.iterrows():
        print(f"    {c.replace('HALLMARK_',''):34s} max|shift| {r.max_shift:.4f} SD units"
              f"   min corr with full-cohort fit {r.min_corr:.6f}")
    summary["sensitivity"] = {
        "panel_4set_vs_5set": {"shared_scores_max_abs_diff": d_panels,
                               "finding": "identical by construction"},
        "representation": {
            "comparison": "frozen median aggregation vs frozen representative-IQR "
                          "(discovery_443 selection); both are 21,815-gene "
                          "indices over the same 453 samples, so the panels "
                          "resolve identically and only the probe->gene rule differs",
            "per_pathway": repdf.to_dict("records"),
            "max_abs_delta_correlation_matrix": d_corr,
            "condition_number_median": rel["5set_sensitivity"]["condition_number"],
            "condition_number_repiqr": float(ew_r[0] / ew_r[-1]),
            "scope_limit": "a representation sensitivity check only. Differences "
                           "are differences in the probe->gene summarisation, not "
                           "biological differences.",
            "sensitivity_only": True},
        "fold_composition": {
            "folds": "Phase 8A outcome-blind 5-way split, regenerated from seed "
                     "20261102; NOT newly created and NOT tuned",
            "per_fold": folddf.to_dict("records"),
            "summary": {c: {"max_abs_shift": float(r.max_shift),
                            "min_corr_with_full_cohort": float(r.min_corr)}
                        for c, r in agg.iterrows()}},
    }

    # =====================================================================
    section("FIGURES")
    # =====================================================================
    plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 150, "font.size": 8,
                         "axes.grid": True, "grid.alpha": .25,
                         "axes.spines.top": False, "axes.spines.right": False})
    short = [c.replace("HALLMARK_", "")[:18] for c in sc5.columns]

    fig, ax = plt.subplots(1, 3, figsize=(14, 3.7))
    parts = ax[0].violinplot([sc5[c].values for c in sc5.columns], showmedians=True)
    ax[0].set_xticks(range(1, 6), short, rotation=30, ha="right", fontsize=6)
    ax[0].set_ylabel("score (training-fold SD units)")
    ax[0].set_title("M2 score distributions (5-set, full-cohort fit)")
    ax[1].bar(range(5), dist.loc[sc5.columns, "sd"].values, color="#4878a8", label="observed SD")
    ax[1].plot(range(5), dist.loc[sc5.columns, "sd_if_members_independent"].values,
               "o--", color="#c0392b", ms=4, label="SD if members independent")
    ax[1].set_xticks(range(5), short, rotation=30, ha="right", fontsize=6)
    ax[1].set_ylabel("SD"); ax[1].legend(frameon=False, fontsize=6)
    ax[1].set_title("Dispersion: data vs the averaging transform")
    ax[2].bar(range(5), dist.loc[sc5.columns, "implied_mean_inter_gene_r"].values,
              color="#7b5ea8")
    ax[2].set_xticks(range(5), short, rotation=30, ha="right", fontsize=6)
    ax[2].set_ylabel("implied mean inter-gene r"); ax[2].set_title("Set coherence")
    fig.tight_layout(); fig.savefig(FIG / "H1_score_distributions.png"); plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(14, 3.7))
    im = ax[0].imshow(cm5.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax[0].set_xticks(range(5), short, rotation=45, ha="right", fontsize=6)
    ax[0].set_yticks(range(5), short, fontsize=6); ax[0].grid(False)
    for i in range(5):
        for j in range(5):
            ax[0].text(j, i, f"{cm5.values[i,j]:.2f}", ha="center", va="center", fontsize=5.5)
    ax[0].set_title("5-set correlation (Pearson)"); fig.colorbar(im, ax=ax[0], fraction=.046)
    ew4 = np.linalg.eigvalsh(sc4.corr().values)[::-1]
    ax[1].plot(range(1, 5), ew4, "o-", label=f"4-set (cond {ew4[0]/ew4[-1]:.2f})")
    ax[1].plot(range(1, 6), np.linalg.eigvalsh(cm5.values)[::-1], "s-",
               label=f"5-set (cond {np.linalg.eigvalsh(cm5.values)[::-1][0]/np.linalg.eigvalsh(cm5.values)[::-1][-1]:.2f})")
    ax[1].axhline(1, color="#888", ls=":", lw=1)
    ax[1].set_xlabel("component"); ax[1].set_ylabel("eigenvalue")
    ax[1].set_title("Correlation-matrix spectrum"); ax[1].legend(frameon=False, fontsize=6)
    ax[2].scatter(A @ beta, sc5[pi].values, s=12, alpha=.75, color="#2a7f62")
    lim = [min(sc5[pi].min(), (A@beta).min()), max(sc5[pi].max(), (A@beta).max())]
    ax[2].plot(lim, lim, color="#c0392b", lw=1, ls="--")
    ax[2].set_xlabel("PI3K predicted from the 4 primary scores")
    ax[2].set_ylabel("PI3K score"); ax[2].set_title(f"PI3K redundancy (R²={r2:.3f})")
    fig.tight_layout(); fig.savefig(FIG / "H2_relationships.png"); plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(14, 3.7))
    for c, m in zip(sc5.columns, "os^Dv"):
        ax[0].scatter(sc5[c].values, sc5r[c].values, s=8, alpha=.6, marker=m,
                      label=c.replace("HALLMARK_", "")[:16])
    ax[0].set_xlabel("median representation"); ax[0].set_ylabel("representative-IQR")
    ax[0].set_title("Representation sensitivity"); ax[0].legend(frameon=False, fontsize=5.5)
    for c in sc5.columns:
        sub = folddf[folddf.pathway == c]
        ax[1].plot(sub.fold, sub.max_abs_shift_vs_full_cohort, "o-", ms=3,
                   label=c.replace("HALLMARK_", "")[:16])
    ax[1].set_xlabel("fold"); ax[1].set_ylabel("max |shift| vs full-cohort fit (SD)")
    ax[1].set_title("Fold-composition sensitivity"); ax[1].legend(frameon=False, fontsize=5.5)
    bm = part["combat_mod"].values
    w = 0.35
    for off, (mod, col) in zip((-w/2, w/2), (("modcan1", "#c0392b"), ("modcan2", "#4878a8"))):
        ax[2].bar(np.arange(5) + off, [sc5[c].values[bm == mod].mean() for c in sc5.columns],
                  w, color=col, label=f"{mod} (n={(bm==mod).sum()})")
    ax[2].set_xticks(range(5), short, rotation=30, ha="right", fontsize=6)
    ax[2].set_ylabel("mean score"); ax[2].legend(frameon=False, fontsize=6)
    ax[2].set_title("Mean score by combat_mod (interpretation unresolved)")
    fig.tight_layout(); fig.savefig(FIG / "H3_sensitivity.png"); plt.close(fig)
    print(f"  wrote {len(list(FIG.glob('*.png')))} figures")

    # =====================================================================
    section("INTEGRITY")
    # =====================================================================
    after = {p.name: sha(p) for p in FROZEN_FILES}
    changed = [k for k in before if before[k] != after[k]]
    assert not changed, f"FROZEN ARTIFACT MODIFIED: {changed}"
    prior_after = {str(p.relative_to(paths.ROOT)): sha(p) for p in READONLY_PRIOR}
    pchanged = [k for k in prior_before if prior_before[k] != prior_after[k]]
    assert not pchanged, f"PRIOR-PHASE ARTIFACT MODIFIED: {pchanged}"
    print(f"  frozen artifacts unchanged        : {len(FROZEN_FILES)}/{len(FROZEN_FILES)}")
    print(f"  Phase 8A/8B artifacts unchanged   : {len(READONLY_PRIOR)}/{len(READONLY_PRIOR)}")
    for n in ("m2_panel_4set.json", "m2_panel_5set.json"):
        ok, why = resolve_hash(PANELS / n, after[n]); print(f"  amendment chain {n:22s}: {why}")

    summary["frozen_artifact_hashes"] = after
    summary["prior_phase_artifact_hashes"] = prior_after
    summary["frozen_artifacts_modified"] = False
    summary["panel_membership_changed"] = False
    summary["panel_roles_changed"] = False
    summary["environment"] = environment_report()
    (OUT / "phase8c_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(f"\n  outputs -> {OUT.relative_to(paths.ROOT)}/ "
          f"({len(list(OUT.glob('*.csv')))} CSV + 1 JSON + {len(list(FIG.glob('*.png')))} figures)")
    print("\nPHASE 8C COMPLETE -- no modelling, no exclusions, no outcome access.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
