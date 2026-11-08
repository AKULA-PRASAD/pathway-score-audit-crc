"""Phase 8A: outcome-blind exploratory data analysis of the approved S3 cohort.

Characterises the deposited GSE39582 expression matrix at the frozen gene-level
representation BEFORE any survival modelling.

HARD CONSTRAINTS ENFORCED IN CODE:
  * NO outcome column is opened. data/processed/analysis/endpoint.csv is
    blocked at the I/O layer by _guard_outcome_files() -- any attempt raises.
  * NO survival modelling, NO recurrence association testing.
  * NO feature selection. Gene-level statistics are reported for ALL 21815
    genes; no gene list is filtered, ranked-and-cut, or carried forward.
  * NO sample exclusion. All 453 S3 samples are retained (Phase 7C D2).
  * NO normalization, NO ComBat, NO batch correction. Matrix used as deposited.
  * Frozen Phase 2 / Phase 3 / M2 / A-19 artifacts are READ-ONLY and their
    hashes are verified before and after the run.

Usage:
    .venv/bin/python scripts/32_phase8a_eda.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                      # noqa: E402
from src.reproducibility import environment_report, set_seeds   # noqa: E402

ANALYSIS = paths.PROCESSED_DIR / "analysis"
FROZEN = paths.PROCESSED_DIR / "frozen"
PANELS = paths.ROOT / "config" / "pathways"
OUT = paths.ROOT / "reports" / "phase8a"
FIG = OUT / "figures"

SEED = 20261102

# Files that contain the recurrence endpoint. Phase 8A is outcome-blind.
BLOCKED = ("endpoint.csv",)


def _guard_outcome_files() -> None:
    """Make reading any outcome file a hard error for the rest of this run."""
    import pandas as _pd

    for fn_name in ("read_csv", "read_parquet", "read_table"):
        orig = getattr(_pd, fn_name)

        def blocked(path, *a, _orig=orig, **kw):
            if any(b in str(path) for b in BLOCKED):
                raise PermissionError(
                    f"PHASE 8A IS OUTCOME-BLIND: refusing to read {path}")
            return _orig(path, *a, **kw)

        setattr(_pd, fn_name, blocked)


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


FROZEN_FILES = [
    FROZEN / "frozen_gene_to_probeset_median.csv",
    FROZEN / "frozen_gene_to_probe_discovery443.csv",
    FROZEN / "frozen_representation.json",
    FROZEN / "selection_population_discovery443.txt",
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
]


def section(t: str) -> None:
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ----------------------------------------------------------------------------
# frozen M2 scoring -- implemented EXACTLY as frozen, train_ids is MANDATORY
# ----------------------------------------------------------------------------
# The production scorer lives in src/pathway_scoring.py (Phase 8C rec. 4).
# It is re-exported here under its original name so this script runs exactly
# as before and existing callers that reach for `score_panel` on this module
# keep working. There is no second copy of the implementation.
from src.pathway_scoring import score_panel                  # noqa: E402,F401


def main() -> int:
    _guard_outcome_files()
    set_seeds(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import stats
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    print("=" * 78)
    print("PHASE 8A -- OUTCOME-BLIND EXPLORATORY DATA ANALYSIS")
    print("=" * 78)
    print("  outcome columns: BLOCKED at the I/O layer (endpoint.csv)")
    print("  normalization / batch correction: NONE")
    print("  sample exclusion: NONE")

    hashes_before = {p.name: sha(p) for p in FROZEN_FILES}

    # ---------------- inputs ----------------
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    assert len(s3) == 453, f"expected 453 S3 samples, got {len(s3)}"

    Xall = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")
    X = Xall[s3].astype(np.float64)
    part = pd.read_csv(ANALYSIS / "partition_and_batch.csv").set_index("sample_id").reindex(s3)
    cov = pd.read_csv(ANALYSIS / "covariates.csv", keep_default_na=False,
                      dtype=str).set_index("sample_id").reindex(s3)

    print(f"\n  frozen median representation : {X.shape[0]} genes x {X.shape[1]} samples")
    print(f"  analysis cohort              : S3, n={len(s3)} (all retained)")
    print(f"  non-finite values            : {int((~np.isfinite(X.values)).sum())}")

    summary: dict = {
        "phase": "8A",
        "title": "outcome-blind exploratory data analysis",
        "cohort": "S3",
        "n_samples": len(s3),
        "n_genes": int(X.shape[0]),
        "representation": "frozen median aggregation (Phase 2B)",
        "matrix_used_as_deposited": True,
        "normalization_applied": "none",
        "batch_correction_applied": "none",
        "samples_excluded": 0,
        "outcome_columns_opened": False,
        "survival_modelling_performed": False,
        "recurrence_association_tested": False,
        "feature_selection_performed": False,
        "seed": SEED,
    }

    V = X.values

    # =====================================================================
    section("1. EXPRESSION-LEVEL ANALYSIS")
    # =====================================================================
    # --- 1a. per-sample distributions ---
    qs = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100]
    Q = np.percentile(V, qs, axis=0)                      # len(qs) x n_samples
    samp = pd.DataFrame(Q.T, index=s3, columns=[f"p{q}" for q in qs])
    samp["mean"] = V.mean(axis=0)
    samp["sd"] = V.std(axis=0, ddof=1)
    samp["iqr"] = samp["p75"] - samp["p25"]
    samp["skew"] = stats.skew(V, axis=0)
    samp["kurtosis"] = stats.kurtosis(V, axis=0)
    samp["dataset"] = part["dataset"].values
    samp["combat_batch"] = part["combat_batch"].values
    samp["combat_mod"] = part["combat_mod"].values
    samp.index.name = "sample_id"
    samp.round(6).to_csv(OUT / "sample_distribution_summary.csv")

    print(f"  per-sample median : {samp['p50'].min():.3f} .. {samp['p50'].max():.3f} "
          f"(spread {samp['p50'].max()-samp['p50'].min():.3f})")
    print(f"  per-sample IQR    : {samp['iqr'].min():.3f} .. {samp['iqr'].max():.3f}")
    print(f"  per-sample mean   : {samp['mean'].mean():.3f} +/- {samp['mean'].std(ddof=1):.3f}")

    summary["expression"] = {
        "sample_median": {"min": float(samp.p50.min()), "max": float(samp.p50.max()),
                          "sd_across_samples": float(samp.p50.std(ddof=1))},
        "sample_iqr": {"min": float(samp.iqr.min()), "max": float(samp.iqr.max()),
                       "sd_across_samples": float(samp.iqr.std(ddof=1))},
        "global": {"min": float(V.min()), "max": float(V.max()),
                   "mean": float(V.mean()), "sd": float(V.std(ddof=1))},
    }

    # --- 1b. gene-level variability (ALL genes; no filtering) ---
    gm = V.mean(axis=1)
    gsd = V.std(axis=1, ddof=1)
    gq = np.percentile(V, [25, 50, 75], axis=1)
    gene = pd.DataFrame({
        "gene_symbol": X.index, "mean": gm, "sd": gsd, "median": gq[1],
        "iqr": gq[2] - gq[0],
        "mad": np.median(np.abs(V - gq[1][:, None]), axis=1),
        "cv": gsd / np.where(np.abs(gm) < 1e-12, np.nan, gm),
        "min": V.min(axis=1), "max": V.max(axis=1),
    }).set_index("gene_symbol")
    gene.round(6).to_csv(OUT / "gene_variability.csv")

    sd_dec = np.percentile(gsd, np.arange(0, 101, 10))
    print(f"  gene SD deciles   : {', '.join(f'{x:.3f}' for x in sd_dec)}")
    print(f"  genes with SD=0   : {int((gsd == 0).sum())}")
    # mean-variance relationship (heteroscedasticity on the log2 scale)
    rho_mean_sd = stats.spearmanr(gm, gsd).statistic
    print(f"  Spearman(mean,SD) : {rho_mean_sd:+.4f}")

    summary["gene_variability"] = {
        "n_genes": int(len(gene)),
        "n_zero_variance": int((gsd == 0).sum()),
        "sd_deciles": [float(x) for x in sd_dec],
        "sd_median": float(np.median(gsd)),
        "spearman_mean_vs_sd": float(rho_mean_sd),
        "note": "reported for ALL genes; no gene is filtered or carried forward",
    }

    # bimodality of gene means (expressed / not-expressed structure)
    hist, edges = np.histogram(gm, bins=80)
    summary["gene_mean_distribution"] = {
        "deciles": [float(x) for x in np.percentile(gm, np.arange(0, 101, 10))],
        "histogram_counts": hist.tolist(),
        "histogram_edges": [float(x) for x in edges],
    }

    # --- 1c. sample-to-sample correlation ---
    def corr_matrix(M: np.ndarray) -> np.ndarray:
        C = M - M.mean(axis=0, keepdims=True)
        C /= np.sqrt((C ** 2).sum(axis=0, keepdims=True))
        return C.T @ C

    P = corr_matrix(V)                                   # Pearson, sample x sample
    R = corr_matrix(np.apply_along_axis(stats.rankdata, 0, V))   # Spearman
    iu = np.triu_indices(len(s3), k=1)
    for nm, M in (("pearson", P), ("spearman", R)):
        v = M[iu]
        print(f"  pairwise {nm:9s}: min {v.min():.4f}  q1 {np.percentile(v,25):.4f}  "
              f"median {np.median(v):.4f}  max {v.max():.4f}")
        summary.setdefault("sample_correlation", {})[nm] = {
            "min": float(v.min()), "p1": float(np.percentile(v, 1)),
            "q1": float(np.percentile(v, 25)), "median": float(np.median(v)),
            "q3": float(np.percentile(v, 75)), "max": float(v.max()),
            "mean": float(v.mean()),
        }
    pd.DataFrame(P, index=s3, columns=s3).round(6).to_csv(
        OUT / "sample_correlation_pearson.csv")

    mean_r = pd.Series((P.sum(axis=1) - 1.0) / (len(s3) - 1), index=s3,
                       name="mean_pearson_r_to_others")
    mean_r.round(6).to_csv(OUT / "sample_mean_correlation.csv")

    # --- 1d. distributional heterogeneity ---
    # Each sample's quantile function vs the cohort-median quantile function.
    grid = np.arange(1, 100)
    QF = np.percentile(V, grid, axis=0)                   # 99 x n_samples
    ref = np.median(QF, axis=1)
    dev = np.abs(QF - ref[:, None])
    het = pd.DataFrame({
        "max_abs_quantile_deviation": dev.max(axis=0),
        "mean_abs_quantile_deviation": dev.mean(axis=0),
        "area_between_quantile_curves": dev.sum(axis=0) / len(grid),
    }, index=s3)
    het.index.name = "sample_id"
    het["dataset"] = part["dataset"].values
    het["combat_batch"] = part["combat_batch"].values
    het.round(6).to_csv(OUT / "distributional_heterogeneity.csv")
    print(f"  max |quantile deviation| from cohort-median curve: "
          f"{het.max_abs_quantile_deviation.min():.3f} .. "
          f"{het.max_abs_quantile_deviation.max():.3f}")
    summary["distributional_heterogeneity"] = {
        "statistic": "sup-norm between a sample's quantile function (percentiles 1-99) "
                     "and the cohort-median quantile function, log2 units",
        "min": float(het.max_abs_quantile_deviation.min()),
        "median": float(het.max_abs_quantile_deviation.median()),
        "max": float(het.max_abs_quantile_deviation.max()),
        "sd_across_quantile_grid_of_cohort": float(np.std(QF, axis=1, ddof=1).mean()),
    }

    # --- 1e. positive control: sex-linked expression vs recorded sex -------
    # Outcome-blind identity check. XIST (X-inactivation) and Y-linked genes
    # are the standard array-level sex control.
    #
    # The two axes are classified INDEPENDENTLY rather than by a single
    # Y-minus-XIST margin, because the two common tumour-biological causes of
    # apparent discordance -- loss of Y (LOY) in male tumours and XIST
    # silencing in female tumours -- move only ONE axis, while a sample or
    # annotation mix-up moves BOTH. A margin rule cannot tell them apart.
    #
    # Thresholds: midpoint between the two 1-D 2-means centres of each axis,
    # computed from the data and reported. Descriptive -- REPORT ONLY, never
    # triggers exclusion (Phase 7C D2 retains all 453 samples).
    ymark = [g for g in ("RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY", "USP9Y", "UTY")
             if g in X.index]
    sexchk = None
    if "XIST" in X.index and ymark:
        def split_1d(v):
            lo, hi = v.min(), v.max()
            c = np.array([lo + (hi - lo) / 4, hi - (hi - lo) / 4])
            for _ in range(200):
                lab = np.abs(v[:, None] - c[None, :]).argmin(axis=1)
                nc = np.array([v[lab == k].mean() if (lab == k).any() else c[k]
                               for k in (0, 1)])
                if np.allclose(nc, c):
                    break
                c = nc
            return float(c.mean()), [float(x) for x in c]

        yscore = X.loc[ymark].mean(axis=0).values
        xist = X.loc["XIST"].values
        ty, cy = split_1d(yscore)
        tx, cx = split_1d(xist)
        y_present = yscore > ty
        xist_high = xist > tx
        rec = cov["sex"].replace({"": np.nan})
        ok = rec.notna().values
        exp_y = (rec.values == "Male")           # expected Y present
        exp_x = (rec.values == "Female")         # expected XIST high

        both = ok & (y_present != exp_y) & (xist_high != exp_x)
        yonly = ok & (y_present != exp_y) & (xist_high == exp_x)
        xonly = ok & (y_present == exp_y) & (xist_high != exp_x)
        conc = ok & (y_present == exp_y) & (xist_high == exp_x)
        ids = np.asarray(s3)

        sexchk = {
            "y_markers_used": ymark,
            "thresholds": {"y_score": ty, "y_cluster_centres": cy,
                           "xist": tx, "xist_cluster_centres": cx,
                           "derivation": "midpoint of 1-D 2-means centres, "
                                         "deterministic init at the quartiles of the range"},
            "n_with_recorded_sex": int(ok.sum()),
            "n_concordant_both_axes": int(conc.sum()),
            "concordance_both_axes": round(float(conc.sum()) / max(int(ok.sum()), 1), 6),
            "discordant_both_axes": sorted(ids[both].tolist()),
            "discordant_y_axis_only": sorted(ids[yonly].tolist()),
            "discordant_xist_axis_only": sorted(ids[xonly].tolist()),
            "interpretation": {
                "both_axes": "female-like on BOTH axes when Male is recorded (or the "
                             "reverse) -- the pattern a sample or annotation mix-up "
                             "would produce; also reachable by LOY plus XIST "
                             "re-expression in the same tumour",
                "y_axis_only": "consistent with loss of Y (LOY), common in male "
                               "solid tumours and in ageing male tissue",
                "xist_axis_only": "consistent with XIST silencing / loss of the "
                                  "inactive X, documented in several carcinomas",
            },
            "adjudication": "NOT ADJUDICATED. Array expression alone cannot separate "
                            "a mix-up from LOY / XIST loss; genotype-level data would "
                            "be required. No sample is excluded or relabelled.",
            "caveat": "the deposited matrix is ComBat-corrected with a 2-level mod "
                      "covariate that does not encode sex, so any confounding between "
                      "sex and ComBat batch was not protected during correction",
            "action": "REPORT ONLY -- never triggers exclusion",
        }
        detail = pd.DataFrame({
            "sample_id": ids, "y_score": yscore, "xist": xist,
            "y_present": y_present, "xist_high": xist_high,
            "recorded_sex": rec.values,
            "combat_batch": part["combat_batch"].values,
        })
        detail["classification"] = np.select(
            [both, yonly, xonly, conc],
            ["discordant_both_axes", "discordant_y_only",
             "discordant_xist_only", "concordant"], default="sex_not_recorded")
        detail.round(6).to_csv(OUT / "sex_positive_control.csv", index=False)
        print(f"  sex control (XIST vs {'/'.join(ymark)}; thresholds "
              f"Y>{ty:.2f}, XIST>{tx:.2f}):")
        print(f"    concordant on both axes : {int(conc.sum())}/{int(ok.sum())} "
              f"({conc.sum()/max(ok.sum(),1)*100:.1f}%)")
        print(f"    discordant on BOTH axes : {int(both.sum())} {sorted(ids[both].tolist())}")
        print(f"    Y axis only (LOY-like)  : {int(yonly.sum())} {sorted(ids[yonly].tolist())}")
        print(f"    XIST axis only          : {int(xonly.sum())} {sorted(ids[xonly].tolist())}")
        print("    NOT ADJUDICATED; no sample excluded or relabelled")
    summary["sex_positive_control"] = sexchk

    # =====================================================================
    section("2. PCA AND SAMPLE STRUCTURE")
    # =====================================================================
    # Documented method: gene-wise centring, NO scaling (matches QC D-S1
    # convention), thin SVD of the centred genes x samples matrix. Fully
    # deterministic -- no randomised solver, no RNG. Sign convention: the
    # loading of largest magnitude is forced positive.
    Xc = V - V.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    for j in range(U.shape[1]):
        k = np.argmax(np.abs(U[:, j]))
        if U[k, j] < 0:
            U[:, j] *= -1
            Vt[j, :] *= -1
    ev = S ** 2
    varexp = ev / ev.sum()
    npc = 20
    scores = pd.DataFrame((Vt[:npc].T * S[:npc]), index=s3,
                          columns=[f"PC{i}" for i in range(1, npc + 1)])
    scores.index.name = "sample_id"
    out_sc = scores.copy()
    out_sc["dataset"] = part["dataset"].values
    out_sc["combat_batch"] = part["combat_batch"].values
    out_sc["combat_mod"] = part["combat_mod"].values
    out_sc.round(6).to_csv(OUT / "pca_scores_gene_level.csv")

    pd.DataFrame({
        "pc": np.arange(1, npc + 1),
        "eigenvalue": ev[:npc],
        "variance_explained": varexp[:npc],
        "cumulative": np.cumsum(varexp)[:npc],
    }).round(8).to_csv(OUT / "pca_variance_explained.csv", index=False)

    print("  variance explained (gene-level, centred not scaled):")
    print("   " + ", ".join(f"PC{i+1} {varexp[i]*100:.2f}%" for i in range(5)))
    print(f"   cumulative PC1-5 {varexp[:5].sum()*100:.2f}%  "
          f"PC1-10 {varexp[:10].sum()*100:.2f}%  PC1-20 {varexp[:20].sum()*100:.2f}%")
    n80 = int(np.searchsorted(np.cumsum(varexp), 0.80) + 1)
    n50 = int(np.searchsorted(np.cumsum(varexp), 0.50) + 1)
    print(f"   PCs for 50% / 80% of variance: {n50} / {n80}")

    # scaled sensitivity PCA
    sd_g = V.std(axis=1, ddof=1)
    nz = sd_g > 0
    Xs = (V[nz] - V[nz].mean(axis=1, keepdims=True)) / sd_g[nz][:, None]
    Ss = np.linalg.svd(Xs, full_matrices=False, compute_uv=False)
    ves = (Ss ** 2) / (Ss ** 2).sum()
    print("  sensitivity (gene-standardised): " +
          ", ".join(f"PC{i+1} {ves[i]*100:.2f}%" for i in range(5)))

    summary["pca"] = {
        "method": "gene-wise centring, no scaling; thin SVD (numpy.linalg.svd, "
                  "full_matrices=False) of the 21815 x 453 centred matrix; "
                  "deterministic, no RNG; sign fixed by largest-|loading| positive",
        "variance_explained_pc1_20": [float(x) for x in varexp[:20]],
        "cumulative_pc5": float(varexp[:5].sum()),
        "cumulative_pc10": float(varexp[:10].sum()),
        "cumulative_pc20": float(varexp[:20].sum()),
        "n_pcs_for_50pct": n50,
        "n_pcs_for_80pct": n80,
        "sensitivity_gene_standardised_pc1_5": [float(x) for x in ves[:5]],
        "qc_phase7c_probe_level_pc1_5": [0.09726447463666461, 0.08085061104061134,
                                         0.06755747764633219, 0.04953452972118323,
                                         0.0381817276108352],
        "full_cohort_rotation_is_descriptive_only": True,
        "note": "this full-cohort rotation is EXPLORATORY. Phase 6B s3.2.2 requires "
                "any PCA rotation used as a model input to be re-estimated inside "
                "each training fold. It is not a model input.",
    }

    # --- 2b. association of PCs with technical and biological annotation ---
    ann_cat = {
        "combat_batch": part["combat_batch"], "dataset": part["dataset"],
        "combat_mod": part["combat_mod"],
        "sex": cov["sex"], "tnm.stage": cov["tnm.stage"],
        "tumor.location": cov["tumor.location"],
        "chemotherapy.adjuvant": cov["chemotherapy.adjuvant"],
        "mmr.status": cov["mmr.status"], "cimp.status": cov["cimp.status"],
        "cin.status": cov["cin.status"], "kras.mutation": cov["kras.mutation"],
        "braf.mutation": cov["braf.mutation"], "tp53.mutation": cov["tp53.mutation"],
    }
    rows = []
    for nm, s in ann_cat.items():
        v = s.astype(str).replace({"": np.nan}).dropna()
        for i in range(1, 6):
            y = scores.loc[v.index, f"PC{i}"].values
            groups = [y[v.values == g] for g in sorted(v.unique())]
            groups = [g for g in groups if len(g) >= 2]
            if len(groups) < 2:
                continue
            H = stats.kruskal(*groups).statistic
            p = stats.kruskal(*groups).pvalue
            n, k = len(v), len(groups)
            eps2 = (H - k + 1) / (n - k) if n > k else np.nan
            rows.append({"annotation": nm, "type": "categorical", "pc": f"PC{i}",
                         "n_available": n, "n_levels": k, "statistic_H": H,
                         "p_value": p, "epsilon_squared": eps2})
    age = pd.to_numeric(cov["age.at.diagnosis_(year)"].replace({"": np.nan}),
                        errors="coerce").dropna()
    for i in range(1, 6):
        r = stats.spearmanr(age.values, scores.loc[age.index, f"PC{i}"].values)
        rows.append({"annotation": "age.at.diagnosis", "type": "continuous",
                     "pc": f"PC{i}", "n_available": len(age), "n_levels": np.nan,
                     "statistic_H": np.nan, "p_value": r.pvalue,
                     "epsilon_squared": np.nan, "spearman_rho": r.statistic})
    assoc = pd.DataFrame(rows)
    assoc.to_csv(OUT / "pc_annotation_associations.csv", index=False)
    top = assoc.dropna(subset=["epsilon_squared"]).nlargest(8, "epsilon_squared")
    print("\n  strongest PC/annotation associations (effect size eps^2, descriptive):")
    for _, r in top.iterrows():
        print(f"    {r.annotation:22s} {r.pc}  eps^2={r.epsilon_squared:6.3f}  "
              f"p={r.p_value:.2e}  (n={r.n_available}, {int(r.n_levels)} levels)")

    summary["pc_associations_top"] = [
        {"annotation": r.annotation, "pc": r.pc,
         "epsilon_squared": float(r.epsilon_squared), "p_value": float(r.p_value),
         "n": int(r.n_available)} for _, r in top.iterrows()]

    # --- 2c. clustering ---
    D = 1.0 - P
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2
    Z = linkage(squareform(D, checks=False), method="average")
    clus = {}
    for k in (2, 3, 4, 5, 6):
        lab = fcluster(Z, t=k, criterion="maxclust")
        sizes = np.bincount(lab)[1:].tolist()
        clus[k] = sizes
        print(f"  hierarchical k={k}: sizes {sizes}")
    lab4 = fcluster(Z, t=4, criterion="maxclust")
    ct = pd.crosstab(pd.Series(lab4, index=s3, name="cluster_k4"),
                     part["combat_batch"])
    ct.to_csv(OUT / "cluster_by_batch.csv")
    summary["clustering"] = {
        "method": "average linkage on 1 - Pearson sample correlation",
        "cluster_sizes": {str(k): v for k, v in clus.items()},
    }

    bt = part["combat_batch"]
    batch_pc = scores[[f"PC{i}" for i in range(1, 6)]].groupby(bt.values).agg(["mean", "std", "count"])
    batch_pc.round(4).to_csv(OUT / "pc_by_batch.csv")

    # --- 2d. technical vs biological structure -----------------------------
    # The deposit's ComBat model covariate `combat_mod` (modcan1/modcan2) was
    # recorded by Phase 7A as UNDOCUMENTED (open question Q1). Test it against
    # the deposited annotation, then ask whether the batch-associated PC
    # structure survives conditioning on it.
    known = cov["mmr.status"].replace({"": np.nan}).dropna()
    ct_mod = pd.crosstab(part.loc[known.index, "combat_mod"], known)
    ct_mod.to_csv(OUT / "combat_mod_vs_mmr.csv")
    off_diag = int(ct_mod.values.sum() - np.trace(
        ct_mod.reindex(index=["modcan1", "modcan2"],
                       columns=["dMMR", "pMMR"]).values))
    print(f"\n  combat_mod vs deposited MMR status (S3, n={len(known)} with MMR recorded):")
    print("   " + ct_mod.to_string().replace("\n", "\n   "))
    print(f"   discordant samples: {off_diag}")

    # same test over every tumour array in the SDRF, not just S3
    full_conc = None
    try:
        from src.external_annotation import load_batch_assignments
        sdrf = load_batch_assignments()
        meta = pd.read_csv(paths.SAMPLE_METADATA, keep_default_na=False, dtype=str)
        j = meta.merge(sdrf, on="sample_id", how="inner")
        jk = j[j["mmr.status"].isin(["dMMR", "pMMR"])]
        ct_full = pd.crosstab(jk["combat_mod"], jk["mmr.status"])
        full_conc = {
            "n_arrays_with_mmr": int(len(jk)),
            "table": {k: v for k, v in ct_full.to_dict().items()},
            "discordant": int(len(jk) - np.trace(ct_full.reindex(
                index=["modcan1", "modcan2"], columns=["dMMR", "pMMR"]).values)),
        }
        print(f"   full SDRF ({len(jk)} tumour arrays with MMR recorded): "
              f"{full_conc['discordant']} discordant")
    except Exception as e:                                   # pragma: no cover
        print(f"   full-series check unavailable: {e}")

    # batch composition and the conditional test
    ct_bm = pd.crosstab(part["combat_batch"], part["combat_mod"])
    ct_bm["pct_modcan1"] = (ct_bm["modcan1"] / ct_bm.sum(axis=1) * 100).round(2)
    ct_bm.to_csv(OUT / "batch_by_combat_mod.csv")
    chi = stats.chi2_contingency(pd.crosstab(part["combat_batch"], part["combat_mod"]))
    print(f"\n  batch x combat_mod: chi2={chi.statistic:.1f} df={chi.dof} "
          f"p={chi.pvalue:.3e}  (batch05 {ct_bm.loc['batch05','pct_modcan1']:.1f}% modcan1 "
          f"vs cohort {(part.combat_mod=='modcan1').mean()*100:.1f}%)")

    cond = []
    for strat in ("modcan2", "modcan1"):
        sel = part["combat_mod"] == strat
        bb = part.loc[sel, "combat_batch"]
        for i in range(1, 6):
            y = scores.loc[sel.values, f"PC{i}"].values
            g = [y[bb.values == b] for b in sorted(bb.unique())]
            g = [x for x in g if len(x) >= 2]
            if len(g) < 2:
                continue
            H = stats.kruskal(*g)
            n, k = sum(len(x) for x in g), len(g)
            cond.append({"stratum": strat, "n": n, "pc": f"PC{i}", "n_batches": k,
                         "statistic_H": H.statistic, "p_value": H.pvalue,
                         "epsilon_squared": (H.statistic - k + 1) / (n - k)})
    condf = pd.DataFrame(cond)
    condf.to_csv(OUT / "pc_batch_conditional_on_combat_mod.csv", index=False)
    marg = assoc[(assoc.annotation == "combat_batch")].set_index("pc")
    print("\n  PC ~ batch, marginal vs stratified by combat_mod:")
    print(f"    {'PC':5s} {'marginal eps2':>14s} {'p':>10s} | "
          f"{'pMMR-only eps2':>15s} {'p':>10s}")
    for i in range(1, 6):
        pc = f"PC{i}"
        s = condf[(condf.stratum == "modcan2") & (condf.pc == pc)]
        if len(s):
            print(f"    {pc:5s} {marg.loc[pc,'epsilon_squared']:14.4f} "
                  f"{marg.loc[pc,'p_value']:10.4f} | "
                  f"{s.epsilon_squared.iloc[0]:15.4f} {s.p_value.iloc[0]:10.4f}")

    summary["technical_vs_biological"] = {
        "combat_mod_meaning": {
            "phase7a_status": "UNDOCUMENTED (open question Q1)",
            "s3_table": {k: v for k, v in ct_mod.to_dict().items()},
            "s3_n_with_mmr_recorded": int(len(known)),
            "s3_discordant": off_diag,
            "full_series": full_conc,
            "observation": "combat_mod is perfectly concordant with deposited MMR "
                           "status; zero discordant arrays",
            "inference_status": "HYPOTHESIS -- concordance does not prove the "
                                "depositors' mod matrix WAS mmr.status; it is "
                                "consistent with it and with any variable that "
                                "partitions these arrays identically",
            "not_used_to_impute": "the 43 S3 samples with MMR unrecorded are NOT "
                                  "imputed from combat_mod; A-19 complete-case is frozen",
        },
        "batch_composition": {
            "chi2": float(chi.statistic), "df": int(chi.dof),
            "p_value": float(chi.pvalue),
            "batch05_pct_modcan1": float(ct_bm.loc["batch05", "pct_modcan1"]),
            "cohort_pct_modcan1": float((part.combat_mod == "modcan1").mean() * 100),
        },
        "pc_batch_conditional": condf.to_dict("records"),
        "interpretation": "the marginal PC1/PC3 batch association does not survive "
                          "stratification by combat_mod. This is consistent with "
                          "composition, not residual technical batch structure. It "
                          "does NOT prove absence of batch effect: ComBat was already "
                          "applied to the deposit, so batch variance was removed by "
                          "construction and this analysis cannot recover it.",
    }

    # =====================================================================
    section("3. FROZEN PATHWAY-SCORE ANALYSIS (M2)")
    # =====================================================================
    panels = {
        "4set_primary": json.loads((PANELS / "m2_panel_4set.json").read_text()),
        "5set_sensitivity": json.loads((PANELS / "m2_panel_5set.json").read_text()),
    }
    path_summary = {}

    for pname, panel in panels.items():
        print(f"\n  --- {pname} ({panel['n_sets']} sets, "
              f"{panel['union_unique_genes']} union genes) ---")

        # (a) DESCRIPTIVE full-cohort scores -- EXPLORATORY ONLY
        sc_full, diag = score_panel(X, panel, train_ids=s3)
        for nm, d in diag.items():
            print(f"    {nm:48s} source {d['n_source']:3d} -> used {d['n_genes_used']:3d}"
                  f"  cov {d['frac_of_source_resolved']*100:5.1f}%"
                  f"  zero-var {d['n_zero_variance_dropped']}"
                  f"  rule {'MET' if d['coverage_rule_met'] else 'FAILED -> DROPPED'}")
        o = sc_full.copy()
        o.index.name = "sample_id"
        o["dataset"] = part["dataset"].values
        o["combat_batch"] = part["combat_batch"].values
        o.round(6).to_csv(OUT / f"pathway_scores_{pname}_EXPLORATORY_full_cohort.csv")

        cmat = sc_full.corr(method="pearson")
        cmat.round(6).to_csv(OUT / f"pathway_score_correlation_{pname}.csv")
        off = cmat.values[np.triu_indices(len(cmat), k=1)]
        print(f"    inter-score |r|: max {np.abs(off).max():.3f}  "
              f"median {np.median(np.abs(off)):.3f}  n(|r|>0.5)={int((np.abs(off)>0.5).sum())}")

        ew = np.linalg.eigvalsh(cmat.values)[::-1]
        pv = ew / ew.sum()
        print(f"    score-matrix condition number {ew[0]/ew[-1]:.3f}  "
              f"min eigenvalue {ew[-1]:.4f}  PCs for 80% {int(np.searchsorted(np.cumsum(pv),0.80)+1)}")

        # (b) FOLD-ISOLATION DEMONSTRATION -- outcome-blind illustrative splits.
        # NOT the Phase 6B resampling scheme (which is defined at modelling
        # time and may use outcome information for stratification).
        rng = np.random.default_rng(SEED)
        order = np.array(sorted(s3))
        rng.shuffle(order)
        folds = np.array_split(order, 5)
        fold_recs, maxdiff = [], 0.0
        for fi, test in enumerate(folds, 1):
            train = [x for x in sorted(s3) if x not in set(test)]
            sc_tr, dg = score_panel(X, panel, train_ids=train)
            # leakage invariant: parameters must not depend on held-out rows
            sc_del, _ = score_panel(X.drop(columns=list(test)), panel,
                                    train_ids=train)
            bit_identical = bool(np.array_equal(
                sc_tr.loc[train].values, sc_del.loc[train].values))
            d = float(np.abs(sc_tr.loc[list(test)].values
                             - sc_full.loc[list(test)].values).max())
            maxdiff = max(maxdiff, d)
            fold_recs.append({
                "fold": fi, "n_train": len(train), "n_test": len(test),
                "params_bit_identical_when_heldout_deleted": bit_identical,
                "max_abs_score_shift_vs_full_cohort_fit": d,
                "all_sets_met_coverage_rule": all(v["coverage_rule_met"]
                                                  for v in dg.values()),
                "zero_variance_dropped": {k: v["n_zero_variance_dropped"]
                                          for k, v in dg.items()},
            })
        assert all(f["params_bit_identical_when_heldout_deleted"] for f in fold_recs)
        print(f"    fold isolation: 5/5 folds bit-identical when held-out "
              f"samples are deleted from the input")
        print(f"    max |score shift| held-out, train-fold fit vs full-cohort fit: "
              f"{maxdiff:.4f} SD units")
        pd.DataFrame(fold_recs).to_csv(
            OUT / f"fold_isolation_{pname}.csv", index=False)

        # (c) ddof: pinned to 1 by amendment M2-A1. The ddof=0 run below is a
        # deliberate deviation, kept to quantify what the pre-amendment
        # ambiguity was worth; it is never used as a score.
        pinned = int(panel["scoring"]["sigma_ddof"])
        assert pinned == 1, f"panel pins sigma_ddof={pinned}, expected 1"
        sc_d0, _ = score_panel(X, panel, train_ids=s3, ddof=0)
        ddof_ratio = float(np.nanmax(np.abs(sc_d0.values / np.where(
            np.abs(sc_full.values) < 1e-12, np.nan, sc_full.values))))
        expected = float(np.sqrt(len(s3) / (len(s3) - 1)))

        path_summary[pname] = {
            "n_sets": panel["n_sets"],
            "coverage": {k: {kk: vv for kk, vv in v.items()
                             if kk != "zero_variance_genes"}
                         for k, v in diag.items()},
            "all_coverage_rules_met": all(v["coverage_rule_met"] for v in diag.values()),
            "scores_dropped": [k for k, v in diag.items() if v["score_dropped_for_this_fold"]],
            "full_cohort_scores_are_exploratory_only": True,
            "full_cohort_scores_used_as_model_input": False,
            "inter_score_correlation": {
                "max_abs_r": float(np.abs(off).max()),
                "median_abs_r": float(np.median(np.abs(off))),
                "n_gt_0.5": int((np.abs(off) > 0.5).sum()),
                "n_gt_0.7": int((np.abs(off) > 0.7).sum()),
                "condition_number": float(ew[0] / ew[-1]),
                "min_eigenvalue": float(ew[-1]),
            },
            "score_distribution": {
                c: {"mean": float(sc_full[c].mean()), "sd": float(sc_full[c].std(ddof=1)),
                    "min": float(sc_full[c].min()), "max": float(sc_full[c].max())}
                for c in sc_full.columns},
            "fold_isolation": {
                "splits": "5 outcome-blind deterministic splits, seed 20261102",
                "is_phase6b_resampling_scheme": False,
                "all_folds_bit_identical_when_heldout_deleted": True,
                "max_abs_score_shift_vs_full_cohort_fit": maxdiff,
            },
            "ddof_sensitivity": {
                "freeze_record_pins_ddof": True,
                "pinned_sigma_ddof": pinned,
                "pinned_by_amendment": "M2-A1",
                "read_from_specification_not_library_default": True,
                "ratio_ddof0_over_ddof1_max": ddof_ratio,
                "expected_constant_factor_sqrt_n_over_n_minus_1": expected,
                "effect": "a single constant rescaling of every score within a fold; "
                          "rank order and correlations unchanged, numeric values are not",
                "amendment_resolves_phase8a_U5": True,
            },
        }

    summary["pathways"] = path_summary

    # pathway scores vs PCs / batch (descriptive)
    sc4, _ = score_panel(X, panels["4set_primary"], train_ids=s3)
    rel = []
    for c in sc4.columns:
        for i in range(1, 6):
            r = stats.spearmanr(sc4[c].values, scores[f"PC{i}"].values)
            rel.append({"pathway": c, "pc": f"PC{i}", "spearman_rho": r.statistic,
                        "p_value": r.pvalue})
        g = [sc4[c].values[bt.values == b] for b in sorted(bt.unique())]
        H = stats.kruskal(*g)
        n, k = len(sc4), len(g)
        rel.append({"pathway": c, "pc": "combat_batch", "spearman_rho": np.nan,
                    "p_value": H.pvalue,
                    "epsilon_squared": (H.statistic - k + 1) / (n - k)})
    pd.DataFrame(rel).to_csv(OUT / "pathway_vs_pc_and_batch.csv", index=False)
    print("\n  pathway score vs PC (Spearman, |rho| >= 0.30):")
    for r in rel:
        if not np.isnan(r.get("spearman_rho", np.nan)) and abs(r["spearman_rho"]) >= 0.30:
            print(f"    {r['pathway']:48s} {r['pc']}  rho={r['spearman_rho']:+.3f}")

    # =====================================================================
    section("FIGURES")
    # =====================================================================
    plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 150,
                         "font.size": 8, "axes.grid": True,
                         "grid.alpha": 0.25, "axes.spines.top": False,
                         "axes.spines.right": False})

    # F1 expression distributions
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    o = np.argsort(samp["p50"].values)
    ax[0].fill_between(np.arange(len(s3)), Q[qs.index(5)][o], Q[qs.index(95)][o],
                       alpha=.25, label="p5-p95", color="#4878a8")
    ax[0].fill_between(np.arange(len(s3)), Q[qs.index(25)][o], Q[qs.index(75)][o],
                       alpha=.45, label="IQR", color="#4878a8")
    ax[0].plot(np.arange(len(s3)), Q[qs.index(50)][o], lw=.8, color="#102a43",
               label="median")
    ax[0].set_xlabel("sample (ordered by median)")
    ax[0].set_ylabel("log2 expression")
    ax[0].set_title(f"Per-sample distribution, all {len(s3)} S3 samples")
    ax[0].legend(frameon=False, fontsize=7)
    for b in sorted(bt.unique()):
        sel = bt.values == b
        ax[1].scatter(samp["p50"].values[sel], samp["iqr"].values[sel], s=10,
                      alpha=.8, label=f"{b} (n={sel.sum()})")
    ax[1].set_xlabel("sample median"); ax[1].set_ylabel("sample IQR")
    ax[1].set_title("Median vs IQR by ComBat batch")
    ax[1].legend(frameon=False, fontsize=5.5, ncol=2)
    ax[2].hist(gm, bins=80, color="#4878a8", alpha=.85)
    ax[2].set_xlabel("gene mean (log2)"); ax[2].set_ylabel("genes")
    ax[2].set_title(f"Gene-mean distribution, {len(gene)} genes")
    fig.tight_layout(); fig.savefig(FIG / "F1_expression_distributions.png"); plt.close(fig)

    # F2 gene variability
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    ax[0].hist(gsd, bins=90, color="#7b5ea8", alpha=.85)
    ax[0].set_xlabel("gene SD (log2 units)"); ax[0].set_ylabel("genes")
    ax[0].set_title("Gene-level variability")
    ax[1].hexbin(gm, gsd, gridsize=60, cmap="viridis", bins="log", mincnt=1)
    ax[1].set_xlabel("gene mean"); ax[1].set_ylabel("gene SD")
    ax[1].set_title(f"Mean-variance relationship (Spearman {rho_mean_sd:+.3f})")
    ax[2].plot(np.arange(len(gsd)) / len(gsd) * 100,
               np.sort(gsd)[::-1], lw=1.2, color="#7b5ea8")
    ax[2].set_xlabel("percentile of genes (most variable first)")
    ax[2].set_ylabel("gene SD"); ax[2].set_title("Variability concentration")
    fig.tight_layout(); fig.savefig(FIG / "F2_gene_variability.png"); plt.close(fig)

    # F3 sample correlation
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    ax[0].hist(P[iu], bins=90, color="#2a7f62", alpha=.85)
    ax[0].set_xlabel("pairwise Pearson r"); ax[0].set_ylabel("sample pairs")
    ax[0].set_title("Sample-to-sample correlation")
    ax[1].scatter(np.arange(len(s3)), np.sort(mean_r.values), s=7, color="#2a7f62")
    ax[1].set_xlabel("sample (ordered)"); ax[1].set_ylabel("mean r to other samples")
    ax[1].set_title("Per-sample mean correlation")
    ordc = np.argsort(fcluster(Z, t=4, criterion="maxclust"), kind="stable")
    im = ax[2].imshow(P[np.ix_(ordc, ordc)], cmap="RdYlBu_r", vmin=np.percentile(P[iu], 1),
                      vmax=1.0, interpolation="nearest")
    ax[2].set_title("Correlation matrix (clustered order)")
    ax[2].set_xticks([]); ax[2].set_yticks([]); ax[2].grid(False)
    fig.colorbar(im, ax=ax[2], fraction=.046)
    fig.tight_layout(); fig.savefig(FIG / "F3_sample_correlation.png"); plt.close(fig)

    # F4 PCA
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.6))
    ax[0].bar(np.arange(1, 21), varexp[:20] * 100, color="#4878a8")
    ax[0].set_xlabel("PC"); ax[0].set_ylabel("% variance")
    ax[0].set_title(f"Scree (PC1-20 = {varexp[:20].sum()*100:.1f}%)")
    axb = ax[0].twinx(); axb.plot(np.arange(1, 21), np.cumsum(varexp[:20]) * 100,
                                  color="#c0392b", lw=1.2, marker="o", ms=2.5)
    axb.set_ylabel("cumulative %", color="#c0392b"); axb.grid(False)
    for a, key, ttl in ((ax[1], bt.values, "ComBat batch"),
                        (ax[2], part["dataset"].values, "deposit partition"),
                        (ax[3], cov["cimp.status"].replace({"": "unknown"}).values,
                         "CIMP status")):
        for g in sorted(set(key)):
            sel = key == g
            a.scatter(scores.PC1[sel], scores.PC2[sel], s=11, alpha=.8,
                      label=f"{g} ({sel.sum()})")
        a.set_xlabel(f"PC1 ({varexp[0]*100:.1f}%)")
        a.set_ylabel(f"PC2 ({varexp[1]*100:.1f}%)")
        a.set_title(f"PC1-PC2 by {ttl}")
        a.legend(frameon=False, fontsize=5.5, ncol=2)
    fig.tight_layout(); fig.savefig(FIG / "F4_pca.png"); plt.close(fig)

    # F5 pathway scores
    nsets = sc4.shape[1]
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.8))
    ax[0].boxplot([sc4[c].values for c in sc4.columns],
                  tick_labels=[c.replace("HALLMARK_", "")[:16] for c in sc4.columns],
                  showfliers=True, flierprops=dict(markersize=2))
    ax[0].set_ylabel("score (training-fold SD units)")
    ax[0].set_title("M2 4-set scores (EXPLORATORY full-cohort fit)")
    ax[0].tick_params(axis="x", rotation=30, labelsize=6)
    cm4 = sc4.corr()
    im = ax[1].imshow(cm4.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax[1].set_xticks(range(nsets), [c.replace("HALLMARK_", "")[:14] for c in cm4.columns],
                     rotation=45, ha="right", fontsize=6)
    ax[1].set_yticks(range(nsets), [c.replace("HALLMARK_", "")[:14] for c in cm4.columns],
                     fontsize=6)
    ax[1].grid(False); ax[1].set_title("Inter-score correlation (4-set)")
    for i in range(nsets):
        for j in range(nsets):
            ax[1].text(j, i, f"{cm4.values[i,j]:.2f}", ha="center", va="center",
                       fontsize=5.5, color="black")
    fig.colorbar(im, ax=ax[1], fraction=.046)
    emt = [c for c in sc4.columns if "EPITHELIAL" in c][0]
    e2f = [c for c in sc4.columns if "E2F" in c][0]
    sc_emt = sc4[emt].values
    for g in sorted(set(bt.values)):
        sel = bt.values == g
        ax[2].scatter(sc_emt[sel], sc4[e2f].values[sel], s=11, alpha=.8, label=g)
    ax[2].set_xlabel("EMT score"); ax[2].set_ylabel("E2F_TARGETS score")
    ax[2].set_title(f"EMT vs E2F (r = {np.corrcoef(sc_emt, sc4[e2f])[0,1]:+.2f})")
    ax[2].legend(frameon=False, fontsize=5.5, ncol=2)
    fig.tight_layout(); fig.savefig(FIG / "F5_pathway_scores.png"); plt.close(fig)

    # F6 pathway vs PC
    fig, ax = plt.subplots(1, nsets, figsize=(3.4 * nsets, 3.4))
    for a, c in zip(np.atleast_1d(ax), sc4.columns):
        sm = a.scatter(scores.PC1, scores.PC2, c=sc4[c], s=12, cmap="coolwarm")
        a.set_xlabel("PC1"); a.set_ylabel("PC2")
        a.set_title(c.replace("HALLMARK_", ""), fontsize=7)
        fig.colorbar(sm, ax=a, fraction=.046)
    fig.tight_layout(); fig.savefig(FIG / "F6_pathway_vs_pca.png"); plt.close(fig)

    print(f"  wrote {len(list(FIG.glob('*.png')))} figures to {FIG.relative_to(paths.ROOT)}")

    # =====================================================================
    section("INTEGRITY")
    # =====================================================================
    hashes_after = {p.name: sha(p) for p in FROZEN_FILES}
    changed = [k for k in hashes_before if hashes_before[k] != hashes_after[k]]
    print(f"  frozen artifacts checked : {len(FROZEN_FILES)}")
    print(f"  frozen artifacts changed : {len(changed)} {changed if changed else ''}")
    assert not changed, f"FROZEN ARTIFACT MODIFIED: {changed}"

    summary["frozen_artifact_hashes"] = hashes_after
    summary["frozen_artifacts_modified"] = False
    summary["environment"] = environment_report()
    (OUT / "phase8a_summary.json").write_text(
        json.dumps(summary, indent=2, default=str) + "\n")

    print(f"\n  outputs -> {OUT.relative_to(paths.ROOT)}/")
    for f in sorted(OUT.glob("*")):
        if f.is_file():
            print(f"    {f.name}")
    print("\nPHASE 8A COMPLETE -- no modelling, no exclusions, no outcome access.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
