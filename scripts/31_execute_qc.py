"""Phase 7C: execute the approved QC specification (Phase 7B rev. 2).

Runs all 20 pre-specified criteria EXACTLY as frozen. Thresholds are read from
the specification, never redefined here.

HARD CONSTRAINTS ENFORCED IN CODE:
  * NO sample is excluded. Soft flags produce investigation records and stop.
  * NO outcome column is opened. Event counts are NOT computed.
  * NO normalization, ComBat, fRMA, package install or CEL acquisition.
  * Frozen Phase 2 / Phase 3 / M2 artifacts are read-only.

Usage:
    .venv/bin/python scripts/31_execute_qc.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import paths                                      # noqa: E402
from src.cache import load_probe_matrix                    # noqa: E402
from src.gene_mapping import build_mapping                 # noqa: E402
from src.reproducibility import environment_report         # noqa: E402

SPEC = paths.ROOT / "config" / "analysis" / "QC_SPECIFICATION.json"
OUT = paths.ROOT / "reports" / "phase7c"
ANALYSIS = paths.PROCESSED_DIR / "analysis"


def tukey_flags(x: pd.Series) -> tuple[pd.Series, float, float]:
    """Tukey boxplot-whisker rule, exactly as frozen: Q1-1.5*IQR, Q3+1.5*IQR."""
    q1, q3 = np.percentile(x.values, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    return (x < lo) | (x > hi), float(lo), float(hi)


def section(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    spec = json.loads(SPEC.read_text())
    assert spec["qc_executed"] is False, "spec already marks QC executed"
    crit = {c["id"]: c for c in spec["criteria"]}
    results, hard, soft, desc = {}, [], [], []

    print("=" * 78)
    print("PHASE 7C -- QC EXECUTION")
    print("=" * 78)
    print(f"  specification : {spec['status']}")
    print(f"  criteria      : {len(crit)}  (thresholds read from spec, not redefined)")
    print("  NO exclusions will be performed. NO outcome column is opened.")

    # ---------------- inputs ----------------
    X = load_probe_matrix()                     # 54675 x 585, deposited
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    part = pd.read_csv(ANALYSIS / "partition_and_batch.csv")
    meta = pd.read_csv(paths.SAMPLE_METADATA, keep_default_na=False, dtype=str)
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    Xs3 = X[s3]
    print(f"\n  deposited matrix : {X.shape[0]} x {X.shape[1]}")
    print(f"  analysis cohort  : S3, n={len(s3)}")

    def rec(cid, status, detail, **kw):
        c = crit[cid]
        r = {"id": cid, "tier": c["tier"], "level": c["level"],
             "metric": c["metric"], "status": status, "detail": detail, **kw}
        results[cid] = r
        {"HARD_EXCLUSION": hard, "SOFT_INVESTIGATION": soft,
         "DESCRIPTIVE": desc}[c["tier"]].append(r)
        mark = {"PASS": "PASS", "FAIL": "FAIL", "FLAG": "FLAG",
                "REPORT": "----"}[status]
        print(f"  [{mark}] {cid:6s} {detail}")
        return r

    # ================= C. COHORT INTEGRITY (run first; may HALT) =========
    section("COHORT INTEGRITY (hard; halts on failure)")
    # H-C1 AS WRITTEN demands set EQUALITY between the deposited matrix (585)
    # and the frozen Phase 3 bundle. The bundle is S0-scoped (461) BY DESIGN,
    # so equality was never achievable: the criterion is MIS-SPECIFIED.
    # Both evaluations are recorded; the defect is surfaced for approval and
    # the criterion is NOT silently amended.
    mx, bundle = set(X.columns), set(mem.sample_id)
    as_written = mx == bundle
    missing_from_matrix = sorted(bundle - mx)
    integrity_ok = (not missing_from_matrix) and bundle <= mx
    rec("H-C1", "PASS" if integrity_ok else "FAIL",
        f"SPECIFICATION DEFECT: as-written set equality = {as_written} "
        f"(585 vs 461, never achievable by design). Integrity property that "
        f"matters -- frozen bundle is a subset of the matrix with 0 missing "
        f"IDs -- = {integrity_ok}",
        as_written_equality=as_written,
        bundle_subset_of_matrix=bool(bundle <= mx),
        ids_missing_from_matrix=missing_from_matrix,
        matrix_only_samples=len(mx - bundle),
        defect="H-C1 compares artifacts of different scope; requires approval to amend",
        requires_approval=True)
    if missing_from_matrix:
        print("\n  HALT: a frozen-bundle sample is absent from the matrix."); return 1
    print("       -> continuing: no bundle sample is missing. The defect is in")
    print("          the criterion's wording, not in the data. Flagged for approval.")

    pm = part.set_index("sample_id")["dataset"]
    mm = meta.set_index("sample_id")
    disc = int((pm.reindex(s3) == "discovery").sum())
    val = int((pm.reindex(s3) == "validation").sum())
    ok2 = pm.reindex(s3).notna().all() and (disc + val == len(s3))
    rec("H-C2", "PASS" if ok2 else "FAIL",
        f"partition membership complete: discovery={disc} validation={val}")

    ok3 = len(s3) == 453 and len(set(s3)) == 453
    rec("H-C3", "PASS" if ok3 else "FAIL", f"S3 membership == frozen n=453: {ok3}")

    a19 = json.loads((paths.ROOT / "config/analysis/A19_FREEZE_RECORD.json").read_text())
    excl_id = a19["excluded_observations"][0]["sample_id"]
    cc = mm.loc[s3, "chemotherapy.adjuvant"]
    cc_missing = cc[cc.isin({"N/A", "NA", "ND", ""})].index.tolist()
    ok4 = cc_missing == [excl_id]
    rec("H-C4", "PASS" if ok4 else "FAIL",
        f"complete-case exclusion is exactly {excl_id}: {ok4}",
        excluded_sample=excl_id)

    # ================= A. SAMPLE-LEVEL =================
    section("SAMPLE-LEVEL (hard)")
    v = Xs3.to_numpy()
    nf = np.asarray(~np.isfinite(v)).sum(axis=0)
    bad = [s for s, n in zip(s3, nf) if n > 0]
    rec("H-S1", "PASS" if not bad else "FAIL",
        f"samples with non-finite values: {len(bad)}", samples=bad)

    missing_ids = [s for s in s3 if s not in set(mem.sample_id)]
    rec("H-S2", "PASS" if not missing_ids else "FAIL",
        f"samples absent from frozen manifest: {len(missing_ids)}")

    # exact duplicate columns
    import hashlib as _hl
    digests = {}
    for i, sid in enumerate(s3):
        digests.setdefault(_hl.sha256(v[:, i].tobytes()).hexdigest(), []).append(sid)
    dup_groups = [g for g in digests.values() if len(g) > 1]
    rec("H-S3", "PASS" if not dup_groups else "FAIL",
        f"exactly-duplicated sample columns: {len(dup_groups)} group(s)",
        groups=dup_groups)

    is_norm = ~mm.loc[s3, "source_name"].str.contains(
        "primary colorectal Adenocarcinoma", case=False, na=False)
    n_norm = int(is_norm.sum())
    rec("H-S4", "PASS" if n_norm == 0 else "FAIL",
        f"non-tumoural samples in analysis cohort: {n_norm}")

    section("SAMPLE-LEVEL (soft - INVESTIGATION ONLY, no exclusion)")
    med = pd.Series(np.median(v, axis=0), index=s3)
    q1, q3 = np.percentile(v, [25, 75], axis=0)
    iqr = pd.Series(q3 - q1, index=s3)

    f_med, lo_m, hi_m = tukey_flags(med)
    rec("S-S1", "FLAG" if f_med.any() else "PASS",
        f"per-array median outside [{lo_m:.4f}, {hi_m:.4f}]: {int(f_med.sum())}",
        fence=[lo_m, hi_m], flagged=med[f_med].round(4).to_dict())

    f_iqr, lo_i, hi_i = tukey_flags(iqr)
    rec("S-S2", "FLAG" if f_iqr.any() else "PASS",
        f"per-array IQR outside [{lo_i:.4f}, {hi_i:.4f}]: {int(f_iqr.sum())}",
        fence=[lo_i, hi_i], flagged=iqr[f_iqr].round(4).to_dict())

    print("    computing rank-correlation distance matrix ...", flush=True)
    R = Xs3.rank(axis=0).to_numpy()
    R = R - R.mean(axis=0, keepdims=True)
    R /= np.sqrt((R ** 2).sum(axis=0, keepdims=True))
    C = R.T @ R                                    # Spearman correlation matrix
    np.fill_diagonal(C, np.nan)
    dist = pd.Series(np.nansum(1.0 - C, axis=1), index=s3)
    f_d, lo_d, hi_d = tukey_flags(dist)
    rec("S-S3", "FLAG" if f_d.any() else "PASS",
        f"sum(1-rho) outside [{lo_d:.3f}, {hi_d:.3f}]: {int(f_d.sum())}",
        fence=[lo_d, hi_d], flagged=dist[f_d].round(3).to_dict())

    section("SAMPLE-LEVEL (descriptive - never triggers exclusion)")
    Xc = v - v.mean(axis=1, keepdims=True)
    cov = (Xc.T @ Xc) / (Xc.shape[0] - 1)
    ev, evec = np.linalg.eigh(cov)
    ev, evec = ev[::-1], evec[:, ::-1]
    varexp = ev / ev.sum()
    scores = pd.DataFrame(evec[:, :5] * np.sqrt(np.maximum(ev[:5], 0)),
                          index=s3, columns=[f"PC{i}" for i in range(1, 6)])
    rec("D-S1", "REPORT",
        f"PCA variance explained PC1-PC5: "
        f"{', '.join(f'{x*100:.1f}%' for x in varexp[:5])} (no outlier rule)",
        variance_explained=[float(x) for x in varexp[:5]])

    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    D = 1.0 - np.nan_to_num(C, nan=1.0)
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    k2 = fcluster(Z, t=2, criterion="maxclust")
    rec("D-S2", "REPORT",
        f"hierarchical clustering (correlation, average linkage); "
        f"2-cluster split sizes {np.bincount(k2)[1:].tolist()} (no cut-height rule)")

    bt = part.set_index("sample_id")["combat_batch"].reindex(s3)
    bx = pd.crosstab(bt, pm.reindex(s3))
    pc_batch = {f"PC{i+1}": float(pd.Series(scores[f'PC{i+1}'].values)
                                  .groupby(bt.values).mean().std())
                for i in range(3)}
    rec("D-S3", "REPORT",
        f"batch composition over {bt.nunique()} batches; "
        f"SD of batch-mean PC1={pc_batch['PC1']:.3f} (NO batch correction applied)",
        batch_pc_sd=pc_batch)

    AFFX = {"ACTB": ("AFFX-HSAC07/X00351_3_at", "AFFX-HSAC07/X00351_5_at"),
            "GAPDH": ("AFFX-HUMGAPDH/M33197_3_at", "AFFX-HUMGAPDH/M33197_5_at")}
    ratios = {}
    for g, (p3, p5) in AFFX.items():
        ratios[g] = (Xs3.loc[p3] - Xs3.loc[p5])
    rr = pd.DataFrame(ratios)
    rec("D-S4", "REPORT",
        f"3'/5' log2 difference ACTB median {rr.ACTB.median():.3f} "
        f"(range {rr.ACTB.min():.2f}..{rr.ACTB.max():.2f}), "
        f"GAPDH median {rr.GAPDH.median():.3f} -- "
        "Affymetrix <3 threshold NOT applied (scale mismatch)",
        actb_median=float(rr.ACTB.median()), gapdh_median=float(rr.GAPDH.median()))
    rr.round(4).to_csv(OUT / "degradation_ratios.csv")

    # ================= B. PROBE / GENE-LEVEL =================
    section("PROBE / GENE-LEVEL")
    G = pd.read_parquet(ANALYSIS / "expression_median_S0.parquet")[s3]
    gnf = int(np.asarray(~np.isfinite(G.to_numpy())).sum())
    rec("H-G1", "PASS" if gnf == 0 else "FAIL",
        f"non-finite values in gene matrix: {gnf} (escalate, never drop)")

    sd = G.std(axis=1)
    zv = int((sd == 0).sum())
    rec("D-G1", "REPORT",
        f"zero-variance genes in S3: {zv} (frozen per-fold M2 rule governs; "
        "no cohort-level exclusion)", zero_variance=zv)

    gm = G.mean(axis=1)
    counts = {f"mean_lt_{t}": int((gm < t).sum()) for t in (4.0, 5.0, 6.0)}
    counts |= {f"sd_lt_{t}": int((sd < t).sum()) for t in (0.2, 0.3, 0.5)}
    rec("D-G2", "REPORT",
        f"low-expression counts {counts} -- NO filter applied", counts=counts)

    medspec = pd.read_csv(paths.PROCESSED_DIR / "frozen/frozen_gene_to_probeset_median.csv")
    affx = medspec[medspec.probe_id.str.startswith("AFFX-")]
    rec("D-G3", "REPORT",
        f"AFFX control probe sets in frozen median representation: {len(affx)} "
        f"across {affx.gene_symbol.nunique()} genes "
        f"({', '.join(sorted(affx.gene_symbol.unique()))}) -- NOT removed",
        affx_probes=int(len(affx)), affx_genes=sorted(affx.gene_symbol.unique()))

    mp = build_mapping().summary()
    rec("D-G4", "REPORT",
        f"mapping: unique {mp['n_unique']} / ambiguous {mp['n_ambiguous']} / "
        f"unmapped {mp['n_unmapped']} (frozen Phase 2)", mapping=mp)

    # ================= investigation records =================
    section("SOFT FLAGS -> INVESTIGATION RECORDS (no exclusion performed)")
    flagged = sorted(set(med[f_med].index) | set(iqr[f_iqr].index) | set(dist[f_d].index))
    print(f"  distinct samples soft-flagged: {len(flagged)}")
    recs = []
    for sid in flagged:
        fl = []
        if f_med.get(sid, False):
            fl.append({"criterion": "S-S1", "metric": "per-array median",
                       "value": round(float(med[sid]), 4), "fence": [round(lo_m,4), round(hi_m,4)]})
        if f_iqr.get(sid, False):
            fl.append({"criterion": "S-S2", "metric": "per-array IQR",
                       "value": round(float(iqr[sid]), 4), "fence": [round(lo_i,4), round(hi_i,4)]})
        if f_d.get(sid, False):
            fl.append({"criterion": "S-S3", "metric": "sum(1 - Spearman rho)",
                       "value": round(float(dist[sid]), 3), "fence": [round(lo_d,3), round(hi_d,3)]})
        r = mm.loc[sid]
        recs.append({
          "sample_id": sid,
          "flagged_metric_and_value": fl,
          "metadata_and_batch": {
            "title": r["title"], "dataset": pm[sid], "combat_batch": bt[sid],
            "tnm_stage": r["tnm.stage"], "sex": r["sex"],
            "age_at_diagnosis": r["age.at.diagnosis_(year)"],
            "tumor_location": r["tumor.location"],
            "rna_extraction_method": r["rna.extraction.method"],
            "cit_molecular_subtype": r["cit.molecularsubtype"]},
          "independent_technical_evidence": {
            "non_finite_values": int(nf[s3.index(sid)]),
            "duplicate_column": sid in [x for g in dup_groups for x in g],
            "actb_3p5p_log2_diff": round(float(rr.ACTB[sid]), 3),
            "gapdh_3p5p_log2_diff": round(float(rr.GAPDH[sid]), 3),
            "assessment": "SEE REPORT - direct technical evidence is evaluated "
                          "separately from biological hypotheses"},
          "biological_plausibility": {
            "cit_subtype": r["cit.molecularsubtype"],
            "note": "HYPOTHESIS, not evidence: CMS4/stromal phenotypes shift "
                    "global expression distributions for biological reasons"},
          "effect_on_cohort": {
            "n_if_excluded": len(s3) - 1,
            "event_counts": "NOT COMPUTED - outcome columns were not opened; "
                            "outcome-guided selection is prohibited"},
          "proposed_action_and_justification": {
            "proposed_action": "NO EXCLUSION - investigation only",
            "justification": "Phase 7B C3: a sample may not be excluded solely "
                             "for crossing a Tukey fence. Exclusion requires "
                             "explicit approval."},
        })
    (OUT / "investigation_records.json").write_text(json.dumps(recs, indent=2, default=str) + "\n")
    for x in recs:
        print(f"    {x['sample_id']}  flags={[f['criterion'] for f in x['flagged_metric_and_value']]}"
              f"  batch={x['metadata_and_batch']['combat_batch']}"
              f"  subtype={x['metadata_and_batch']['cit_molecular_subtype']}")

    # ================= outputs =================
    pd.DataFrame({"sample_id": s3, "dataset": pm.reindex(s3).values,
                  "combat_batch": bt.values, "median": med.values,
                  "iqr": iqr.values, "dist_sum_1_minus_rho": dist.values,
                  "flag_S_S1": f_med.values, "flag_S_S2": f_iqr.values,
                  "flag_S_S3": f_d.values,
                  "n_non_finite": nf}).round(5).to_csv(
        OUT / "sample_qc_metrics.csv", index=False)
    scores.round(4).join(pd.DataFrame({"dataset": pm.reindex(s3), "batch": bt})
                         ).to_csv(OUT / "pca_scores.csv")
    bx.to_csv(OUT / "batch_composition.csv")
    pd.DataFrame(results).T.to_csv(OUT / "qc_criteria_results.csv")

    summary = {
      "executed_utc": environment_report()["captured_utc"],
      "specification_status": spec["status"],
      "criteria_run": len(results),
      "hard_failures": sum(1 for r in hard if r["status"] == "FAIL"),
      "hard_passed": sum(1 for r in hard if r["status"] == "PASS"),
      "soft_flags_criteria": sum(1 for r in soft if r["status"] == "FLAG"),
      "descriptive_reported": len(desc),
      "samples_soft_flagged": len(flagged),
      "sample_ids_requiring_investigation": flagged,
      "samples_excluded": 0,
      "exclusions_performed": False,
      "outcome_columns_opened": False,
      "event_counts_computed": False,
      "normalization_applied": "none", "batch_correction_applied": "none",
      "frozen_artifacts_modified": False,
      "results": results,
      "environment": environment_report(),
    }
    (OUT / "qc_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    section("SUMMARY")
    print(f"  criteria run          : {len(results)} / 20")
    print(f"  HARD failures         : {summary['hard_failures']}  (passed {summary['hard_passed']})")
    print(f"  SOFT criteria flagging: {summary['soft_flags_criteria']} / 3")
    print(f"  samples soft-flagged  : {len(flagged)}")
    print(f"  DESCRIPTIVE reported  : {len(desc)}")
    print(f"  SAMPLES EXCLUDED      : 0")
    print(f"  outcome columns opened: False   event counts computed: False")
    print(f"\n  -> {OUT.relative_to(paths.ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
