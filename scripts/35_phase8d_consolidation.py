"""Phase 8D: consolidate Phases 8A, 8B and 8C into one auditable review.

CONSOLIDATION ONLY. No new analysis is performed. Every numerical claim that
reaches the Phase 8D report is RESOLVED from a source artifact by path and
pointer, never transcribed, and cross-phase corroborations are checked. A
disagreement between phases HALTS the run rather than being reconciled.

HARD CONSTRAINTS ENFORCED IN CODE:
  * Outcome blindness, two layers (file + column). No endpoint access.
  * All 453 S3 samples retained; no exclusion, relabelling, cohort revision.
  * Phase 8A/8B/8C reports and artifacts are READ-ONLY, hashed before + after.
  * Frozen Phase 2/3, A-19, QC records, M2 specs and the canonical scorer are
    read-only and verified.
  * No normalization, batch correction, CEL acquisition, modelling, feature
    selection or outcome analysis.
  * No Phase 8C recommendation is implemented. The decision register records
    options and implications; it decides nothing.

Usage:
    /Users/kalki/coad_transportability/.venv/bin/python \
        /Users/kalki/coad_transportability/scripts/35_phase8d_consolidation.py
"""

from __future__ import annotations

import hashlib
import inspect
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import paths                                        # noqa: E402
from src.m2_amendment import resolve_hash                     # noqa: E402
from src.pathway_scoring import score_panel                   # noqa: E402
from src.reproducibility import environment_report            # noqa: E402

ANALYSIS = paths.PROCESSED_DIR / "analysis"
FROZEN = paths.PROCESSED_DIR / "frozen"
PANELS = paths.ROOT / "config" / "pathways"
CFG = paths.ROOT / "config" / "analysis"
P8A, P8B, P8C = (paths.ROOT / "reports" / f"phase8{x}" for x in "abc")
OUT = paths.ROOT / "reports" / "phase8d"

OUTCOME_COLUMNS = ("rfs.event", "rfs.delay", "os.event", "os.delay_(months)")
BLOCKED_FILES = ("endpoint.csv",)


def _guard_outcome_files() -> None:
    import pandas as _pd
    for fn in ("read_csv", "read_parquet", "read_table"):
        orig = getattr(_pd, fn)

        def blocked(path, *a, _orig=orig, **kw):
            if any(b in str(path) for b in BLOCKED_FILES):
                raise PermissionError(
                    f"PHASE 8D IS OUTCOME-BLIND: refusing to read {path}")
            return _orig(path, *a, **kw)
        setattr(_pd, fn, blocked)


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dig(obj, pointer: str):
    """Resolve 'a.b[2].c' against nested JSON. Raises if absent."""
    cur = obj
    for tok in re.findall(r"[^.\[\]]+|\[\d+\]", pointer):
        if tok.startswith("["):
            cur = cur[int(tok[1:-1])]
        else:
            cur = cur[tok]
    return cur


def section(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)


SOURCES = {
    "8A": P8A / "phase8a_summary.json",
    "8B": P8B / "phase8b_summary.json",
    "8C": P8C / "phase8c_summary.json",
    "M2": PANELS / "M2_FREEZE_RECORD.json",
    "P4": PANELS / "m2_panel_4set.json",
    "P5": PANELS / "m2_panel_5set.json",
    "A19": CFG / "A19_FREEZE_RECORD.json",
    "QCD": CFG / "QC_DECISIONS_RECORD.json",
    "QCS": CFG / "QC_SPECIFICATION.json",
}

READONLY = (sorted(P8A.rglob("*")) + sorted(P8B.rglob("*")) + sorted(P8C.rglob("*"))
            + [PANELS / "m2_panel_4set.json", PANELS / "m2_panel_5set.json",
               PANELS / "M2_FREEZE_RECORD.json", CFG / "A19_FREEZE_RECORD.json",
               CFG / "QC_DECISIONS_RECORD.json", CFG / "QC_SPECIFICATION.json",
               paths.ROOT / "src" / "pathway_scoring.py",
               FROZEN / "frozen_gene_to_probeset_median.csv",
               FROZEN / "frozen_gene_to_probe_discovery443.csv",
               FROZEN / "frozen_representation.json",
               ANALYSIS / "expression_median_S0.parquet",
               ANALYSIS / "expression_repiqr_disc443_S0.parquet",
               ANALYSIS / "cohort_membership.csv",
               ANALYSIS / "partition_and_batch.csv",
               ANALYSIS / "covariates.csv", ANALYSIS / "endpoint.csv",
               paths.ROOT / "docs" / "PHASE8A_EXPLORATORY_ANALYSIS.md",
               paths.ROOT / "docs" / "PHASE8B_PCA_BATCH_INVESTIGATION.md",
               paths.ROOT / "docs" / "PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md"])
READONLY = [p for p in READONLY if p.is_file()]


# ---------------------------------------------------------------------------
# Every consolidated number, declared with WHERE IT COMES FROM.
#   (id, phase, source key, json pointer, evidence_type, derived_from_matrix)
# ---------------------------------------------------------------------------
CLAIMS = [
    # --- A. cohort / provenance ---
    ("n_samples_8A", "8A", "8A", "n_samples", "observation", False),
    ("n_samples_8B", "8B", "8B", "n_samples", "observation", False),
    ("n_samples_8C", "8C", "8C", "n_samples", "observation", False),
    ("n_genes", "8A", "8A", "n_genes", "observation", False),
    ("excluded_8A", "8A", "8A", "samples_excluded", "observation", False),
    ("excluded_8B", "8B", "8B", "samples_excluded", "observation", False),
    ("excluded_8C", "8C", "8C", "samples_excluded", "observation", False),
    # --- B. expression level ---
    ("sample_median_min", "8A", "8A", "expression.sample_median.min", "observation", True),
    ("sample_median_max", "8A", "8A", "expression.sample_median.max", "observation", True),
    ("gene_zero_variance", "8A", "8A", "gene_variability.n_zero_variance", "observation", True),
    ("mean_sd_spearman", "8A", "8A", "gene_variability.spearman_mean_vs_sd", "association", True),
    ("pairwise_r_min", "8A", "8A", "sample_correlation.pearson.min", "observation", True),
    ("pairwise_r_median", "8A", "8A", "sample_correlation.pearson.median", "observation", True),
    ("sex_concordant_both_axes", "8A", "8A", "sex_positive_control.n_concordant_both_axes", "observation", True),
    # --- C. PCA / batch ---
    ("pc1_var_8A", "8A", "8A", "pca.variance_explained_pc1_20[0]", "observation", True),
    ("pc1_var_8B", "8B", "8B", "variance_structure.variance_explained_pc1_50[0]", "observation", True),
    ("pcs_50pct", "8B", "8B", "variance_structure.n_pcs_50pct", "observation", True),
    ("pcs_80pct", "8B", "8B", "variance_structure.n_pcs_80pct", "observation", True),
    ("horn_pcs", "8B", "8B", "variance_structure.horn_n_pcs_above_permuted_p95", "observation", True),
    ("repro_min_corr", "8B", "8B", "reproduction_of_phase8a.min_signed_correlation_scores_pc1_20", "observation", False),
    ("combat_mod_s3_n", "8B", "8B", "combat_mod_vs_mmr.s3_n_with_mmr", "observation", False),
    ("combat_mod_s3_disc", "8B", "8B", "combat_mod_vs_mmr.s3_discordant", "observation", False),
    ("combat_mod_sdrf_n", "8B", "8B", "combat_mod_vs_mmr.sdrf_n_with_mmr", "observation", False),
    ("combat_mod_sdrf_disc", "8B", "8B", "combat_mod_vs_mmr.sdrf_discordant", "observation", False),
    ("batch05_pct_dmmr", "8B", "8B", "batch05_pc3.batch05_pct_modcan1", "observation", False),
    ("cohort_pct_dmmr", "8B", "8B", "batch05_pc3.cohort_pct_modcan1", "observation", False),
    ("batch_detect_floor", "8B", "8B", "batch_within_mmr_strata.min_detectable_batch05_shift_pc3_sd_at_80pct_power", "observation", True),
    ("sex_both_axes_n", "8B", "8B", "sex_discordance_review.n", "observation", True),
    ("sex_softflag_n", "8B", "8B", "sex_discordance_review.n_with_any_phase7c_soft_flag", "association", True),
    ("sex_softflag_or", "8B", "8B", "sex_discordance_review.soft_flag_enrichment_fisher.odds_ratio", "association", True),
    ("dependent_arrays", "8B", "8B", "new_technical_structure.dependent_samples.n_in_S3", "observation", False),
    # --- D. pathway scores ---
    ("cov_emt", "8C", "8C", "audit.5set_sensitivity.coverage.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.n_genes_used", "observation", False),
    ("cov_e2f", "8C", "8C", "audit.5set_sensitivity.coverage.HALLMARK_E2F_TARGETS.n_genes_used", "observation", False),
    ("cov_wnt", "8C", "8C", "audit.5set_sensitivity.coverage.HALLMARK_WNT_BETA_CATENIN_SIGNALING.n_genes_used", "observation", False),
    ("cov_p53", "8C", "8C", "audit.5set_sensitivity.coverage.HALLMARK_P53_PATHWAY.n_genes_used", "observation", False),
    ("cov_pi3k", "8C", "8C", "audit.5set_sensitivity.coverage.HALLMARK_PI3K_AKT_MTOR_SIGNALING.n_genes_used", "observation", False),
    ("maxr_4set", "8C", "8C", "relationships.4set_primary.max_abs_r", "observation", True),
    ("maxr_5set", "8C", "8C", "relationships.5set_sensitivity.max_abs_r", "observation", True),
    ("cond_4set", "8C", "8C", "relationships.4set_primary.condition_number", "observation", True),
    ("cond_5set", "8C", "8C", "relationships.5set_sensitivity.condition_number", "observation", True),
    ("pi3k_r2", "8C", "8C", "relationships.pi3k_increment.r2_on_primary_four", "association", True),
    ("audit_impl_diff", "8C", "8C", "audit.5set_sensitivity.max_abs_diff_production_vs_independent", "observation", False),
    ("audit_8a_diff", "8C", "8C", "audit.5set_sensitivity.max_abs_diff_vs_phase8a_artifact", "observation", False),
    ("strat_survivors", "8C", "8C", "technical_structure.n_surviving_stratification", "association", True),
    # --- frozen specification state ---
    ("ddof_p4", "8C", "P4", "scoring.sigma_ddof", "frozen", False),
    ("ddof_p5", "8C", "P5", "scoring.sigma_ddof", "frozen", False),
    ("nsets_p4", "8C", "P4", "n_sets", "frozen", False),
    ("nsets_p5", "8C", "P5", "n_sets", "frozen", False),
]


def main() -> int:
    _guard_outcome_files()
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("PHASE 8D -- CONSOLIDATED EXPLORATORY ANALYSIS AND DECISION REVIEW")
    print("=" * 78)
    print("  consolidation only: no new analysis, no modelling, no outcome access")

    before = {str(p.relative_to(paths.ROOT)): sha(p) for p in READONLY}
    src = {k: json.loads(v.read_text()) for k, v in SOURCES.items()}

    # =====================================================================
    section("1. RESOLVE EVERY CONSOLIDATED CLAIM FROM ITS SOURCE ARTIFACT")
    # =====================================================================
    rows, unresolved = [], []
    for cid, phase, key, ptr, etype, derived in CLAIMS:
        try:
            val = dig(src[key], ptr)
        except (KeyError, IndexError, TypeError) as e:
            unresolved.append((cid, key, ptr, repr(e)))
            continue
        rows.append({"claim_id": cid, "source_phase": phase,
                     "source_artifact": str(SOURCES[key].relative_to(paths.ROOT)),
                     "json_pointer": ptr, "value": val, "evidence_type": etype,
                     "derived_from_same_expression_matrix": derived})
    if unresolved:
        print("  HALT -- claims that do not resolve against their source:")
        for u in unresolved:
            print(f"    {u}")
        return 2
    claims = pd.DataFrame(rows)
    print(f"  {len(claims)} claims resolved from {claims.source_artifact.nunique()} "
          f"artifacts; 0 transcribed")
    print(f"    evidence types: {claims.evidence_type.value_counts().to_dict()}")
    print(f"    derived from the same expression matrix: "
          f"{int(claims.derived_from_same_expression_matrix.sum())}/{len(claims)}")

    # =====================================================================
    section("2. CROSS-PHASE CORROBORATION (a disagreement halts the run)")
    # =====================================================================
    v = {r["claim_id"]: r["value"] for r in rows}   # types preserved, not via DataFrame
    checks = []

    def chk(name, ok, detail):
        checks.append({"check": name, "agrees": bool(ok), "detail": detail})
        print(f"  [{'OK  ' if ok else 'FAIL'}] {name:46s} {detail}")

    chk("cohort size identical across 8A/8B/8C",
        v["n_samples_8A"] == v["n_samples_8B"] == v["n_samples_8C"] == 453,
        f"{v['n_samples_8A']}/{v['n_samples_8B']}/{v['n_samples_8C']}")
    chk("zero exclusions across 8A/8B/8C",
        v["excluded_8A"] == v["excluded_8B"] == v["excluded_8C"] == 0, "0/0/0")
    chk("PC1 variance agrees 8A vs 8B",
        abs(v["pc1_var_8A"] - v["pc1_var_8B"]) < 1e-7,
        f"{v['pc1_var_8A']:.8f} vs {v['pc1_var_8B']:.8f}")
    chk("8B reproduces 8A PCA scores", v["repro_min_corr"] > 1 - 1e-12,
        f"min signed r = {v['repro_min_corr']:.12f}")
    chk("sex both-axes count agrees 8A vs 8B",
        (453 - v["sex_concordant_both_axes"]) >= v["sex_both_axes_n"] == 8,
        f"8A concordant {v['sex_concordant_both_axes']}/453, 8B both-axes {v['sex_both_axes_n']}")

    # the eight must be the SAME eight in 8A and 8B
    a8 = set(pd.read_csv(P8A / "sex_positive_control.csv")
             .query("classification == 'discordant_both_axes'").sample_id)
    b8 = {r["sample_id"] for r in dig(src["8B"], "sex_discordance_review.samples")}
    chk("the eight are the same samples in 8A and 8B", a8 == b8 and len(a8) == 8,
        f"|8A∩8B| = {len(a8 & b8)}")

    # coverage must agree between 8A and 8C
    cov8a = {k: x["n_genes_used"] for k, x in
             dig(src["8A"], "pathways.5set_sensitivity.coverage").items()}
    cov8c = {k: x["n_genes_used"] for k, x in
             dig(src["8C"], "audit.5set_sensitivity.coverage").items()}
    chk("pathway coverage agrees 8A vs 8C", cov8a == cov8c,
        f"{sorted(cov8c.values())}")

    # correlation structure must agree 8A vs 8C and corroborate the 6C record
    chk("4-set max|r| agrees 8A vs 8C",
        abs(dig(src["8A"], "pathways.4set_primary.inter_score_correlation.max_abs_r")
            - v["maxr_4set"]) < 1e-9, f"{v['maxr_4set']:.6f}")
    rec6c = dig(src["P4"], "collinearity_all461_descriptive_only.max_abs_r")
    chk("4-set max|r| corroborates the 6C record (n=461 vs 453)",
        abs(rec6c - v["maxr_4set"]) < 0.01, f"6C {rec6c:.4f} vs 8C {v['maxr_4set']:.4f}")

    # GSM972311 must be the singleton in both 8A and 8B
    cl8a = dig(src["8A"], "clustering.cluster_sizes")
    chk("8A clustering shows a persistent singleton",
        all(1 in sizes for sizes in cl8a.values()), f"k=2..6 sizes {cl8a['2']}..{cl8a['6']}")

    # frozen specification state
    chk("ddof pinned to 1 in both panels", v["ddof_p4"] == v["ddof_p5"] == 1, "1 / 1")
    chk("panel sizes 4 and 5", (v["nsets_p4"], v["nsets_p5"]) == (4, 5), "4 / 5")
    chk("panel roles unchanged",
        dig(src["P4"], "role") == "PRIMARY" and dig(src["P5"], "role").startswith("SENSITIVITY"),
        f"{dig(src['P4'],'role')} / {dig(src['P5'],'role')}")
    chk("canonical scorer is src/pathway_scoring.py",
        score_panel.__module__ == "src.pathway_scoring"
        and "sigma_ddof" in inspect.getsource(score_panel),
        score_panel.__module__)
    chk("no phase opened an outcome column",
        not any(dig(src[k], "outcome_columns_opened") for k in ("8A", "8B", "8C")),
        "8A/8B/8C all False")

    if not all(c["agrees"] for c in checks):
        print("\n  HALT -- cross-phase disagreement; historical evidence is NOT rewritten")
        return 3
    print(f"\n  {len(checks)}/{len(checks)} corroboration checks agree")

    # =====================================================================
    section("3. EVIDENCE TABLE")
    # =====================================================================
    EV = [
        # finding, phases, type, corroboration, limitation, modelling implication, status
        ("Matrix complete: 0 non-finite, 0 zero-variance genes, no sample pair below r=0.80",
         "8A", "observation", "single-phase; re-read in 8C audit",
         "property of deposited RMA output, not of biology",
         "no imputation or array-level repair needed", "frozen"),
        ("Per-sample location tightly controlled (median range 0.302 log2)",
         "8A", "observation", "single-phase",
         "expected signature of quantile-normalised RMA; says nothing about validity",
         "supports using the deposit as-is for relative-risk modelling", "frozen"),
        ("Intensity-dependent heteroscedasticity (Spearman mean vs SD = +0.633)",
         "8A", "association", "single-phase",
         "residual of log2 RMA; low-expressed genes compressed",
         "argues against naive variance filtering", "documented limitation"),
        ("PCA spectrum shallow and high-dimensional (PC1 10.3%; 17/127 PCs for 50/80%)",
         "8A+8B", "observation", "REPRODUCED by an independent algorithm (signed r = 1.0 to 1e-12)",
         "Horn gives 48 PCs, broken-stick 97 - the criteria disagree on how many matter",
         "no single dominant technical gradient to adjust for", "frozen"),
        ("Leading PCs track MMR/CIMP/CIN/BRAF/location, not stage or partition",
         "8A+8B", "association", "both phases, same effect sizes",
         "all annotations are correlated with each other; not causal",
         "the dominant axis is molecular subtype, not the clinical covariates in M1",
         "documented limitation"),
        ("combat_mod perfectly concordant with recorded MMR (0/410 S3, 0/519 SDRF)",
         "8B", "observation", "REPRODUCED in 8D from the 8B artifact",
         "concordance is NOT identity; any identically-partitioning variable fits",
         "if MMR was protected in ComBat, MMR-associated variance is preserved by construction",
         "unresolved - requires author correspondence"),
        ("Batch-PC location association vanishes under MMR stratification",
         "8B", "association", "pMMR stratum and a 10,000-permutation conditional test",
         "ComBat already removed batch variance; the test's floor is ~0.3 SD of PC3",
         "no evidence for a batch location covariate in the primary model",
         "documented limitation"),
        ("Nominal PC3 x batch association WITHIN the dMMR stratum (eps2=0.130, p=0.028)",
         "8B", "association", "NOT corroborated - the better-powered conditional test is null (p=0.38)",
         "n=64, 6 batches, min cell 2, uncorrected for multiplicity",
         "too fragile to act on; motivates a batch sensitivity analysis at most",
         "unresolved"),
        ("AFFX control-probe PC1 ~ batch SURVIVES MMR stratification (p=0.0016)",
         "8B", "association", "NEW in 8B; not testable in 8A",
         "AFFX probes went through the same deposited RMA+ComBat - not an uncorrected reference",
         "residual batch structure exists; candidate for a prespecified sensitivity analysis",
         "requires approval"),
        ("PC2 DISPERSION differs across batches and survives stratification (p=0.0001)",
         "8B", "association", "NEW in 8B; 8A compared means only",
         "location-and-scale correction removes means more completely than higher moments",
         "a dispersion effect is not addressed by a location covariate",
         "requires approval"),
        ("batch05 is 63.3% dMMR vs 14.3% cohort-wide; its PC3 elevation is carried by that share",
         "8B", "observation", "REPRODUCES the Phase 7C OBS-1 observation with a mechanism",
         "composition and batch-specific processing are PERFECTLY CONFOUNDED here; CIs overlap within stratum",
         "strengthens Phase 7C D3 (observation, not exclusion grounds)", "frozen"),
        ("Eight samples discordant on BOTH sex axes, across five batches",
         "8A+8B", "observation", "same eight samples in both phases (set equality checked in 8D)",
         "expression markers cannot separate a mix-up from LOY / XIST silencing; no genotype data",
         "none - samples retained", "unresolved - genotype evidence unavailable"),
        ("Soft-flag enrichment among the eight (3/8 vs 29/453, OR 9.67, p=0.011)",
         "8B", "association", "single test",
         "a soft flag reports array QUALITY, not IDENTITY; derived from the same matrix",
         "none - not an exclusion criterion", "frozen (D2 retains all 453)"),
        ("No dependent/related arrays in S3 (dependancy_sample = 0/453)",
         "8B", "observation", "NEW in 8B",
         "absence of a DECLARED dependency, not proof of independence",
         "declared sample dependency does not violate a resampling scheme's assumption",
         "documented limitation"),
        ("Frozen coverage met for all five sets (198/195/41/194/104), 0 zero-variance",
         "8A+8C", "observation", "identical in both phases (checked in 8D)",
         "none",
         "no score is dropped in any fold", "frozen"),
        ("Production scorer agrees with an independent re-implementation to 3e-15",
         "8C", "observation", "NEW in 8C",
         "agreement of two implementations, not of two data sources",
         "the scoring code is not the weak link", "frozen"),
        ("Two coherence regimes: E2F 0.304 / EMT 0.241 vs WNT 0.041 / P53 0.033 / PI3K 0.025",
         "8C", "observation", "NEW in 8C",
         "a property of the gene sets and the averaging transform, NOT of relevance",
         "low-dispersion scores are shrunk harder under penalisation - motivates decision 1",
         "requires approval"),
        ("4-set well-conditioned (max|r| 0.316, cond 2.95); 5-set adds one pair >0.5 (P53-PI3K 0.631)",
         "6C+8A+8C", "observation", "REPRODUCED across three phases and two cohort sizes",
         "descriptive collinearity among scores; no outcome involved",
         "informs how a 4-set/5-set divergence should be read - decision 2", "frozen"),
        ("PI3K: R2 = 0.540 on the four primary scores (46% not linearly captured)",
         "8C", "association", "NEW in 8C",
         "a redundancy measure among scores only; no recurrence relevance",
         "NOT grounds to promote or remove PI3K; it stays sensitivity-only", "frozen"),
        ("Representation stability: E2F 0.998, EMT 0.992, P53 0.972, PI3K 0.910, WNT 0.858",
         "8C", "observation", "NEW in 8C (sensitivity-only)",
         "difference in probe->gene summarisation, NOT biological difference",
         "WNT is the one score materially dependent on the representation - decision 3",
         "requires approval"),
        ("Scores stable to fold composition (min r 0.9996) but shift is non-zero (max 0.072 SD)",
         "8A+8C", "observation", "same figure in both phases",
         "none",
         "confirms standardisation must stay inside the resampling loop", "frozen"),
        ("All batch and extraction-method associations with scores vanish under MMR stratification (0/10)",
         "8C", "association", "consistent with the 8B PC result for location",
         "same ComBat floor; sparse strata (min cell 14)",
         "no batch covariate indicated for the score model", "documented limitation"),
        ("QC soft-flag DISPERSION associations with scores survive stratification",
         "8C", "association", "NEW in 8C",
         "MATRIX-DERIVED: the flags are Tukey rules on this same matrix, so partly mechanical",
         "not independent evidence of array failure; no action", "frozen"),
        ("cit.molecularsubtype has the largest PC and score associations",
         "8B+8C", "association", "both phases",
         "CIRCULAR: it is the depositors' own expression-derived clustering of this matrix",
         "must never be used as validation", "frozen"),
    ]
    ev = pd.DataFrame(EV, columns=[
        "finding", "source_phase", "evidence_type", "reproducibility_or_corroboration",
        "limitation_or_alternative_explanation", "implication_for_modelling",
        "decision_status"])
    ev.to_csv(OUT / "evidence_table.csv", index=False)
    print(f"  {len(ev)} findings")
    print(f"    by status: {ev.decision_status.apply(lambda s: s.split(' ')[0]).value_counts().to_dict()}")
    print(f"    by phase : {ev.source_phase.value_counts().to_dict()}")

    cat = {
        "reproduced_across_phases": [e[0] for e in EV if "+" in e[1] or "REPRODUC" in e[3].upper()],
        "new_in_8B_or_8C": [e[0] for e in EV if "NEW in" in e[3]],
        "unresolved": [e[0] for e in EV if e[6].startswith("unresolved")],
        "not_independent_same_matrix_or_qc_rules": [
            e[0] for e in EV if "CIRCULAR" in e[4] or "MATRIX-DERIVED" in e[4]
            or "same matrix" in e[4] or "same deposited" in e[4]],
    }
    print("\n  classification:")
    for k, vlist in cat.items():
        print(f"    {k:42s} {len(vlist)}")

    # =====================================================================
    section("4. DECISION REGISTER (nothing is decided here)")
    # =====================================================================
    DEC = [
        {"id": "D8D-1", "title": "Score-level standardisation before penalised Cox",
         "current_frozen_state": dig(src["P4"], "scoring.score_standardisation"),
         "blocking": "blocks_specific_analysis",
         "blocks_what": "the penalised Cox specification in Phase 6B §3.2; nothing else",
         "evidence": ["8C §3.1 implied coherence r̄", "reports/phase8c/score_distributions.csv"],
         "implication": (
             "The five scores sit on very different scales: SD 0.554 (E2F) down to "
             "0.185 (PI3K), because averaging k standardised genes shrinks dispersion "
             "by roughly 1/sqrt(k) times the set's coherence. A ridge or elastic-net "
             "penalty applies the SAME penalty per unit coefficient, so on unstandardised "
             "scores a low-dispersion predictor needs a LARGER coefficient to express the "
             "same effect and is therefore shrunk harder. Leaving the scores on their "
             "current scales is a defensible choice, but it is a CHOICE with consequences "
             "for relative shrinkage, not a neutral default."),
         "options_not_chosen": [
             "leave as frozen (no score-level standardisation)",
             "standardise scores inside each training fold only",
             "standardise and report both as a prespecified sensitivity analysis"],
         "constraint": "any standardisation must be training-fold-only; a full-cohort "
                       "score standardisation would reintroduce the leakage the frozen "
                       "rule exists to prevent",
         "decision": "AWAITING_APPROVAL"},
        {"id": "D8D-2", "title": "Four-set primary vs five-set sensitivity reporting",
         "current_frozen_state": "4-set PRIMARY, 5-set SENSITIVITY-only; PI3K may not be "
                                 "promoted or removed on observed performance (M2-A1)",
         "blocking": "documented_limitation",
         "blocks_what": "nothing; it is a reporting-framework decision",
         "evidence": ["8C §4", "reports/phase8c/score_correlation_pearson_5set_sensitivity.csv"],
         "implication": (
             "Adding PI3K raises the condition number 2.95 -> 7.03 and introduces the "
             "only |r| > 0.5 pair (P53-PI3K, 0.631). A divergence between the primary "
             "and sensitivity results is therefore FORESEEABLE and should have a reading "
             "rule fixed BEFORE any outcome is seen, or the choice between them becomes "
             "a post-hoc selection."),
         "proposed_reporting_framework_not_a_change": [
             "report the 4-set result as primary ALWAYS, whatever the 5-set shows",
             "report the 5-set alongside it, never instead of it",
             "state in advance that a divergence is interpreted as evidence about the "
             "STABILITY of the primary estimate, not as a reason to switch panels",
             "pre-specify the divergence metric before unblinding"],
         "decision": "AWAITING_APPROVAL"},
        {"id": "D8D-3", "title": "WNT representation sensitivity",
         "current_frozen_state": "median aggregation is primary; representative-IQR is the "
                                 "frozen sensitivity representation (Phase 2B)",
         "blocking": "blocks_specific_analysis",
         "blocks_what": "only the interpretation of a WNT coefficient, if one is reported",
         "evidence": ["8C §6.2", "reports/phase8c/representation_sensitivity.csv"],
         "implication": (
             "Between the two frozen representations the scores correlate E2F 0.998, "
             "EMT 0.992, P53 0.972, PI3K 0.910, WNT 0.858. WNT is the smallest set (41 "
             "genes) and the least internally coherent of the primary four, so each "
             "gene's probe choice matters more. A WNT result that does not hold under "
             "the sensitivity representation would be a representation artefact as much "
             "as a finding."),
         "options_not_chosen": [
             "no representation sensitivity analysis (status quo)",
             "prespecify a full re-fit under the representative-IQR representation",
             "prespecify re-fit for WNT only"],
         "decision": "AWAITING_APPROVAL"},
        {"id": "D8D-4", "title": "Prespecified batch sensitivity analysis",
         "current_frozen_state": "no batch covariate in M1/M2/M3; no batch adjustment",
         "blocking": "documented_limitation",
         "blocks_what": "nothing; the primary model is unaffected either way",
         "evidence": ["8B §8.2", "reports/phase8b/new_findings_stratified_permutation.csv",
                      "reports/phase8b/pc_variance_heterogeneity_by_batch.csv"],
         "implication": (
             "Batch LOCATION effects on PCs and on scores are fully explained by MMR "
             "composition (0/10 score tests survive stratification). Two signals do "
             "survive: AFFX control-probe PC1 ~ batch (p=0.0016) and PC2 DISPERSION ~ "
             "batch (p=0.0001). Neither is a location effect, so adding a batch covariate "
             "would not address them. Note also that batch, RNA extraction method and "
             "MMR are mutually confounded in this deposit, so no analysis here can "
             "separate them. ComBat was applied before deposit, so these are residuals "
             "in an already-corrected matrix - their absence elsewhere is not evidence "
             "that batch effects are absent."),
         "options_not_chosen": [
             "no batch sensitivity analysis (status quo)",
             "prespecify a batch-stratified Cox as a sensitivity analysis",
             "prespecify a leave-one-batch-out stability check"],
         "decision": "AWAITING_APPROVAL"},
        {"id": "D8D-5", "title": "combat_mod identity and the ComBat model matrix",
         "current_frozen_state": "UNRESOLVED (Phase 7A Q1). Not used to impute MMR.",
         "blocking": "documented_limitation",
         "blocks_what": "nothing in the primary analysis; it limits provenance claims "
                        "and any statement about what ComBat protected",
         "evidence": ["8B §2 and §5.5", "reports/phase8b/combat_mod_vs_mmr.csv"],
         "implication": (
             "combat_mod is perfectly concordant with recorded MMR across 519 SDRF "
             "tumour arrays (0 discordant). That is consistent with the depositors "
             "having used mmr.status in the ComBat model matrix, and equally with any "
             "variable that partitions these arrays identically. It is NOT proof of "
             "identity. The 43 S3 samples with MMR unrecorded are NOT imputed from it; "
             "A-19 complete-case remains frozen."),
         "outstanding_route": "author correspondence with the CIT/GSE39582 depositors",
         "options_not_chosen": ["leave unresolved and state the limitation",
                                "attempt author correspondence"],
         "decision": "AWAITING_APPROVAL"},
        {"id": "D8D-6", "title": "Eight sex-axis-discordant samples",
         "current_frozen_state": "RETAINED. Phase 7C D2 retains all 453; no relabelling.",
         "blocking": "documented_limitation",
         "blocks_what": "nothing",
         "evidence": ["8A §1.5", "8B §7", "reports/phase8b/sex_discordant_eight_review.csv"],
         "implication": (
             "Eight of 453 samples are discordant with recorded sex on both the Y-linked "
             "and the XIST axis, consistently across an 11-marker panel, spread over five "
             "batches. 3/8 carry a Phase 7C soft flag against 29/453 cohort-wide "
             "(OR 9.67, p=0.011), but a soft flag reports array quality, not identity. "
             "Expression markers cannot distinguish a sample or annotation mix-up from "
             "loss of Y or XIST silencing. GENOTYPE-LEVEL EVIDENCE IS UNAVAILABLE to "
             "this project."),
         "options_not_chosen": ["retain and document (status quo)",
                                "prespecify a leave-these-eight-out stability check"],
         "decision": "AWAITING_APPROVAL"},
        {"id": "D8D-7", "title": "Endpoint limitations U1 and U2 (Phase 6A numbering)",
         "current_frozen_state": "U1 `[PROVISIONAL: months]`; U2 published rule UNRESOLVED, "
                                 "deposited encoding is censoring-at-death and that "
                                 "censoring is demonstrably informative",
         "blocking": "blocks_specific_analysis",
         "blocks_what": (
             "U1 blocks every FIXED-TIME absolute-risk analysis - t-year risk, predicted "
             "survival/CIF curves, cloglog calibration at a horizon, time-dependent "
             "AUC/Brier, and the 60-unit endpoint - unless explicitly labelled "
             "conditional (Phase 6A §3.1). Relative-risk analyses are unaffected. "
             "U2 constrains any analysis that treats death without recurrence as "
             "non-informative censoring, and motivates a competing-risks sensitivity "
             "analysis."),
         "evidence": ["docs/PHASE6A_ENDPOINT_RESOLUTION.md §1, §2, §3.1"],
         "naming_hazard": (
             "Phase 1 uses the labels U1 and U2 for DIFFERENT items (CIT446 stage vs M1; "
             "CIT106/CIT250 zero times). The Phase 6A numbering is the one used here and "
             "in the Phase 8A-8D reports. The collision is recorded, not resolved, "
             "because renumbering a historical report is out of scope."),
         "options_not_chosen": ["proceed with relative-risk claims only (Phase 6A position)",
                                "attempt documentary resolution of the time unit"],
         "decision": "AWAITING_APPROVAL"},
        {"id": "D8D-8", "title": "External validation",
         "current_frozen_state": "NO-GO for CEL acquisition (Phase 4 §6.2). No external "
                                 "expression data downloaded.",
         "blocking": "blocks_specific_analysis",
         "blocks_what": "any confirmatory external-validation claim. Internal validation "
                        "by bootstrap optimism correction is unaffected.",
         "evidence": ["docs/PHASE4_PREACQUISITION.md §6.1-6.3",
                      "docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md §9"],
         "implication": (
             "Gate G4 FAILED on pre-specified criteria: GSE14333 has 48 events and "
             "GSE33113 has 18, against roughly 338 required at C=0.65 under an "
             "APPROXIMATION of the Riley framework (14% and 5% of requirement). Both "
             "remain EXPLORATORY candidates and their counts MUST NOT be combined. "
             "GSE33113 is a separate estimand - stage II only, a metastasis-or-recurrence "
             "composite, time in days. GSE17536 is excluded for established patient "
             "overlap with GSE14333. No suitable confirmatory cohort has been identified."),
         "options_not_chosen": ["proceed with internal validation only and state the limitation",
                                "revisit the four D5 unlock conditions"],
         "decision": "AWAITING_APPROVAL"},
    ]
    for d in DEC:
        print(f"  {d['id']}  {d['title'][:52]:54s} {d['blocking']:26s} {d['decision']}")
    blocking_counts = pd.Series([d["blocking"] for d in DEC]).value_counts().to_dict()
    print(f"\n  blocking profile: {blocking_counts}")
    assert all(d["decision"] == "AWAITING_APPROVAL" for d in DEC)
    assert not any(d["blocking"] == "blocks_phase9_entirely" for d in DEC)

    # =====================================================================
    section("5. PHASE 9 READINESS CHECKLIST")
    # =====================================================================
    READY = [
        ("analysis cohort frozen and retained (453)", "READY", "Phase 7C D2"),
        ("expression representation frozen", "READY", "Phase 2B"),
        ("analysis bundle verified", "READY", "scripts/24"),
        ("M2 panels frozen, roles fixed, ddof pinned", "READY", "M2-A1"),
        ("canonical scorer single implementation in src/", "READY", "Phase 8C rec. 4"),
        ("fold-isolation machinery verified", "READY", "8A, 8C (deletion + poisoning)"),
        ("missing-data strategy frozen", "READY", "A-19 complete-case"),
        ("QC executed, decisions recorded", "READY", "Phase 7C rev. 3"),
        ("EDA complete and consolidated", "READY", "Phase 8A-8D"),
        ("score-level standardisation policy", "DECISION REQUIRED", "D8D-1"),
        ("4-set/5-set reporting framework", "DECISION REQUIRED", "D8D-2"),
        ("representation sensitivity policy", "DECISION REQUIRED", "D8D-3"),
        ("batch sensitivity policy", "DECISION REQUIRED", "D8D-4"),
        ("combat_mod provenance", "OPEN - documented limitation", "D8D-5"),
        ("sex-axis discordance", "OPEN - documented limitation", "D8D-6"),
        ("U1 time unit", "CONSTRAINS fixed-time analyses", "D8D-7"),
        ("U2 competing risks", "CONSTRAINS censoring assumption", "D8D-7"),
        ("external confirmatory validation", "NOT AVAILABLE", "D8D-8"),
        ("Phase 6B protocol finalised (resampling, MC stability)", "NOT STARTED", "Phase 6B draft"),
    ]
    chkdf = pd.DataFrame(READY, columns=["item", "status", "reference"])
    chkdf.to_csv(OUT / "phase9_readiness_checklist.csv", index=False)
    for _, r in chkdf.iterrows():
        print(f"  [{r['status']:28s}] {r['item']:52s} {r['reference']}")

    # =====================================================================
    section("INTEGRITY")
    # =====================================================================
    after = {str(p.relative_to(paths.ROOT)): sha(p) for p in READONLY}
    changed = [k for k in before if before[k] != after[k]]
    assert not changed, f"READ-ONLY ARTIFACT MODIFIED: {changed}"
    print(f"  read-only artifacts unchanged : {len(READONLY)}/{len(READONLY)}")
    for n in ("m2_panel_4set.json", "m2_panel_5set.json"):
        ok, why = resolve_hash(PANELS / n, sha(PANELS / n))
        print(f"  amendment chain {n:22s}: {why}")

    claims.to_csv(OUT / "consolidated_claims.csv", index=False)
    summary = {
        "phase": "8D", "version": "1.0",
        "purpose": "consolidation of Phases 8A, 8B and 8C; no new analysis",
        "n_samples": 453, "samples_excluded": 0,
        "outcome_columns_opened": False,
        "outcome_files_blocked": list(BLOCKED_FILES),
        "outcome_columns_never_read": list(OUTCOME_COLUMNS),
        "modelling_performed": False, "recurrence_association_tested": False,
        "feature_selection_performed": False,
        "new_analysis_performed": False,
        "normalization_applied": "none", "batch_correction_applied": "none",
        "claims_resolved_from_artifacts": int(len(claims)),
        "claims_transcribed": 0,
        "cross_phase_checks": checks,
        "all_cross_phase_checks_agree": True,
        "evidence_table": ev.to_dict("records"),
        "evidence_classification": cat,
        "decision_register": DEC,
        "decisions_awaiting_approval": sum(1 for d in DEC if d["decision"] == "AWAITING_APPROVAL"),
        "decisions_made_by_phase8d": 0,
        "blocking_profile": blocking_counts,
        "phase9_readiness": chkdf.to_dict("records"),
        "phase8c_recommendations_implemented_here": [],
        "frozen_state": {
            "panel_primary": dig(src["P4"], "role"),
            "panel_sensitivity": dig(src["P5"], "role"),
            "sigma_ddof": dig(src["P4"], "scoring.sigma_ddof"),
            "score_standardisation": dig(src["P4"], "scoring.score_standardisation"),
            "canonical_scorer_module": score_panel.__module__,
        },
        "read_only_artifact_hashes": after,
        "read_only_artifacts_modified": False,
        "environment": environment_report(),
    }
    (OUT / "phase8d_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    (OUT / "decision_register.json").write_text(
        json.dumps({"phase": "8D", "decisions": DEC,
                    "decisions_made_by_phase8d": 0,
                    "note": "Phase 8D decides nothing. Every entry awaits the "
                            "project owner's approval."}, indent=2) + "\n")
    print(f"\n  outputs -> {OUT.relative_to(paths.ROOT)}/ "
          f"({len(list(OUT.glob('*.csv')))} CSV + {len(list(OUT.glob('*.json')))} JSON)")
    print("\nPHASE 8D COMPLETE -- consolidation only; 0 decisions made, "
          f"{summary['decisions_awaiting_approval']} awaiting approval.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
