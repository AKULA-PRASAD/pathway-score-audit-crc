"""Phase 9: modeling protocol PLANNING -- v2.0 finalization.

No model is fitted; no outcome is read. v2.0 resolves D8D-1..D8D-4 in full,
reassesses the v1.0 sparse-batch pooling rule (and REVERSES it), adds the
penalisation-scope decision D9-1 that D8D-1 exposes, and finalises the age
spline and the validation details.

Emits the machine-readable half of the Phase 9 deliverables:
  * decision register for D8D-1..D8D-4 (recommendations, PENDING USER APPROVAL)
  * leakage-control matrix
  * model comparison plan
  * validation / resampling specification
  * modeling-readiness checklist
  * proposed M2 amendment text (NOT applied)

HARD CONSTRAINTS ENFORCED IN CODE:
  * endpoint.csv blocked at the pandas I/O layer. Documented event counts are
    RESOLVED BY REGEX from the governing reports, never recomputed from data.
  * No model is fitted. No outcome column is joined to any expression data.
  * All frozen records are read-only and hashed before and after.
  * Nothing is applied: every recommendation carries status PENDING_USER_APPROVAL.

Usage:
    /Users/kalki/coad_transportability/.venv/bin/python \
        /Users/kalki/coad_transportability/scripts/36_phase9_protocol_plan.py
"""

from __future__ import annotations

import hashlib
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

DOCS = paths.ROOT / "docs"
PANELS = paths.ROOT / "config" / "pathways"
CFG = paths.ROOT / "config" / "analysis"
ANALYSIS = paths.PROCESSED_DIR / "analysis"
OUT = paths.ROOT / "reports" / "phase9"

BLOCKED_FILES = ("endpoint.csv",)
PENDING = "PENDING_USER_APPROVAL"


def _guard_outcome_files() -> None:
    import pandas as _pd
    for fn in ("read_csv", "read_parquet", "read_table"):
        orig = getattr(_pd, fn)

        def blocked(path, *a, _orig=orig, **kw):
            if any(b in str(path) for b in BLOCKED_FILES):
                raise PermissionError(
                    f"PHASE 9 PLANNING IS OUTCOME-BLIND: refusing to read {path}")
            return _orig(path, *a, **kw)
        setattr(_pd, fn, blocked)


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# NOTE: reports/phase9/ is EXCLUDED -- it is this script's own output directory.
# Including it would make the integrity check guard the files the script writes,
# which both misstates what is being protected and fails on any content change.
READONLY = sorted(
    [p for p in (paths.ROOT / "reports").rglob("*")
     if p.is_file() and "phase9" not in p.parts] +
    # PHASE9_MODELING_PROTOCOL_DRAFT.md is EXCLUDED: it is this phase's own
    # deliverable, not a prior-phase record to be protected.
    [p for p in DOCS.glob("PHASE*.md") if "PHASE9" not in p.name] +
    [PANELS / "m2_panel_4set.json", PANELS / "m2_panel_5set.json",
     PANELS / "M2_FREEZE_RECORD.json", CFG / "A19_FREEZE_RECORD.json",
     CFG / "QC_DECISIONS_RECORD.json", CFG / "QC_SPECIFICATION.json",
     paths.ROOT / "src" / "pathway_scoring.py",
     ANALYSIS / "cohort_membership.csv", ANALYSIS / "covariates.csv",
     ANALYSIS / "expression_median_S0.parquet",
     ANALYSIS / "expression_repiqr_disc443_S0.parquet"])


def section(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)


def cite(relpath: str, pattern: str, label: str):
    """Resolve a documented figure from its governing report by regex.

    Phase 9 must not recompute an event count -- that is outcome data. Figures
    already established and approved in earlier phases are quoted FROM THOSE
    REPORTS, with the match recorded so the citation is checkable.
    """
    txt = (paths.ROOT / relpath).read_text()
    m = re.search(pattern, txt)
    if not m:
        raise LookupError(f"{label}: pattern not found in {relpath}")
    return {"label": label, "value": m.group(1), "source": relpath,
            "matched_text": m.group(0)[:120]}


def main() -> int:
    _guard_outcome_files()
    OUT.mkdir(parents=True, exist_ok=True)
    print("=" * 78)
    print("PHASE 9 -- MODELING PROTOCOL PLANNING (no model fitted, no outcome read)")
    print("=" * 78)

    before = {str(p.relative_to(paths.ROOT)): sha(p) for p in READONLY}

    # =====================================================================
    section("0. GOVERNING RECORDS -- filename verification")
    # =====================================================================
    expected = {
        "phase8d_report": "docs/PHASE8D_CONSOLIDATED_EDA_REVIEW.md",
        "phase8d_register": "reports/phase8d/decision_register.json",
        "phase8d_checklist": "reports/phase8d/phase9_readiness_checklist.csv",
        "phase6b_protocol": "docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md",
        "phase6a_endpoint": "docs/PHASE6A_ENDPOINT_RESOLUTION.md",
        "phase5_precision": "docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md",
        "phase6c_panel": "docs/PHASE6C_PATHWAY_PANEL_REVIEW.md",
        "phase6d_missing": "docs/PHASE6D_MISSINGNESS_AUDIT.md",
        "panel_4set": "config/pathways/m2_panel_4set.json",
        "panel_5set": "config/pathways/m2_panel_5set.json",
        "m2_record": "config/pathways/M2_FREEZE_RECORD.json",
        "scorer": "src/pathway_scoring.py",
    }
    missing = {k: v for k, v in expected.items() if not (paths.ROOT / v).exists()}
    assert not missing, f"governing record missing: {missing}"
    print(f"  {len(expected)} governing records verified present")
    print("  NOTE: the brief named docs/PHASE6C_HALLMARK_PANEL_REVIEW.md, which does")
    print("        NOT exist. The actual Phase 6C report is")
    print("        docs/PHASE6C_PATHWAY_PANEL_REVIEW.md and is the record used here.")

    # =====================================================================
    section("1. DOCUMENTED FIGURES -- cited from reports, never recomputed")
    # =====================================================================
    CITED = [
        cite("docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md",
             r"Primary cohort \| \*\*S3\*\* \((\d+) samples, \d+ events\)", "s3_n"),
        cite("docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md",
             r"Primary cohort \| \*\*S3\*\* \(\d+ samples, (\d+) events\)", "s3_events"),
        cite("docs/PHASE6A_ENDPOINT_RESOLUTION.md",
             r"\*\*S-1\*\* \| CIF \+ Fine–Gray subdistribution \| \d+ rec / (\d+) competing",
             "competing_events"),
        cite("docs/PHASE6A_ENDPOINT_RESOLUTION.md",
             r"\*\*S-2\*\* \| Composite DFS \| \*\*(\d+)\*\*", "composite_dfs_events"),
        cite("docs/PHASE6A_ENDPOINT_RESOLUTION.md",
             r"\*\*S-3\*\* \| 60-unit administrative horizon \| (\d+)", "s3_horizon_events"),
    ]
    for c in CITED:
        print(f"  {c['label']:22s} = {c['value']:>5s}   <- {c['source']}")
    S3_N = int(next(c["value"] for c in CITED if c["label"] == "s3_n"))
    EV = int(next(c["value"] for c in CITED if c["label"] == "s3_events"))
    assert S3_N == 453

    # complete-case counts come from the Phase 6D audit, not recomputed
    cc = {m["model"]: m for m in
          json.loads((paths.ROOT / "reports/phase6d/definitive_audit.json").read_text())["models"]}
    a19 = json.loads((CFG / "A19_FREEZE_RECORD.json").read_text())
    dropped = [e["sample_id"] for e in a19["excluded_observations"]]
    # Both count sets, side by side -- the v2.0 wording was easy to misread.
    _pb = pd.read_csv(ANALYSIS / "partition_and_batch.csv").set_index("sample_id")
    _mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    _s3 = _mem.loc[_mem.in_S3.astype(bool), "sample_id"].tolist()
    _cov = pd.read_csv(ANALYSIS / "covariates.csv", keep_default_na=False,
                       dtype=str).set_index("sample_id")
    _cc = [s for s in _s3 if all(_cov.loc[s, v] != "" for v in
                                 ["age.at.diagnosis_(year)", "sex", "tnm.stage",
                                  "tumor.location", "chemotherapy.adjuvant"])]
    BATCH_COUNTS = {
        "s3_cohort_453": {k: int(v) for k, v in
                          _pb.loc[_s3, "combat_batch"].value_counts().sort_index().items()},
        "complete_case_452": {k: int(v) for k, v in
                              _pb.loc[_cc, "combat_batch"].value_counts().sort_index().items()},
        "difference": "one sample, GSM972262 in batch05, dropped by A-19 for "
                      "missing chemotherapy.adjuvant",
        "join_integrity": "0 missing, 0 duplicated, 0 blank, 0 ambiguous; the "
                          "bundle agrees with the ArrayExpress SDRF for all 453",
        "cohort_unchanged": "all 453 S3 samples remain in the cohort"}
    print(f"\n  batch counts, BOTH sets (v2.1 clarification):")
    print(f"    {'batch':9s} {'S3 cohort (453)':>16s} {'complete-case (452)':>21s}")
    for k in BATCH_COUNTS["s3_cohort_453"]:
        print(f"    {k:9s} {BATCH_COUNTS['s3_cohort_453'][k]:16d} "
              f"{BATCH_COUNTS['complete_case_452'].get(k, 0):21d}")

    print(f"\n  complete-case model sizes (Phase 6D audit, A-19 complete-case):")
    for k, v in cc.items():
        print(f"    {k:42s} p={v['n_predictors']:2d}  n={v['complete_cases']:3d}  lost={v['lost']}")
    print(f"    A-19 drops: {dropped}")
    print(f"\n  EVENTS IN THE COMPLETE-CASE SET ARE NOT DETERMINED HERE.")
    print(f"    S3 carries {EV} events (cited). A-19 drops {len(dropped)} sample(s);")
    print(f"    whether that sample is an event is OUTCOME DATA and was not read.")
    print(f"    Complete-case events therefore lie in [{EV-len(dropped)}, {EV}].")
    epv_lo = (EV - len(dropped)) / cc["M2 primary (4-set) + M1"]["n_predictors"]
    epv_hi = EV / cc["M2 primary (4-set) + M1"]["n_predictors"]
    print(f"    M2 primary EPV therefore lies in [{epv_lo:.1f}, {epv_hi:.1f}].")

    # =====================================================================
    section("2. DECISION REGISTER -- recommendations, NOTHING APPLIED")
    # =====================================================================
    p4 = json.loads((PANELS / "m2_panel_4set.json").read_text())
    current_std = p4["scoring"]["score_standardisation"]
    # M2-A2 is applied, so the panel's current value is the NEW text. The
    # amendment's "from" must be the superseded wording, taken from history.
    _m2a2 = next((a for a in p4.get("amendment_history", [])
                  if a["amendment_id"] == "M2-A2"), None)
    superseded_std = (_m2a2["superseded_score_standardisation"] if _m2a2
                      else current_std)

    DEC = [
        {
            "id": "D8D-1", "title": "Score-level standardisation before penalised Cox",
            "decision_status": "APPROVED 2026-11-03 (Option B) -- amendment M2-A2 APPLIED",
            "current_frozen_state": current_std,
            "frozen_state_before_m2a2": superseded_std,
            "options_considered": [
                "A. leave scores on their current scales (status quo)",
                "B. standardise scores using TRAINING-FOLD mean and SD only",
                "C. standardise and report both as a prespecified sensitivity analysis"],
            "recommendation": "B",
            "recommendation_text":
                "Standardise the pathway scores to unit variance using training-fold "
                "parameters only, inside every resampling boundary, and retain the "
                "unstandardised scores unchanged as the auditable source features.",
            "rationale": [
                "The frozen justification -- 'penalisation handles scale' -- is not "
                "correct as stated. Ridge and elastic-net penalties are NOT "
                "scale-invariant: the penalty lambda*sum(beta^2) is applied per unit "
                "coefficient, so a predictor with smaller dispersion needs a larger "
                "coefficient to express the same effect and is therefore shrunk harder.",
                "Phase 8C measured a 3x dispersion range across the five scores "
                "(SD 0.554 E2F to 0.185 PI3K). Under the current state the penalty "
                "would treat them unequally for a reason unrelated to their content.",
                "Standardised coefficients are interpretable on a common footing -- "
                "log-hazard per 1 SD of score -- which is what a cross-pathway "
                "comparison requires.",
                "The coherence information is NOT lost: implied within-set gene "
                "correlation is already recorded in reports/phase8c/score_distributions.csv "
                "and is reported alongside, where it belongs, rather than being "
                "smuggled into the penalty.",
            ],
            "limitations_of_the_recommendation": [
                "Standardising deliberately removes the dispersion difference from the "
                "coefficient scale. That difference is real information about set "
                "coherence; it must be reported separately or it is simply lost.",
                "It adds one more quantity estimated inside each fold, which adds "
                "Monte-Carlo variance and must be covered by the MC-1..MC-4 checks.",
                "Any comparison with an unstandardised fit is not like-for-like and "
                "must not be presented as one.",
                "Low dispersion is NOT low biological relevance. Standardisation must "
                "not be read as a claim that the five scores are equally informative; "
                "no outcome association has been tested.",
            ],
            "scaling_guard_rules": {
                "zero_variance_score_in_a_training_fold":
                    "if the training-fold SD of a score is 0 (or < 1e-12), the "
                    "standardisation is undefined. That score is DROPPED for that "
                    "fold and the drop is recorded per fold -- the same convention "
                    "the frozen gene-level rule already uses for zero-variance "
                    "genes. It is never silently set to 0 or 1.",
                "non_finite_score":
                    "a non-finite training-fold mean or SD aborts that replicate, "
                    "which is recorded as a non-estimable replicate (see the "
                    "validation specification) rather than being imputed away.",
                "expected_frequency":
                    "neither condition is expected: Phase 8C observed 0 "
                    "zero-variance genes and 0 non-finite values across 453 "
                    "samples. The rules exist so the behaviour is defined, not "
                    "because the conditions are anticipated.",
                "recorded_in": "per-fold diagnostics, reported with the model"},
            "clinical_predictor_scaling": {
                "question": "does standardising the pathway scores settle the "
                            "scaling of the clinical predictors?",
                "answer": "NO. It settles only the pathway block.",
                "detail":
                    "The clinical block is four binaries plus a 2-column age "
                    "spline basis on the YEAR scale. Those columns have wildly "
                    "different dispersions from each other and from a "
                    "standardised score, so if they are penalised they inherit "
                    "exactly the scale-dependence D8D-1 is fixing -- in the "
                    "opposite direction, since a large-valued basis column "
                    "carries a small coefficient and is barely penalised.",
                "resolution": "this is NOT resolved inside D8D-1. It is raised as "
                              "a separate decision D9-1 (penalisation scope), "
                              "because the clean answer is to leave the clinical "
                              "block UNPENALISED rather than to standardise it."},
            "coefficient_interpretation_after_standardisation": {
                "scale": "log cause-specific hazard ratio per 1 training-fold SD "
                         "of the pathway score",
                "comparability": "coefficients become comparable ACROSS pathways "
                                 "because each is expressed per its own SD",
                "what_it_does_NOT_mean":
                    "a larger standardised coefficient does NOT mean a pathway is "
                    "more biologically important. Standardisation removes the "
                    "dispersion difference from the coefficient scale, so the "
                    "coefficient no longer carries information about how "
                    "coordinated the gene set is.",
                "dispersion_vs_relevance":
                    "score dispersion is a joint property of the gene set's "
                    "internal coherence and of the averaging transform "
                    "(SD = 1/sqrt(k) if members were independent). A "
                    "low-dispersion score such as P53 (implied inter-gene "
                    "r = 0.033) is a weakly coordinated AVERAGE, which is a "
                    "statement about co-expression in this tissue under this "
                    "preprocessing -- NOT a statement about biological relevance "
                    "and NOT a statement about recurrence, which has not been "
                    "tested. Phase 8C's coherence table is reported alongside the "
                    "coefficients precisely so the two are not conflated."},
            "leakage_constraints_if_approved": [
                "mean and SD estimated on TRAINING samples only",
                "training-derived parameters applied UNCHANGED to held-out samples",
                "scaling recomputed independently inside EVERY resampling, "
                "bootstrap and inner-CV tuning operation -- including the inner "
                "tuning loop, not only the outer split",
                "held-out samples never influence the scaling parameters",
                "unstandardised scores preserved UNCHANGED as the auditable "
                "source features, written to disk and reported",
                "zero-variance and non-finite conditions handled by the explicit "
                "rules above and recorded per fold",
            ],
            "proposed_amendment_id": "M2-A2",
            "amendment_status": "APPROVED 2026-11-03 AND APPLIED",
            "amendment_applied": True,
            "applied_amendment": {
                "target_files": ["config/pathways/m2_panel_4set.json",
                                 "config/pathways/m2_panel_5set.json",
                                 "config/pathways/M2_FREEZE_RECORD.json"],
                "field": "scoring.score_standardisation",
                "from": superseded_std,
                "current_on_disk": current_std,
                "to": ("TRAINING-FOLD ONLY. Each pathway score is centred and scaled to "
                       "unit variance using the mean and sample SD (ddof=1) of the "
                       "TRAINING samples of that fold, applied unchanged to held-out "
                       "samples. Recomputed inside every resampling replicate. The "
                       "unstandardised scores are preserved as the auditable source "
                       "features. This is a SECOND standardisation, applied to the "
                       "scores, and does not alter the frozen gene-level "
                       "standardisation (sigma_ddof=1) that produces them."),
                "membership_changed": False, "roles_changed": False,
                "outcome_data_used": False,
                "permitted_under": "modification_policy: non-outcome correction, "
                                   "recorded with reason and timestamp, prior version "
                                   "retained, re-freeze",
            },
            "status": PENDING,
        },
        {
            "id": "D8D-2", "title": "Four-set primary vs five-set sensitivity interpretation",
            "current_frozen_state": "4-set PRIMARY, 5-set SENSITIVITY-only; PI3K may not "
                                    "be promoted or removed on observed performance (M2-A1)",
            "recommendation": "Adopt the reporting framework below. NO role change.",
            "comparisons_to_be_reported": [
                "M1 (clinical baseline) -- the fixed benchmark",
                "M2 primary (M1 + 4 scores) vs M1: optimism-corrected delta-C with CI "
                "and Monte-Carlo SE -- THE primary comparison",
                "M2 sensitivity (M1 + 5 scores) vs M1: same quantities",
                "M2 sensitivity vs M2 primary: reported as a STABILITY check, never as "
                "a selection contest",
                "M3 exploratory vs M1: reported with its exploratory label attached",
            ],
            "interpretation_rules_fixed_in_advance": [
                "The 4-set result is reported as primary ALWAYS, whatever the 5-set shows.",
                "The 5-set is reported ALONGSIDE, never instead of, the 4-set.",
                "A delta-C between the panels that lies within its Monte-Carlo SE is "
                "declared 'no evidence of a difference' -- not a tie to be broken by a "
                "second metric.",
                "A metric improvement in the 5-set is NOT grounds to switch panels. "
                "Added flexibility and the extra parameter can improve apparent fit; "
                "that is a property of the comparison, not evidence about PI3K.",
                "Coefficient changes for P53 and PI3K are reported together with their "
                "CIs and variance inflation, because r = 0.631 between them makes "
                "individual coefficient instability FORESEEABLE.",
                "PI3K's descriptive redundancy (R2 = 0.540 on the primary four, Phase "
                "8C) is stated whenever its coefficient is discussed.",
                "Divergence is interpreted as evidence about the STABILITY of the "
                "primary estimate, not as a reason to prefer either panel.",
            ],
            "prespecified_divergence_metric": "optimism-corrected delta-C (5-set minus "
                                              "4-set) with its bootstrap CI and "
                                              "Monte-Carlo SE, fixed before unblinding",
            "approval": {
                "approved": True, "date": "2026-11-03",
                "approved_by": "project owner",
                "scope": "panel reporting hierarchy only; no model fitting, no "
                         "outcome analysis authorised",
                "original_proposal_preserved_above": True,
                "basis": ["Phase 8D decision register D8D-2",
                          "Phase 9 protocol v2.0 section 3.2",
                          "Phase 9 Final Audit (v2.1)"]},
            "formalised_decision": {
                "primary_panel": {
                    "panel": "4-set", "role": "PRIMARY (pre-specified)",
                    "members": ["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                "HALLMARK_E2F_TARGETS",
                                "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                                "HALLMARK_P53_PATHWAY"]},
                "sensitivity_panel": {
                    "panel": "5-set", "role": "SENSITIVITY-ONLY",
                    "adds": "HALLMARK_PI3K_AKT_MTOR_SIGNALING"},
                "unchanged_by_this_approval": [
                    "pathway membership", "gene mappings",
                    "scoring definitions (sigma_ddof=1, M2-A1; score "
                    "standardisation, M2-A2)",
                    "the frozen panel specification files"],
                "prohibitions": [
                    "the 5-set panel MUST NOT be promoted to primary on the basis of "
                    "model performance or statistical significance",
                    "the 4-set analysis IS the pre-specified primary result and is "
                    "reported as such REGARDLESS of which panel performs better",
                    "the 5-set analysis is reported as a SENSITIVITY analysis, never "
                    "as the headline result",
                    "pathway-score dispersion or correlation MUST NOT be interpreted "
                    "as proof of biological activity or clinical relevance"],
                "mandatory_reporting_with_the_5_set_analysis": {
                    "p53_pi3k_pearson_r": 0.6314,
                    "collinearity_diagnostics": {
                        "condition_number_4set": 2.949,
                        "condition_number_5set": 7.032,
                        "min_eigenvalue_4set": 0.4811,
                        "min_eigenvalue_5set": 0.2777,
                        "n_pairs_abs_r_gt_0.5_4set": 0,
                        "n_pairs_abs_r_gt_0.5_5set": 1,
                        "pi3k_r2_on_the_primary_four": 0.5398,
                        "also_report": "variance inflation factors and the "
                                       "coefficient CIs for P53 and PI3K TOGETHER"},
                    "source": "reports/phase8c/score_correlation_pearson_5set_"
                              "sensitivity.csv; reports/phase8c/phase8c_summary.json"},
                "why_pi3k_is_sensitivity_only": {
                    "reason": "its correlation with P53 (r = 0.6314) is the only "
                              "|r| > 0.5 pair in the 5-set panel and raises the "
                              "condition number from 2.949 to 7.032, which makes "
                              "individual coefficient instability FORESEEABLE",
                    "explicitly_NOT_the_reason":
                        "this is NOT evidence that PI3K/AKT/mTOR is biologically "
                        "irrelevant. It is a statement about the CONDITIONING of the "
                        "score matrix in this cohort under this preprocessing, and "
                        "nothing more. No outcome association has been tested for "
                        "any pathway."},
                "analysis_hierarchy_vs_observed_performance": {
                    "rule": "the reporting protocol distinguishes the PRE-SPECIFIED "
                            "ANALYSIS HIERARCHY from OBSERVED MODEL PERFORMANCE. The "
                            "hierarchy was fixed before any outcome was seen and does "
                            "not move; performance is a result reported within it.",
                    "operationally": [
                        "primary result = 4-set, always, stated first",
                        "sensitivity result = 5-set, reported alongside and labelled",
                        "a delta-C between panels within its Monte-Carlo SE is "
                        "reported as 'no evidence of a difference', not as a tie to "
                        "be broken by another metric",
                        "divergence is reported transparently and read as evidence "
                        "about the STABILITY of the primary estimate"]}},
            "status": "APPROVED 2026-11-03 -- reporting protocol formalised",
        },
        {
            "id": "D8D-3", "title": "WNT representation sensitivity",
            "current_frozen_state": "median aggregation PRIMARY; representative-IQR is "
                                    "the frozen sensitivity representation (A-05), and "
                                    "Phase 6A already registers estimand S-7 as a "
                                    "representative-IQR re-fit",
            "observed_evidence": {
                "wnt_median_vs_repiqr_pearson": 0.8577,
                "others": {"E2F_TARGETS": 0.9975, "EMT": 0.9920, "P53": 0.9720,
                           "PI3K": 0.9096},
                "source": "reports/phase8c/representation_sensitivity.csv"},
            "recommendation":
                "Do NOT create a new policy. Use the ALREADY-REGISTERED estimand S-7 "
                "(representative-IQR re-fit) and attach to it a prespecified, "
                "WNT-focused comparison.",
            "comparison_specification": {
                "unit_of_comparison": "a complete re-fit of M2 primary under the "
                                      "representative-IQR representation; the same "
                                      "pipeline, the same folds, the same seeds",
                "quantities_compared": [
                    "per-pathway standardised coefficient and its CI, pathway by pathway",
                    "optimism-corrected C-index and delta-C vs M1",
                    "rank correlation of the linear predictor between representations",
                    "the WNT coefficient specifically, since WNT is the only score "
                    "whose two representations correlate below 0.9"],
                "interpretation": [
                    "Agreement means the result does not depend on the probe-to-gene "
                    "rule. It is NOT independent confirmation -- it is the same samples, "
                    "the same genes and the same deposited matrix, summarised differently.",
                    "Disagreement confined to WNT is read as representation sensitivity "
                    "of a 41-gene, low-coherence set (implied inter-gene r = 0.041), "
                    "not as biology.",
                    "Disagreement across several pathways would call the representation "
                    "choice itself into question and is reported as such."],
                "selection_guard": "The MEDIAN representation is the primary result and "
                                   "is reported as primary REGARDLESS of which "
                                   "representation is more favourable. The "
                                   "representative-IQR re-fit is reported in full "
                                   "whatever it shows. Neither may be chosen after "
                                   "seeing outcomes."},
            "approval": {
                "approved": True, "date": "2026-11-03",
                "approved_by": "project owner",
                "approved_scope": "median gene representation PRIMARY vs "
                                  "representative-IQR SENSITIVITY -- broader than the "
                                  "original WNT-focused proposal, which is preserved "
                                  "above unchanged",
                "original_proposal_preserved_above": True,
                "basis": ["Phase 8D decision register D8D-3",
                          "Phase 9 protocol v2.0 section 3.3",
                          "Phase 6A estimand S-7", "Phase 2B representation freeze"]},
            "formalised_decision": {
                "primary_representation": {
                    "name": "median probe aggregation", "role": "PRIMARY",
                    "frozen_in": "Phase 2B (A-04)",
                    "artifact": "data/processed/analysis/expression_median_S0.parquet"},
                "sensitivity_representation": {
                    "name": "representative-IQR (discovery_443 selection)",
                    "role": "SENSITIVITY-ONLY", "frozen_in": "Phase 2B (A-05)",
                    "registered_estimand": "S-7",
                    "artifact": "data/processed/analysis/"
                                "expression_repiqr_disc443_S0.parquet"},
                "prohibitions": [
                    "the primary representation MUST NOT be replaced on the basis of "
                    "observed model performance",
                    "the representation correlations MUST NOT be interpreted as "
                    "evidence that one representation is universally superior",
                    "neither representation may be chosen after seeing outcomes"],
                "mandatory_reporting_when_the_analyses_are_authorised": [
                    "whether key conclusions are DIRECTIONALLY consistent across the "
                    "two representations",
                    "whether they are QUANTITATIVELY consistent",
                    "any differences in pathway scores, coefficients, uncertainty, "
                    "discrimination, or other relevant model outputs",
                    "the median-representation result as primary, reported in full "
                    "regardless of which representation looks more favourable"],
                "representation_agreement_phase8c": {
                    "metric": "Pearson correlation of the pathway score between the "
                              "two frozen representations, over the 453 S3 samples",
                    "HALLMARK_E2F_TARGETS": 0.9975,
                    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": 0.9920,
                    "HALLMARK_P53_PATHWAY": 0.9720,
                    "HALLMARK_PI3K_AKT_MTOR_SIGNALING": 0.9096,
                    "HALLMARK_WNT_BETA_CATENIN_SIGNALING": 0.8577,
                    "weakest_agreement": "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                    "weakest_value": 0.8577,
                    "source": "reports/phase8c/representation_sensitivity.csv",
                    "acknowledgement":
                        "WNT has the WEAKEST representation agreement of the five "
                        "scores (0.858) and is the only one below 0.9. It is the "
                        "smallest set (41 genes) and the least internally coherent of "
                        "the primary four (implied inter-gene r = 0.041), so each "
                        "gene's probe choice matters more. This is acknowledged "
                        "explicitly and carried into the reporting of any WNT result.",
                    "not_evidence_of_superiority":
                        "these correlations say how much the two summarisation rules "
                        "AGREE. They do NOT rank the representations, and a higher "
                        "correlation does not make either one more correct."},
                "robustness_is_not_performance": {
                    "statement": "representation ROBUSTNESS and PREDICTIVE PERFORMANCE "
                                 "are different properties and are reported "
                                 "separately, never substituted for one another.",
                    "robustness": "whether a conclusion survives changing the "
                                  "probe-to-gene summarisation rule. It is a "
                                  "stability property of the analysis pipeline.",
                    "performance": "how well a model discriminates or calibrates. It "
                                   "is an outcome-dependent quantity and NONE has "
                                   "been computed.",
                    "and_neither_is_replication":
                        "agreement between representations is NOT independent "
                        "replication: the same samples, the same genes and the same "
                        "deposited matrix, summarised differently."}},
            "status": "APPROVED 2026-11-03 -- reporting protocol formalised",
        },
        {
            "id": "D8D-4", "title": "Batch sensitivity analysis",
            "decision_status": "APPROVED 2026-11-03 -- design and D8D-4-R1 APPLIED",
            "current_frozen_state": "no batch covariate in M1/M2/M3; B-04 makes any "
                                    "batch correction a separate decision requiring approval",
            "evidence_considered": {
                "affx_pc1_vs_batch_stratified_p": 0.0016,
                "pc2_dispersion_vs_batch_stratified_p": 0.0001,
                "score_tests_surviving_mmr_stratification": "0 of 10",
                "batch_x_mmr_chi2_p": 1.607e-25,
                "source": ["reports/phase8b/new_findings_stratified_permutation.csv",
                           "reports/phase8c/score_vs_annotation_stratified.csv",
                           "reports/phase8b/batch_by_combat_mod_cells.csv"]},
            "recommendation":
                "A NARROW batch sensitivity analysis is justified, as a DIAGNOSTIC plus "
                "one stratified sensitivity model. NO batch covariate is added to M1 or "
                "M2, and batch adjustment never becomes the primary analysis.",
            "eligible_batch_variables": {
                "combat_batch": "ELIGIBLE -- 9 recovered ComBat labels; the only "
                                "variable that denotes processing group",
                "combat_mod": "NOT ELIGIBLE -- perfectly concordant with MMR status, so "
                              "adjusting for it would adjust away the dominant "
                              "biological axis of the cohort",
                "rna.extraction.method": "NOT ELIGIBLE as an adjustment -- almost "
                                         "completely confounded with batch "
                                         "(chi2 = 2265, df = 40), min cell 14; reported "
                                         "descriptively only",
                "dataset (discovery/validation)": "NOT ELIGIBLE -- Phase 1 established "
                                                  "these are not independent cohorts; "
                                                  "expression-level separation is "
                                                  "negligible (eps2 <= 0.005)"},
            "design": [
                "DIAGNOSTIC 1 -- leave-one-batch-out stability: refit M2 primary nine "
                "times, each omitting one ComBat batch, and report the spread of the "
                "standardised coefficients and of delta-C. This is a stability "
                "diagnostic, not an adjusted model, and estimates no batch effect. "
                "CAVEAT, pre-registered: the batches differ ~7-fold in size (14 to "
                "103), so omitting batch03 or batch04 removes ~23% of the cohort "
                "while omitting batch02 removes 3%. The spread therefore conflates "
                "batch identity with sample size, and the omitted n is reported "
                "beside every refit so the two are not confused.",
                "SENSITIVITY MODEL 1 -- batch-STRATIFIED Cox: M2 primary with the "
                "baseline hazard stratified by ComBat batch, estimating NO batch "
                "coefficients. Stratification allows different baseline hazards without "
                "forcing batch to absorb biological variation through a coefficient.",
                "Neither is run before the primary analysis, and neither can replace it.",
            ],
            "sparse_level_handling": {
                "v1_proposal_REVERSED":
                    "Phase 9 v1.0 proposed pooling strata below 20 samples into a "
                    "single 'small-batch' stratum. That proposal is WITHDRAWN. It "
                    "was wrong, and the reason is technical, not cosmetic.",
                "why_the_threshold_was_wrong":
                    "A stratified Cox model estimates NO parameters for its strata. "
                    "The partial likelihood factorises across strata, so a stratum "
                    "consumes no degrees of freedom. The usual 'too few events per "
                    "parameter' concern -- which motivates pooling for a COVARIATE "
                    "-- simply does not apply to STRATIFICATION. A stratum with few "
                    "events contributes little information; a stratum with zero "
                    "events contributes a likelihood factor of 1 and drops out. "
                    "Neither is an estimability failure. The 20-sample threshold "
                    "was a numerical reflex imported from the covariate case.",
                "what_pooling_would_cost":
                    "Pooling forces batch02, batch09 and batch10 -- which have "
                    "nothing in common except being small -- to SHARE one baseline "
                    "hazard. That is the precise assumption the sensitivity "
                    "analysis exists to test, so pooling would conceal exactly the "
                    "between-batch heterogeneity being probed, and would do so in "
                    "the three batches where heterogeneity is least observable. The "
                    "pooled stratum would also be an artefact with no technical "
                    "provenance: 'small' is not a processing group.",
                "options_compared": {
                    "1_keep_all_nine_strata":
                        "RECOMMENDED. No parameters are consumed; nothing becomes "
                        "non-estimable; no heterogeneity is concealed. Cost is "
                        "efficiency, not validity.",
                    "2_prespecified_pooling_on_counts_and_provenance":
                        "REJECTED. No technical provenance links the three small "
                        "batches, so any pooling would rest on sample count alone "
                        "-- the reflex just rejected. Pooling on provenance would "
                        "be defensible IF the SDRF documented shared processing; "
                        "it does not.",
                    "3_omit_the_stratified_model":
                        "HELD IN RESERVE as the pre-registered failure response, "
                        "not as the design."},
                "RECOMMENDED_RULE":
                    "All nine ComBat batches enter as individual strata. NO "
                    "pooling, by count or otherwise.",
                "stratum_sizes_complete_case": {
                    **BATCH_COUNTS["complete_case_452"],
                    "total": sum(BATCH_COUNTS["complete_case_452"].values())},
                "stratum_sizes_s3_cohort_for_contrast": {
                    **BATCH_COUNTS["s3_cohort_453"],
                    "total": sum(BATCH_COUNTS["s3_cohort_453"].values())},
                "estimability_depends_on_outcome_information":
                    "How much information the nine strata actually carry depends on "
                    "EVENTS PER STRATUM, which is outcome data and was NOT read "
                    "here. The DESIGN is fixed outcome-blind; its information "
                    "content can only be measured at execution. This is stated "
                    "rather than papered over.",
                "prespecified_failure_response": {
                    "amendment_id": "D8D-4-R1",
                    "supersedes": "the v2.0 rule 'fewer than 50% of events lie in "
                                  "strata that contribute to the partial "
                                  "likelihood', WITHDRAWN by the Phase 9 final "
                                  "audit as logically near-vacuous (every event "
                                  "lies in a stratum containing an event -- itself) "
                                  "and statistically unfounded (50% has no basis, "
                                  "and labels a PRECISION statement as "
                                  "'estimability')",
                    "principle":
                        "declare non-estimability on MATHEMATICAL and NUMERICAL "
                        "grounds only -- criteria whose DEFINITIONS are "
                        "outcome-blind -- and report information loss as a "
                        "continuous diagnostic rather than gating on an invented "
                        "threshold",
                    "declared_non_estimable_if_and_only_if": [
                        "the stratified fit fails to converge under the declared "
                        "criteria (log-likelihood change < 1e-9, max 100 "
                        "iterations), OR",
                        "the observed information matrix for beta is singular or "
                        "not positive definite (numerical rank < number of "
                        "coefficients), OR",
                        "the separation sentinel trips (any |beta| > 10 or SE > 10)"],
                    "reported_always_never_gated_on": [
                        "per-stratum sample counts (available now, outcome-blind)",
                        "per-stratum event counts (at authorised execution)",
                        "per-stratum count of events with a NON-DEGENERATE "
                        "within-stratum risk set -- the events that actually "
                        "generate a comparison",
                        "the number of strata contributing at least one such event",
                        "convergence status, iteration count, final gradient norm",
                        "RISK-SET COMPARISON RETENTION N_strat / N_unstrat "
                        "(see REPLACEMENT_information_diagnostic) -- descriptive, "
                        "no threshold",
                        "the per-stratum breakdown of risk-set sizes at event times",
                        "the stratified and unstratified coefficient estimates with "
                        "their bootstrap intervals, side by side and unreduced to a "
                        "ratio",
                        "an explicit statement of the limitation beside the estimate"],
                    "REJECTED_the_squared_SE_ratio": {
                        "proposed_in_v2_1":
                            "median over coefficients of "
                            "(SE_unstratified / SE_stratified)^2",
                        "verdict": "REJECTED -- not defensible for THIS model",
                        "reasons": [
                            "ESTIMAND MISMATCH. In a stratified Cox, beta is a "
                            "WITHIN-stratum log hazard ratio; the unstratified beta "
                            "is a different functional of the data, equal only under "
                            "homogeneity of the effect AND proportionality of the "
                            "baseline hazards across strata. Comparing the "
                            "precisions of two different parameters is not an "
                            "information statement about either one.",
                            "PENALISED SEs ARE NOT NOMINAL. Under P6B-A1 the M2 "
                            "pathway coefficients are ridge-penalised, and a "
                            "penalised coefficient's model-based SE does not carry "
                            "the usual Wald interpretation. lambda is tuned "
                            "separately in each fit, so the ratio would partly "
                            "measure a difference in the TUNED PENALTY rather than "
                            "in information.",
                            "NO GUARANTEED DIRECTION. Stratification removes "
                            "between-stratum comparisons (lowering information) but "
                            "also removes confounding by stratum (which can lower "
                            "residual variation). The ratio is therefore not bounded "
                            "by 1, and a value above 1 is uninterpretable as "
                            "'retention'.",
                            "It is a ratio of ESTIMATES with its own sampling "
                            "variability, for which no uncertainty quantification "
                            "was specified."],
                        "not_forced_into_the_protocol": True},
                    "REPLACEMENT_information_diagnostic": {
                        "name": "risk-set comparison retention",
                        "formula":
                            "N_strat / N_unstrat, where N_unstrat = sum over events "
                            "i of (|R(t_i)| - 1) and N_strat = sum over events i of "
                            "(|R_{s(i)}(t_i)| - 1); R(t) is the risk set at t and "
                            "R_s(t) the risk set within the event's own stratum",
                        "why_it_is_defensible": [
                            "MODEL-FREE: a combinatorial property of the design and "
                            "the event times. It requires no fit, no SE, no lambda "
                            "and no convergence.",
                            "BOUNDED BY CONSTRUCTION: R_s(t) is a subset of R(t), so "
                            "the ratio lies in [0, 1] and the direction is "
                            "guaranteed.",
                            "DIRECTLY INTERPRETABLE: 'stratification retained X% of "
                            "the risk-set comparisons the partial likelihood could "
                            "have used'.",
                            "It subsumes the single-event-stratum question "
                            "automatically: an event whose within-stratum risk set "
                            "is itself alone contributes 0 comparisons, exactly "
                            "matching its log(1) = 0 likelihood contribution."],
                        "limitations_stated_not_hidden": [
                            "It counts COMPARISONS, not Fisher information. A "
                            "retained comparison between two subjects with identical "
                            "covariate vectors is counted but contributes nothing to "
                            "the information for beta.",
                            "It is reported alongside, not instead of, the number of "
                            "events with a non-degenerate within-stratum risk set.",
                            "It says nothing about precision of any particular "
                            "coefficient, and is NOT a substitute for the reported "
                            "confidence intervals.",
                            "It requires event times, so it is computed only at "
                            "authorised execution; its DEFINITION is outcome-blind "
                            "and is fixed here."],
                        "no_threshold": "DESCRIPTIVE ONLY. No pass/fail cut-point is "
                                        "placed on it, now or later; a cut-point "
                                        "would reintroduce the arbitrariness the "
                                        "withdrawn 50% rule was removed for."},
                    "diagnostics_are_necessary_not_sufficient": {
                        "statement":
                            "The three estimability criteria are NECESSARY "
                            "conditions for a usable fit. They are NOT sufficient "
                            "for adequate information or precision.",
                        "explicitly": [
                            "A positive-definite information matrix does NOT prove "
                            "adequate information or precision -- it proves only that "
                            "the parameters are locally identified and the Hessian is "
                            "invertible. A model can pass all three checks and still "
                            "have confidence intervals too wide to be useful.",
                            "A numerical failure is NOT the same thing as a wide "
                            "confidence interval, and the two are never reported as "
                            "if they were. A wide interval is a RESULT, reported as "
                            "such; a numerical failure means there is no result.",
                            "A zero-event stratum contributes no event information "
                            "and MUST NOT be treated as a numerical failure.",
                            "A one-event stratum is not automatically informative: "
                            "whether its event has a non-degenerate risk set is "
                            "reported, not assumed."]},
                    "response_if_declared_not_estimable":
                        "the stratified sensitivity model is REPORTED as not "
                        "estimable on this cohort and the leave-one-batch-out "
                        "diagnostic stands alone. We do NOT fall back to pooling: a "
                        "pooling choice made after seeing events-per-stratum would "
                        "be an outcome-informed design change.",
                    "why_events_per_stratum_is_the_wrong_quantity":
                        "a zero-event stratum contributes a likelihood factor of 1 "
                        "and drops out; an event whose within-stratum risk set is "
                        "itself alone contributes log(1) = 0 -- zero score, zero "
                        "information -- while still being counted in any naive "
                        "events-per-stratum tally",
                    "status": "APPROVED 2026-11-03 AND APPLIED"}},
            "interpretation": [
                "Agreement between the stratified and primary fits means the result is "
                "not driven by processing group. It does NOT mean batch effects are "
                "absent: ComBat was applied before deposit, so between-batch variance "
                "was removed by construction.",
                "Disagreement is reported as an unresolved limitation, NOT as grounds "
                "to adopt the adjusted model.",
                "The two surviving Phase 8B signals are a control-probe effect and a "
                "DISPERSION effect. Neither is a location effect, so neither would be "
                "addressed by a batch coefficient -- which is a further reason the "
                "adjustment is not promoted.",
            ],
            "why_never_primary": [
                "batch, RNA extraction method and MMR are mutually confounded in this "
                "deposit; no analysis of it can separate them",
                "a batch covariate would absorb biological variation, since batch05 is "
                "63.3% dMMR against 14.3% cohort-wide",
                "the deposited matrix is already ComBat-corrected; adjusting again "
                "risks removing signal twice",
                "the primary analysis was pre-specified without batch, and changing "
                "that after seeing results would forfeit its pre-specified status"],
            "status": PENDING,
        },
        {
            "id": "D9-1", "title": "Penalisation scope and clinical-predictor scaling",
            "decision_status": "APPROVED 2026-11-03 (Option B) -- amendment P6B-A1 APPLIED",
            "raised_by": "D8D-1 -- standardising the pathway scores exposes the "
                         "question for the clinical block, which it does not answer",
            "current_protocol_text": "Phase 6B §3.1 specifies 'ridge, lambda tuned "
                                     "inside resampling' for M1 without saying WHICH "
                                     "predictors the penalty applies to",
            "is_this_a_frozen_item": False,
            "frozen_status_note":
                "Penalisation scope is NOT a Register A item. A-11 fixes the "
                "predictor LIST; A-15 fixes model ROLES; A-17 the resampling. The "
                "penalty's scope is protocol text, so Phase 9 may finalise it -- but "
                "it changes the 6B draft, so it is raised for approval rather than "
                "settled silently.",
            "problem":
                "Ridge is not scale-invariant. The clinical block is four binaries "
                "plus a 2-column age spline basis on the YEAR scale. If they are "
                "penalised alongside the scores, the basis columns (values in the "
                "tens) carry small coefficients and are barely shrunk, while the "
                "binaries are shrunk comparatively hard -- an arbitrary ordering "
                "driven by units.",
            "options_considered": [
                "A. penalise everything, standardising the clinical block too "
                "(including binaries, which is unusual and makes a 0/1 contrast "
                "depend on the training fold's prevalence)",
                "B. force-include the clinical block UNPENALISED; penalise only the "
                "pathway scores in M2 and only the selected genes in M3",
                "C. leave the scope unspecified (status quo -- not viable, since the "
                "implementation must do something)"],
            "recommendation": "B",
            "recommendation_text":
                "The five clinical predictors (four binaries plus the age spline "
                "basis) are FORCE-INCLUDED and UNPENALISED in M1, M2 and M3. The "
                "ridge/elastic-net penalty applies ONLY to the pathway scores (M2) "
                "and to the candidate genes (M3).",
            "rationale": [
                "An unpenalised coefficient is scale-equivariant: the fit is "
                "invariant to any linear rescaling of that predictor. Option B "
                "therefore makes the clinical scaling question DISAPPEAR rather "
                "than answering it with a convention.",
                "M1 is the benchmark. With 6 df against 135-136 events (EPV ~22) a "
                "classical unpenalised Cox fit is appropriate and is the natural "
                "reference every genomic model is measured against.",
                "It makes delta-C(M2 - M1) attributable to the PATHWAY BLOCK rather "
                "than to differential shrinkage of the clinical terms, which is the "
                "comparison the protocol is built around.",
                "It is already what the M3 specification assumes ('clinical "
                "covariates force-included and unpenalised'), so option B makes the "
                "three models mutually consistent instead of M1/M2 differing from M3."],
            "amendment_status": "APPROVED 2026-11-03 AND APPLIED",
            "amendment_applied": True,
            "applied_amendment": {
                "amendment_id": "P6B-A1",
                "target": "docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md section 3.1",
                "field": "the 'Penalisation' row of the M1 table",
                "from": "ridge, lambda tuned inside resampling",
                "to": "the clinical block (four binary covariates and the age "
                      "spline basis) is FORCE-INCLUDED and UNPENALISED in M1, M2 "
                      "and M3. The ridge penalty applies only to the M2 pathway "
                      "scores and the elastic-net penalty only to the M3 candidate "
                      "genes; lambda (and alpha for M3) are tuned by inner "
                      "cross-validation inside each training fold. M1 therefore has "
                      "no penalised term and no lambda to tune.",
                "register_A_item": False,
                "consequence_to_state": "M1 ceases to be a penalised model; any "
                                        "comparison with the 6B draft text must "
                                        "say so.",
                "applied": True,
                "applied_to": ["docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md section 3.1 "
                               "+ new section 3.1.1",
                               "config/analysis/PHASE6B_AMENDMENT_RECORD.json"],
                "original_wording_preserved": True},
            "audit_findings_phase9_final_audit": {
                "nesting": "M1 is exactly M2 with the pathway coefficients set to "
                           "zero, so the models are properly nested. Had the "
                           "clinical block been penalised with a model-specific "
                           "lambda, M1 would NOT be the limit of M2 and the nesting "
                           "would be broken -- this choice is what MAKES the "
                           "comparison coherent.",
                "clinical_comparability": "IMPROVED -- unpenalised in both models, "
                                          "clinical coefficients differ only "
                                          "through genuine adjustment, not through "
                                          "a shrinkage artefact.",
                "molecular_estimand": "delta-C(M2 - M1) estimates the gain from "
                                      "adding a SHRUNK pathway block to an "
                                      "UNSHRUNK clinical model. The penalty pulls "
                                      "the pathway coefficients toward zero, i.e. "
                                      "toward M1, so the comparison is "
                                      "CONSERVATIVE for the molecular block. This "
                                      "must be stated whenever delta-C is reported.",
                "complementarity_with_D8D_1": "D8D-1 makes the penalty scale-fair "
                                              "across the pathway scores; D9-1 "
                                              "removes the clinical scaling "
                                              "question entirely. Neither alone is "
                                              "sufficient."},
            "limitations_of_the_recommendation": [
                "M1 becomes an unpenalised Cox model, so no lambda is tuned for it. "
                "That is a change from the 6B draft text and must be stated as such "
                "in any comparison with that draft.",
                "Unpenalised clinical terms are not shrunk, so M1's apparent "
                "performance carries the usual optimism -- which is exactly what the "
                "bootstrap correction estimates, and the uniform shrinkage factor is "
                "reported whether or not it is applied.",
                "If a future protocol wished to penalise the clinical block, the "
                "scaling question returns and option A's consequences would have to "
                "be faced then."],
            "status": PENDING,
        },
    ]
    APPLIED = {"D8D-1": "M2-A2", "D9-1": "P6B-A1", "D8D-4": "D8D-4-R1"}
    for d in DEC:
        if d["id"] in APPLIED:
            d["status"] = (f"APPROVED 2026-11-03 -- amendment "
                           f"{APPLIED[d['id']]} APPLIED")
        # D8D-2 / D8D-3 carry their approval in the entry itself: they formalise
        # a REPORTING protocol and require no specification amendment.
        print(f"  {d['id']}  {d['title'][:44]:46s} {d['status']}")
    n_pending = sum(1 for d in DEC if d["status"] == PENDING)
    assert n_pending == 0, [d["id"] for d in DEC if d["status"] == PENDING]
    assert all(d["status"].startswith("APPROVED") for d in DEC)

    # =====================================================================
    section("2b. M1 AGE FUNCTIONAL FORM -- finalised")
    # =====================================================================
    SPLINE = {
        "variable": "age.at.diagnosis_(year)",
        "form": "restricted cubic spline, 3 knots (2 df)",
        "why_3_knots": "with 6 df total against 135-136 events (EPV ~22), 3 knots "
                       "buys curvature without spending the budget; 4-5 knots would "
                       "be justified only at a far larger event count",
        "knot_rule": {
            "placement": "10th, 50th and 90th percentiles of age AMONG THE TRAINING "
                         "SAMPLES OF THAT FOLD",
            "quantile_definition": "numpy.percentile default (linear interpolation, "
                                   "type 7) -- named so two implementations agree",
            "recomputed_per_fold": True,
            "outcome_blind": "knots depend on the age distribution only; the outcome "
                             "is never inspected to place a knot or to choose the form"},
        "cohort_reference_values_outcome_blind": {
            "n_with_age": 453, "n_unique_ages": 89,
            "p10": 49.0, "median": 69.0, "p90": 83.0, "min": 22.0, "max": 97.0,
            "smallest_gap_between_knots_in_full_cohort_years": 14.0,
            "note": "quoted to show duplicate knots are implausible here, NOT to "
                    "place the knots -- placement is per training fold"},
        "reference_age_for_reporting": {
            "value": 69,
            "basis": "the cohort median age, an outcome-blind covariate summary",
            "fixed_now": True,
            "why_fixed_rather_than_per_fold":
                "the Cox model has no intercept, so the reference does not affect "
                "the FIT. It affects only the REPORTED hazard ratio HR(age = a vs "
                "a_ref). A per-fold reference would make reported HRs "
                "non-comparable across folds, so the reference is a single constant "
                "declared here and never changed.",
            "held_out_samples": "held-out ages are mapped through the TRAINING "
                                "fold's knots, unchanged; the reference is the same "
                                "fixed constant everywhere"},
        "degenerate_cases": {
            "duplicate_knots": "if two knots coincide (equal to within 1e-6), the "
                               "distinct knot values are used; if fewer than 3 "
                               "distinct knots remain, that fold FALLS BACK TO "
                               "LINEAR AGE and the fallback is recorded per fold",
            "rank_deficient_basis": "if the training-fold design matrix is "
                                    "rank-deficient (a basis column constant, or "
                                    "collinear to tolerance 1e-10), that fold falls "
                                    "back to linear age and records the fallback",
            "fallback_is_recorded_not_silent": True,
            "expected_frequency": "none expected: 89 distinct ages and 14- and "
                                  "20-year gaps between the full-cohort knots"},
        "prespecified_alternative": "a linear age term, compared by likelihood ratio "
                                    "WITHIN the training fold. The comparison is "
                                    "reported; it never selects the form using "
                                    "held-out data or outcomes",
        "all_transformations_inside_the_training_fold": True,
        "covariate_set_unchanged": ["age", "sex", "stage", "tumour location",
                                    "adjuvant chemotherapy"],
        "optional_covariates_not_added": "T, N, MMR and molecular covariates are NOT "
                                         "added; doing so needs approval (A-11)",
    }
    (OUT / "age_spline_specification.json").write_text(
        json.dumps(SPLINE, indent=2) + "\n")
    print(f"  form: {SPLINE['form']}")
    print(f"  knots: {SPLINE['knot_rule']['placement']}")
    print(f"  reference age (reporting only, fixed): {SPLINE['reference_age_for_reporting']['value']}")
    print(f"  degenerate-case fallback: linear age, recorded per fold")

    # =====================================================================
    section("3. LEAKAGE-CONTROL MATRIX")
    # =====================================================================
    LEAK = [
        ("frozen panel membership (which genes)", "OUTSIDE", "frozen metadata; fixed before any outcome contact", "A-13"),
        ("gene->probe representation map", "OUTSIDE", "frozen in Phase 2B; selection used dataset membership only, never outcome", "A-04/A-06"),
        ("cohort membership (S3)", "OUTSIDE", "frozen in Phase 1", "A-01"),
        ("QC exclusions", "OUTSIDE", "none performed; all 453 retained", "7C D2"),
        ("gene-level mu/sigma for pathway scoring", "INSIDE", "training fold only, sigma_ddof=1 from spec", "A-14 / M2-A1"),
        ("zero-variance gene exclusion", "INSIDE", "per fold, recorded per fold", "A-14"),
        ("pathway-score standardisation", "INSIDE (IF D8D-1 APPROVED)", "training fold only; unstandardised scores retained", "D8D-1 " + PENDING),
        ("expression/variance filter cut-point (M3)", "INSIDE", "value recomputed per training fold; the rule is frozen", "B-02"),
        ("univariate screening (M3)", "INSIDE", "training fold only; NO full-cohort screening at any point", "M3 spec"),
        ("feature selection (M3)", "INSIDE", "strictly nested within the training fold", "M3 spec"),
        ("penalty tuning (lambda, alpha)", "INSIDE", "inner CV within the training fold", "B-07"),
        # NOT "N/A": pandas reads that back as NaN and the row loses its value.
        ("imputation", "NOT_APPLICABLE", "complete-case; no imputation performed", "A-19"),
        ("PCA rotation, if ever used as a predictor", "INSIDE", "re-estimated per training fold; Phase 8A/8B rotations are descriptive only", "6B 3.2.2"),
        ("batch stratum definition (if D8D-4 approved)", "OUTSIDE", "pooling rule fixed from sample counts before unblinding", "D8D-4 " + PENDING),
        ("choice of B for the bootstrap", "OUTSIDE", "driven by MC-1..MC-4; never lowered to suit a result", "B-06"),
        ("random seeds", "OUTSIDE", "declared in advance; >=5 seeds for MC-2", "A-17"),
        ("outcome column join", "INSIDE the analysis only", "joined at model-fitting time; never during preprocessing or scoring", "A-16"),
    ]
    leak = pd.DataFrame(LEAK, columns=["operation", "permitted_location",
                                       "constraint", "authority"])
    leak.to_csv(OUT / "leakage_control_matrix.csv", index=False)
    for _, r in leak.iterrows():
        print(f"  {r.operation:44s} {r.permitted_location:28s} {r.authority}")
    assert not leak[leak.permitted_location.str.startswith("INSIDE")].empty

    # =====================================================================
    section("4. MODEL COMPARISON PLAN")
    # =====================================================================
    CMP = [
        ("M1", "clinical baseline (benchmark)", "age, sex, stage, location, adjuvant chemo",
         5, cc["M1 primary (clinical only)"]["complete_cases"], "CONFIRMATORY benchmark",
         "fixed first; never re-tuned afterwards"),
        ("M2 primary", "M1 + frozen 4-set scores", "M1 + EMT, E2F, WNT, P53",
         cc["M2 primary (4-set) + M1"]["n_predictors"],
         cc["M2 primary (4-set) + M1"]["complete_cases"], "PRIMARY pre-specified",
         "the reported primary result, whatever the sensitivity models show"),
        ("M2 sensitivity", "M1 + frozen 5-set scores", "M2 primary + PI3K",
         cc["M2 sensitivity (5-set) + M1"]["n_predictors"],
         cc["M2 sensitivity (5-set) + M1"]["complete_cases"], "SENSITIVITY only",
         "never replaces the primary; a metric gain is not grounds to switch"),
        ("M3", "transcriptome-wide", "M1 + nested-selected genes from 21,815",
         "data-driven inside folds", cc["M2 primary (4-set) + M1"]["complete_cases"],
         "EXPLORATORY only", "hypothesis-generating; no confirmatory claim"),
    ]
    cmp_df = pd.DataFrame(CMP, columns=["model", "description", "predictors",
                                        "n_parameters", "complete_case_n", "status",
                                        "role_constraint"])
    cmp_df.to_csv(OUT / "model_comparison_plan.csv", index=False)
    for _, r in cmp_df.iterrows():
        print(f"  {r.model:16s} p={str(r.n_parameters):24s} n={r.complete_case_n}  {r.status}")

    # =====================================================================
    section("5. VALIDATION / RESAMPLING SPECIFICATION")
    # =====================================================================
    VAL = {
        "internal_validation_design": {
            "primary": "bootstrap optimism correction (Harrell), the ENTIRE pipeline "
                       "repeated inside every replicate",
            "cross_check": "repeated 10-fold cross-validation",
            "external_validation": "NONE -- no suitable confirmatory cohort exists "
                                   "(Phase 4 G4 FAILED). No analysis here is external "
                                   "validation.",
            "authority": "6B §5"},
        "resampling_unit": {
            "unit": "the PATIENT (one array per patient in S3)",
            "justification": "Phase 8B established 0 declared dependent arrays in S3 "
                             "(dependancy_sample), so array and patient coincide here",
            "source": "reports/phase8b/phase8b_summary.json"},
        "stratification": {
            "bootstrap": "stratified by EVENT STATUS, so every replicate carries a "
                         "comparable number of events",
            "cross_validation": "folds stratified by event status",
            "outcome_blindness_note": "stratifying by event status REQUIRES the outcome "
                                      "and is therefore a MODEL-FITTING-TIME operation. "
                                      "It is specified here and executed only in Phase "
                                      "9 execution; no fold used in Phases 8A-8C was "
                                      "outcome-stratified, and those folds were "
                                      "explicitly labelled illustrative",
            "not_stratified_by": ["batch", "dataset partition", "MMR"],
            "why_not": "stratifying on a confounded covariate would make the replicates "
                       "unrepresentative of the cohort the estimand refers to"},
        "hyperparameter_tuning": {
            "boundary": "inner cross-validation INSIDE each training fold / bootstrap "
                        "replicate; never on held-out data, never on the full cohort",
            "parameters": ["ridge/elastic-net lambda", "elastic-net alpha (M3)"],
            "authority": "B-07"},
        "monte_carlo_stability": {
            "rules": "MC-1 running trace; MC-2 >=5 seeds; MC-3 MC-SE <= 10% of sampling "
                     "SE; MC-4 report MC-SE alongside every corrected estimate",
            "B_initial": 500,
            "B_is_an_output": True,
            "B_rule": "B increases (1000, 2000, ...) until MC-3 passes; it is NEVER "
                      "lowered to suit a result",
            "expected_hardest_case": "M3 nested selection",
            "authority": "6B §5.1"},
        "seed_policy": {
            "declared_seeds": [20261102, 20261103, 20261104, 20261105, 20261106],
            "n_seeds": 5,
            "rule": "seeds are declared HERE, before execution, and are never changed "
                    "after seeing a result. The MC-2 spread across these five seeds IS "
                    "the Monte-Carlo SE.",
            "note": "20261102 is the project seed already used by Phases 8A-8C"},
        "monte_carlo_vs_sampling_uncertainty": {
            "sampling_SE": "uncertainty from having 453 patients rather than infinitely "
                           "many; estimated by the bootstrap distribution itself; does "
                           "NOT shrink with more replicates",
            "monte_carlo_SE": "uncertainty from having run B replicates rather than "
                              "infinitely many; estimated as the SD of the corrected "
                              "estimate across the 5 declared seeds; DOES shrink as "
                              "B increases, roughly as 1/sqrt(B)",
            "reporting_rule": "both are reported for every corrected estimate, labelled "
                              "separately. They are never summed into one interval and "
                              "never conflated.",
            "criterion": "MC-SE <= 10% of sampling SE (MC-3)"},
        "optimism_correction": {
            "estimator": "Harrell's bootstrap optimism correction",
            "procedure": [
                "1. fit the ENTIRE pipeline on the original sample; record apparent "
                "performance P_app",
                "2. for b = 1..B: draw a bootstrap sample (stratified by event "
                "status), refit the ENTIRE pipeline on it -- including scoring "
                "parameters, any score standardisation, filtering, selection and "
                "tuning -- and record P_boot(b) on the bootstrap sample and "
                "P_orig(b) on the ORIGINAL sample using the bootstrap-fitted model",
                "3. optimism(b) = P_boot(b) - P_orig(b)",
                "4. corrected performance = P_app - mean_b optimism(b)"],
            "reported_together": "apparent performance, mean optimism, and corrected "
                                 "performance are reported TOGETHER. Apparent "
                                 "performance is never reported alone (6B §5).",
            "shrinkage": "a uniform shrinkage factor is estimated from the bootstrap "
                         "and reported whether or not it is applied"},
        "convergence_and_non_estimable_replicates": {
            "convergence_criteria": {
                "log_likelihood_change": "< 1e-9 between iterations",
                "max_iterations": 100,
                "separation_sentinel": "any |beta| > 10 or SE > 10 flags the fit"},
            "recording_rule":
                "every replicate that fails to converge, hits the iteration cap, "
                "trips the separation sentinel, or cannot be fitted (rank-deficient "
                "design, undefined scaling) is RECORDED with its seed, replicate "
                "index and failure reason in a per-model replicate log.",
            "exclusion_rule":
                "a non-estimable replicate is EXCLUDED from the optimism average and "
                "COUNTED. The count and percentage are reported beside every "
                "corrected estimate.",
            "pre_registered_threshold":
                "if more than 5% of replicates for a model are non-estimable, the "
                "corrected estimate for that model is reported as UNRELIABLE and the "
                "failure pattern is described. Replicates are never silently "
                "re-drawn to replace failures, which would bias the optimism "
                "estimate toward well-behaved samples.",
            "threshold_fixed_now": True},
        "tuning_nesting": {
            "structure": "inner cross-validation INSIDE each bootstrap replicate / "
                         "outer training fold",
            "inner_folds": 5,
            "what_is_tuned": ["ridge lambda (M2 pathway block)",
                              "elastic-net alpha and lambda (M3)"],
            "what_is_never_tuned_outside": "no hyperparameter is ever selected using "
                                           "held-out data or the full cohort",
            "score_standardisation_if_approved": "recomputed inside the INNER loop "
                                                 "too, not only the outer split"},
        "quantities_estimated": [
            "optimism-corrected Harrell C-index",
            "calibration slope on the LOG-HAZARD scale (unit-free)",
            "calibration-in-the-large on the log-hazard scale",
            "cause-specific hazard ratios with CIs",
            "delta-C vs M1 with CI and MC-SE",
            "uniform shrinkage factor (reported whether or not applied)",
            "effective degrees of freedom for penalised fits"],
        "quantities_NOT_estimated": [
            "t-year / fixed-horizon absolute risk (blocked by U1)",
            "predicted survival or CIF curves (blocked by U1)",
            "time-dependent AUC or Brier score (blocked by U1)",
            "cloglog calibration at a horizon, i.e. Riley's estimator (blocked by U1)",
            "any external-validation quantity (no cohort exists)"],
        "phase5_status": {
            "statement": "The Phase 5 simulation is an APPROXIMATION of the Riley et al. "
                         "(2022) external-validation sample-size framework, NOT a "
                         "faithful implementation of it. The ~338-event figure is an "
                         "approximation and was explicitly NOT a cross-check of Riley's "
                         "published ~365.",
            "source": "docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md §0.1"},
    }
    (OUT / "validation_specification.json").write_text(json.dumps(VAL, indent=2) + "\n")
    print(f"  design: {VAL['internal_validation_design']['primary'][:58]}...")
    print(f"  resampling unit: {VAL['resampling_unit']['unit']}")
    print(f"  stratification: {VAL['stratification']['bootstrap']}")
    print(f"  seeds: {VAL['seed_policy']['declared_seeds']}")
    print(f"  estimated: {len(VAL['quantities_estimated'])} quantities; "
          f"NOT estimated: {len(VAL['quantities_NOT_estimated'])}")

    # =====================================================================
    section("6. ESTIMAND SEPARATION (competing risks / composite)")
    # =====================================================================
    EST = [
        ("PRIMARY", "cause-specific recurrence hazard", EV,
         "cause-specific hazard ratio + discrimination",
         "censors competing death; answers 'rate among those still at risk'", "A-07/A-09"),
        ("S-1", "Fine-Gray subdistribution / CIF",
         next(c["value"] for c in CITED if c["label"] == "competing_events") + " competing",
         "absolute recurrence probability in the presence of death",
         "DIFFERENT estimand; coefficients are NOT comparable with the primary", "A-08"),
        ("S-2", "composite DFS (recurrence or death)",
         next(c["value"] for c in CITED if c["label"] == "composite_dfs_events"),
         "recurrence-or-death hazard",
         "DIFFERENT event definition; a larger event count is not more power for the "
         "primary question", "A-08"),
        ("S-3", "60-unit administrative horizon",
         next(c["value"] for c in CITED if c["label"] == "s3_horizon_events"),
         "same as primary, truncated",
         "CONDITIONAL on U1 -- '60 units' has no defined length while the unit is "
         "unresolved", "A-08/A-10"),
        ("S-7", "representative-IQR representation", EV, "same as primary",
         "aggregation-choice sensitivity; the vehicle for D8D-3", "A-08"),
    ]
    est = pd.DataFrame(EST, columns=["id", "estimand", "events", "quantity",
                                     "distinctness", "authority"])
    est.to_csv(OUT / "estimand_separation.csv", index=False)
    for _, r in est.iterrows():
        print(f"  {r.id:8s} {r.estimand:40s} events={str(r.events):14s} {r.authority}")

    # =====================================================================
    section("7. MODELING-READINESS CHECKLIST")
    # =====================================================================
    READY = [
        ("cohort frozen (S3, 453, no exclusions)", "READY", "A-01 / 7C D2"),
        ("endpoint + estimand set frozen", "READY", "A-07/A-08/A-09"),
        ("representation frozen (median primary)", "READY", "A-04/A-05"),
        ("M2 panels, roles, sigma_ddof=1 frozen", "READY", "A-13/A-14/M2-A1"),
        ("canonical scorer in src/", "READY", "8C rec. 4"),
        ("missing-data strategy frozen", "READY", "A-19 complete-case"),
        ("leakage rules specified", "READY", "6B §4 + this matrix"),
        ("resampling + MC stability rules frozen", "READY", "A-17 / 6B §5.1"),
        ("resampling unit + stratification specified", "READY", "Phase 9 §5"),
        ("seed policy declared", "READY", "Phase 9 §5"),
        ("performance measures restricted to U1-compatible set", "READY", "A-10 / 6A §3.1"),
        ("M1 parameter budget", "READY", "6 df, EPV ~22"),
        ("M3 nested-selection specification", "READY", "Phase 9 §4"),
        ("D8D-1 score standardisation", "APPROVED - APPLIED", "amendment M2-A2, spec revision 3"),
        ("D8D-2 4-set/5-set reporting framework", "APPROVED - FORMALISED", "4-set primary, 5-set sensitivity-only"),
        ("D8D-3 representation reporting", "APPROVED - FORMALISED", "median primary, rep-IQR sensitivity (S-7)"),
        ("D8D-4 batch sensitivity design", "APPROVED - APPLIED", "9 strata, no pooling; D8D-4-R1"),
        ("D9-1 penalisation scope / clinical scaling", "APPROVED - APPLIED", "amendment P6B-A1"),
        ("M1 age functional form", "READY", "3-knot RCS, training-fold knots, ref age 69"),
        ("convergence + non-estimable replicate policy", "READY", "5% threshold, fixed now"),
        ("optimism-correction procedure", "READY", "Harrell bootstrap, 4 steps"),
        ("complete-case EVENT count", "CONFIRM AT EXECUTION", f"in [{EV-len(dropped)}, {EV}]; not read here"),
        ("U1 time unit", "CONSTRAINS fixed-time analyses", "A-10"),
        ("U2 competing risks", "CONSTRAINS censoring assumption; S-1/S-2 registered", "A-08"),
        ("external confirmatory validation", "NOT AVAILABLE", "Phase 4 G4 FAILED"),
        ("protocol approval", "AWAITING USER", "this document"),
    ]
    chk = pd.DataFrame(READY, columns=["item", "status", "reference"])
    chk.to_csv(OUT / "phase9_readiness_checklist.csv", index=False)
    for _, r in chk.iterrows():
        print(f"  [{r['status']:34s}] {r['item']}")

    # =====================================================================
    section("8. CONFLICT SCAN against frozen decisions")
    # =====================================================================
    conflicts = []
    if current_std.startswith("NONE") and DEC[0]["recommendation"] == "B":
        conflicts.append({
            "type": "PROPOSED_AMENDMENT_NOT_A_CONFLICT",
            "detail": "D8D-1 recommends a change to a FROZEN field "
                      "(scoring.score_standardisation). It is presented as proposed "
                      "amendment M2-A2 and is NOT applied. The frozen records are "
                      "unchanged on disk.",
            "resolved_by": "explicit user approval"})
    p5 = json.loads((PANELS / "m2_panel_5set.json").read_text())
    if p4["role"] != "PRIMARY" or not p5["role"].startswith("SENSITIVITY"):
        conflicts.append({"type": "REAL_CONFLICT", "detail": "panel roles drifted"})
    for c in conflicts:
        print(f"  {c['type']}: {c['detail'][:92]}")
    assert not any(c["type"] == "REAL_CONFLICT" for c in conflicts)
    print(f"  real conflicts with frozen decisions: 0")

    # =====================================================================
    section("INTEGRITY")
    # =====================================================================
    after = {str(p.relative_to(paths.ROOT)): sha(p) for p in READONLY}
    changed = [k for k in before if before[k] != after[k]]
    assert not changed, f"READ-ONLY ARTIFACT MODIFIED: {changed}"
    print(f"  read-only artifacts unchanged : {len(READONLY)}/{len(READONLY)}")
    for n in ("m2_panel_4set.json", "m2_panel_5set.json"):
        ok, why = resolve_hash(PANELS / n, after[f"config/pathways/{n}"])
        print(f"  amendment chain {n:22s}: {why}")

    CHANGELOG = [
        {"version": "1.0", "date": "2026-11-03", "status": "superseded",
         "summary": "initial Phase 9 planning draft: recommendations for D8D-1..4, "
                    "leakage matrix, validation spec, estimand separation",
         "frozen_records_modified": "NONE"},
        {"version": "2.0", "date": "2026-11-03", "status": "superseded",
         "changes": [
             "D8D-1 EXPANDED: explicit zero-variance / non-finite scaling guards; "
             "an explicit answer that score standardisation does NOT settle "
             "clinical-predictor scaling; coefficient interpretation and the "
             "dispersion-vs-relevance distinction written out.",
             "D8D-4 REVERSED: the v1.0 rule pooling strata below 20 samples is "
             "WITHDRAWN. Stratification consumes no degrees of freedom, so the "
             "threshold was a reflex imported from the covariate case; pooling "
             "would have concealed the very heterogeneity being probed. All nine "
             "batches are kept as individual strata, with a pre-registered "
             "non-estimability response instead of a fallback to pooling.",
             "D8D-4: leave-one-batch-out caveat added -- batches differ ~7-fold in "
             "size, so the spread conflates batch identity with sample size.",
             "D9-1 ADDED: penalisation scope. D8D-1 exposed that the clinical block "
             "would inherit the same scale-dependence; recommendation is to leave "
             "the clinical block unpenalised, which removes the question rather "
             "than answering it by convention.",
             "AGE SPLINE FINALISED: knot rule, quantile definition, fixed reporting "
             "reference age (69), duplicate-knot and rank-deficiency fallbacks.",
             "VALIDATION COMPLETED: optimism-correction procedure written out in "
             "four steps; convergence criteria; non-estimable replicate recording "
             "with a 5% pre-registered threshold and a no-silent-redraw rule; inner "
             "tuning loop nesting made explicit.",
             "ARTIFACT RENAMED: modeling_readiness_checklist.csv -> "
             "phase9_readiness_checklist.csv (the name the brief specifies).",
         ],
         "frozen_records_modified": "NONE"},
        {"version": "2.1", "date": "2026-11-03", "status": "superseded",
         "issued_by": "Phase 9 Final Protocol Audit (docs/PHASE9_FINAL_AUDIT.md)",
         "changes": [
             "FINDING 1 -- batch counts: NO ERROR OF FACT. The nine counts were "
             "correct for the A-19 complete-case set (452); the audit request's "
             "premise conflated the S3 cohort (453) with that set. Join integrity "
             "verified against the Phase 3 bundle AND the ArrayExpress SDRF: 0 "
             "missing, 0 duplicated, 0 blank, 0 ambiguous. v2.1 now prints BOTH "
             "count columns side by side. All 453 remain in the cohort.",
             "FINDING 2 -- the D8D-4 '50% of events' rule is WITHDRAWN. It was "
             "logically near-vacuous (every event lies in a stratum containing an "
             "event, itself) and statistically unfounded (50% has no basis, and "
             "labelled a precision statement as estimability). Replaced by "
             "D8D-4-R1: estimability declared on convergence, information-matrix "
             "rank and the separation sentinel only, with information loss REPORTED "
             "as a continuous ratio and never thresholded.",
             "FINDING 3 -- D9-1 audited SOUND and retained. Nesting, clinical "
             "comparability and the conservative-deltaC statement recorded; "
             "amendment P6B-A1 drafted for Phase 6B 3.1 (NOT applied).",
             "FINDING 4 -- 33/33 internal consistency checks pass against frozen "
             "Register A.",
         ],
         "frozen_records_modified": "NONE"},
        {"version": "3.0", "date": "2026-11-03", "status": "superseded",
         "issued_by": "user authorisation to apply three pending amendments",
         "changes": [
             "M2-A2 APPLIED -- training-fold-only standardisation of the M2 "
             "pathway scores. Both panels and the freeze record move to "
             "specification revision 3; implementation added at "
             "src/score_standardisation.py. Panel MEMBERSHIP, ROLES and the "
             "gene-level sigma_ddof=1 are UNCHANGED, and src/pathway_scoring.py "
             "is unmodified.",
             "P6B-A1 APPLIED -- Phase 6B section 3.1 penalisation row amended; the "
             "original wording is preserved verbatim in the document and in the new "
             "config/analysis/PHASE6B_AMENDMENT_RECORD.json. Register A is "
             "untouched. The binding delta-C interpretation (conservative, NOT an "
             "unbiased estimate of the molecular contribution) is recorded.",
             "D8D-4-R1 APPLIED -- estimability declared on convergence, "
             "information-matrix rank and separation only. The v2.1 "
             "squared-SE information-retention ratio is REJECTED as not defensible "
             "for this model (estimand mismatch, non-nominal penalised SEs, no "
             "guaranteed direction) and REPLACED by risk-set comparison retention "
             "N_strat/N_unstrat -- model-free, bounded by construction, descriptive "
             "only, no threshold. The diagnostics are stated as NECESSARY but NOT "
             "SUFFICIENT for adequate precision.",
             "D8D-2 and D8D-3 remain PENDING USER APPROVAL -- they were not part of "
             "this authorisation.",
         ],
         "frozen_records_modified":
             "M2 panel specifications and M2 freeze record moved to revision 3 under "
             "the APPROVED amendment M2-A2, with every superseded hash retained in "
             "the amendment chain. A-19, QC records, Phase 2/3 artifacts, panel "
             "membership, cohort, endpoint and Phase 8 findings: NONE"},
    ]
    CHANGELOG.append({
        "version": "3.1", "date": "2026-11-03", "status": "current",
        "issued_by": "user authorisation of D8D-2 and D8D-3",
        "changes": [
            "D8D-2 APPROVED -- four-set Hallmark panel formalised as the "
            "PRE-SPECIFIED PRIMARY M2 analysis (EMT, E2F_TARGETS, "
            "WNT_BETA_CATENIN, P53_PATHWAY); five-set panel remains "
            "SENSITIVITY-ONLY, adding PI3K_AKT_MTOR_SIGNALING. The reporting "
            "protocol now distinguishes the pre-specified analysis HIERARCHY "
            "from OBSERVED PERFORMANCE, mandates reporting the P53-PI3K "
            "correlation (r = 0.6314) and the collinearity diagnostics with any "
            "five-set result, and records that PI3K is sensitivity-only because "
            "of CONDITIONING, explicitly NOT because it is biologically "
            "irrelevant.",
            "D8D-3 APPROVED -- median probe aggregation formalised as the PRIMARY "
            "gene representation; representative-IQR remains SENSITIVITY-ONLY "
            "(registered estimand S-7). WNT's weakest representation agreement "
            "(0.8577, the only score below 0.9) is acknowledged explicitly "
            "alongside E2F 0.9975, EMT 0.9920, P53 0.9720, PI3K 0.9096. "
            "Representation ROBUSTNESS and PREDICTIVE PERFORMANCE are recorded as "
            "different properties, never substituted for one another; agreement "
            "is explicitly NOT independent replication.",
            "No specification amendment was required: both decisions formalise "
            "REPORTING protocol. Panel membership, gene mappings, scoring "
            "definitions, the frozen representations and the M2-A1/M2-A2 "
            "amendment chain are UNCHANGED.",
            "All five Phase 9 decisions are now approved; none remains pending.",
        ],
        "frozen_records_modified": "NONE"})
    (OUT / "protocol_change_log.json").write_text(
        json.dumps({"document": "docs/PHASE9_MODELING_PROTOCOL_DRAFT.md",
                    "revisions": CHANGELOG}, indent=2) + "\n")
    (OUT / "decision_register.json").write_text(json.dumps(
        {"phase": "9", "stage": "PLANNING", "version": "3.1", "decisions": DEC,
         "decisions_applied": 0,
         "note": "Every recommendation is PENDING_USER_APPROVAL. No frozen record "
                 "was modified. The proposed M2-A2 amendment text is recorded but "
                 "NOT applied."}, indent=2) + "\n")
    summary = {
        "phase": "9", "stage": "PLANNING", "version": "3.1",
        "models_fitted": 0, "outcome_columns_opened": False,
        "outcome_files_blocked": list(BLOCKED_FILES),
        "event_counts_recomputed": False,
        "documented_figures_cited": CITED,
        "cohort": {"id": "S3", "n": S3_N, "exclusions": 0,
                   "complete_case_n_M1": cc["M1 primary (clinical only)"]["complete_cases"],
                   "complete_case_n_M2": cc["M2 primary (4-set) + M1"]["complete_cases"],
                   "a19_dropped": dropped,
                   "complete_case_events": f"[{EV-len(dropped)}, {EV}] -- not determined "
                                           f"here; requires outcome access"},
        "epv_M2_primary_range": [round(epv_lo, 2), round(epv_hi, 2)],
        "decisions": DEC,
        "amendments_applied": ["M2-A2", "P6B-A1", "D8D-4-R1"],
        "amendments_applied_count": 3,
        "decisions_still_pending": [d["id"] for d in DEC if d["status"] == PENDING],
        "decisions_approved": [d["id"] for d in DEC if d["status"].startswith("APPROVED")],
        "reporting_protocols_formalised": ["D8D-2", "D8D-3"],
        "all_recommendations_pending": False,
        "models_fitted_remains": 0,
        "leakage_matrix": leak.to_dict("records"),
        "model_comparison_plan": cmp_df.to_dict("records"),
        "validation_specification": VAL,
        "estimand_separation": est.to_dict("records"),
        "readiness": chk.to_dict("records"),
        "age_spline_specification": SPLINE,
        "batch_counts_both_sets": BATCH_COUNTS,
        "change_log": CHANGELOG,
        "conflicts": conflicts,
        "real_conflicts": 0,
        "frozen_state": {
            "panel_primary_role": p4["role"], "panel_sensitivity_role": p5["role"],
            "sigma_ddof": p4["scoring"]["sigma_ddof"],
            "score_standardisation_ON_DISK": current_std,
            "canonical_scorer_module": score_panel.__module__},
        "read_only_artifact_hashes": after,
        "read_only_artifacts_modified": False,
        "environment": environment_report(),
    }
    (OUT / "phase9_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(f"\n  outputs -> {OUT.relative_to(paths.ROOT)}/ "
          f"({len(list(OUT.glob('*.csv')))} CSV + {len(list(OUT.glob('*.json')))} JSON)")
    print(f"\nPHASE 9 PROTOCOL v3.1 -- 0 models fitted. Amendments APPLIED: "
          f"M2-A2, P6B-A1, D8D-4-R1. Reporting protocols FORMALISED: D8D-2, "
          f"D8D-3. Decisions still pending: "
          f"{[d['id'] for d in DEC if d['status'] == PENDING] or 'none'}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
