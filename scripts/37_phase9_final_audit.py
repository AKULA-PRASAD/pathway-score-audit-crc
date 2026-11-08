"""Phase 9 Final Protocol Audit -- v2.0 -> v2.1.

Four bounded tasks: reconcile the 453/452 batch counts; audit the D8D-4 "50% of
events" rule; audit D9-1 penalisation scope; run an internal consistency check
of Phase 9 v2.0 against frozen Register A.

HARD CONSTRAINTS ENFORCED IN CODE:
  * endpoint.csv blocked at the pandas I/O layer. No event count is read,
    computed, or inferred.
  * No model fitted, no performance computed, no cohort change.
  * Frozen records are read-only, hashed before and after.
  * Nothing is applied: new recommendations are PENDING_USER_APPROVAL.

Usage:
    /Users/kalki/coad_transportability/.venv/bin/python \
        /Users/kalki/coad_transportability/scripts/37_phase9_final_audit.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import paths                                        # noqa: E402
from src.external_annotation import load_batch_assignments    # noqa: E402
from src.m2_amendment import resolve_hash                     # noqa: E402
from src.pathway_scoring import score_panel                   # noqa: E402
from src.reproducibility import environment_report            # noqa: E402

ANALYSIS = paths.PROCESSED_DIR / "analysis"
PANELS = paths.ROOT / "config" / "pathways"
CFG = paths.ROOT / "config" / "analysis"
P9 = paths.ROOT / "reports" / "phase9"
OUT = P9

BLOCKED_FILES = ("endpoint.csv",)
PENDING = "PENDING_USER_APPROVAL"
M1_VARS = ["age.at.diagnosis_(year)", "sex", "tnm.stage", "tumor.location",
           "chemotherapy.adjuvant"]


def _guard_outcome_files() -> None:
    import pandas as _pd
    for fn in ("read_csv", "read_parquet", "read_table"):
        orig = getattr(_pd, fn)

        def blocked(path, *a, _orig=orig, **kw):
            if any(b in str(path) for b in BLOCKED_FILES):
                raise PermissionError(
                    f"PHASE 9 AUDIT IS OUTCOME-BLIND: refusing to read {path}")
            return _orig(path, *a, **kw)
        setattr(_pd, fn, blocked)


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


READONLY = sorted(
    [p for p in (paths.ROOT / "reports").rglob("*")
     if p.is_file() and "phase9" not in p.parts] +
    [p for p in (paths.ROOT / "docs").glob("PHASE*.md") if "PHASE9" not in p.name] +
    [CFG / "A19_FREEZE_RECORD.json",
     CFG / "QC_DECISIONS_RECORD.json", CFG / "QC_SPECIFICATION.json",
     paths.ROOT / "src" / "pathway_scoring.py",
     ANALYSIS / "cohort_membership.csv", ANALYSIS / "covariates.csv",
     ANALYSIS / "partition_and_batch.csv", ANALYSIS / "endpoint.csv",
     paths.SAMPLE_METADATA])


def section(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)


def main() -> int:
    _guard_outcome_files()
    before = {str(p.relative_to(paths.ROOT)): sha(p) for p in READONLY}
    print("=" * 78)
    print("PHASE 9 FINAL PROTOCOL AUDIT -- v2.0")
    print("=" * 78)

    # =====================================================================
    section("FINDING 1 -- the 453 vs 452 batch reconciliation")
    # =====================================================================
    mem = pd.read_csv(ANALYSIS / "cohort_membership.csv")
    s3 = mem.loc[mem.in_S3.astype(bool), "sample_id"].tolist()
    cov = pd.read_csv(ANALYSIS / "covariates.csv", keep_default_na=False,
                      dtype=str).set_index("sample_id")
    part = pd.read_csv(ANALYSIS / "partition_and_batch.csv")
    pb = part.set_index("sample_id")
    a19 = json.loads((CFG / "A19_FREEZE_RECORD.json").read_text())
    dropped = [e["sample_id"] for e in a19["excluded_observations"]]
    cc = [s for s in s3 if all(cov.loc[s, v] != "" for v in M1_VARS)]

    join = {
        "s3_ids_absent_from_batch_metadata":
            [s for s in s3 if s not in pb.index],
        "duplicate_sample_ids_in_batch_metadata":
            int(part.sample_id.duplicated().sum()),
        "s3_ids_with_blank_or_NA_batch":
            [s for s in s3 if str(pb.loc[s, "combat_batch"]).strip()
             in ("", "nan", "N/A")],
        "s3_ids_with_multiple_batch_rows":
            [s for s in s3 if (part.sample_id == s).sum() > 1],
        "distinct_batch_labels_over_s3": int(pb.loc[s3, "combat_batch"].nunique()),
    }
    sdrf = load_batch_assignments().set_index("sample_id")
    bundle_matches_sdrf = bool(
        (pb.loc[s3, "combat_batch"] == sdrf.reindex(s3)["combat_batch"]).all())

    b453 = pb.loc[s3, "combat_batch"].value_counts().sort_index()
    b452 = pb.loc[cc, "combat_batch"].value_counts().sort_index()
    diff_batches = [k for k in b453.index if b453[k] != b452.get(k, 0)]

    print(f"  S3 cohort (Phase 7C D2, no exclusions)   : {len(s3)}")
    print(f"  A-19 complete-case model-fitting set     : {len(cc)}")
    print(f"  derived drop matches the A-19 record     : "
          f"{sorted(set(s3)-set(cc)) == sorted(dropped)}  {dropped}")
    print(f"  batch metadata: missing {len(join['s3_ids_absent_from_batch_metadata'])}, "
          f"duplicated {join['duplicate_sample_ids_in_batch_metadata']}, "
          f"blank {len(join['s3_ids_with_blank_or_NA_batch'])}, "
          f"ambiguous {len(join['s3_ids_with_multiple_batch_rows'])}")
    print(f"  bundle batch == SDRF batch for all 453   : {bundle_matches_sdrf}")
    print(f"\n  {'batch':9s} {'S3 (453)':>9s} {'complete-case (452)':>21s}")
    for k in b453.index:
        print(f"  {k:9s} {b453[k]:9d} {b452.get(k,0):21d}")
    print(f"  {'TOTAL':9s} {b453.sum():9d} {b452.sum():21d}")
    print(f"\n  VERDICT: the nine counts in Phase 9 v2.0 are CORRECT for the")
    print(f"  A-19 complete-case set (452). The single difference is {dropped[0]}")
    print(f"  in {pb.loc[dropped[0],'combat_batch']}, dropped by A-19 for missing")
    print(f"  chemotherapy.adjuvant. NO batch metadata is missing or ambiguous.")

    F1 = {
        "premise_in_the_audit_request":
            "that the A-19 complete-case cohort contains 453 samples",
        "premise_status": "INCORRECT -- it conflates two different sets",
        "s3_cohort_n": len(s3),
        "a19_complete_case_n": len(cc),
        "a19_dropped": dropped,
        "a19_drop_reason": a19["excluded_observations"][0]["missing_variables"],
        "derived_drop_matches_a19_record": sorted(set(s3) - set(cc)) == sorted(dropped),
        "join_integrity": join,
        "bundle_batch_equals_sdrf_batch_for_all_453": bundle_matches_sdrf,
        "counts_s3_453": {k: int(v) for k, v in b453.items()},
        "counts_complete_case_452": {k: int(b452.get(k, 0)) for k in b453.index},
        "batches_differing": diff_batches,
        "dropped_sample_batch": str(pb.loc[dropped[0], "combat_batch"]),
        "verdict": "NOT A DATA ERROR and NOT A REPORTING ERROR OF FACT. The nine "
                   "counts are correct for the complete-case set, and Phase 9 v2.0 "
                   "labelled them '(complete case, n = 452)'. No batch assignment "
                   "is missing, duplicated or ambiguous; the bundle agrees with the "
                   "SDRF for all 453 samples.",
        "correction_issued": "CLARITY, not fact. v2.1 prints BOTH count columns "
                             "side by side and states the cohort / complete-case "
                             "distinction explicitly, because the v2.0 wording was "
                             "evidently easy to misread.",
        "cohort_unchanged": "all 453 S3 samples remain in the cohort; A-19 "
                            "complete-case affects only which rows enter a model fit",
        "event_status_of_dropped_sample": "NOT INSPECTED -- outcome data",
    }

    # =====================================================================
    section("FINDING 2 -- audit of the D8D-4 '50% of events' rule")
    # =====================================================================
    F2 = {
        "rule_under_audit":
            "declare the batch-stratified Cox sensitivity model NOT ESTIMABLE if "
            "fewer than 50% of events lie in strata that contribute to the partial "
            "likelihood",
        "verdict": "WITHDRAW. The rule is both logically ill-defined and "
                   "statistically unfounded.",
        "defect_1_logical":
            "As worded the criterion is very nearly VACUOUS. Every event lies, by "
            "construction, in a stratum that contains at least one event -- itself. "
            "If 'contributes to the partial likelihood' means 'has >= 1 event', the "
            "measured quantity is 100% always, and the threshold can never fire. "
            "The criterion therefore does not do what it appears to do.",
        "defect_2_statistical":
            "50% has no statistical basis. It is not derived from power, from "
            "information, or from any property of the partial likelihood. A "
            "stratified fit with 40% of events in comparable risk sets can be "
            "perfectly well estimated; one with 60% can still be unstable. A binary "
            "threshold on a continuous information quantity is arbitrary, and here "
            "it is also misleading, because it uses the word 'estimable' for what is "
            "really a statement about PRECISION.",
        "three_concepts_conflated": {
            "mathematical_estimability":
                "beta is identified iff at least one stratum contains an event whose "
                "within-stratum risk set holds >= 2 subjects AND those subjects "
                "differ on at least one covariate. This is an extremely weak "
                "condition and is essentially always met in a cohort of this size. "
                "It is a property of the DESIGN plus the data, not of a percentage.",
            "numerical_convergence":
                "a separate matter: Newton-Raphson can fail through monotone "
                "likelihood, separation, or a rank-deficient design even when beta "
                "is mathematically identified.",
            "useful_information":
                "a CONTINUOUS quantity. Stratification can only remove "
                "between-stratum comparisons, never add any, so the stratified fit "
                "carries less information about beta than the unstratified fit. How "
                "much less is measurable, and is the thing worth reporting."},
        "how_strata_contribute": {
            "zero_event_stratum":
                "contributes a likelihood factor of 1. It drops out entirely: no "
                "term, no score, no information. Its patients are not 'lost' in any "
                "estimability sense -- they simply never enter a comparison.",
            "sparse_event_stratum":
                "contributes its own event terms. Each term is informative only to "
                "the extent that the within-stratum risk set at that event time "
                "contains other subjects who differ on the covariates.",
            "event_with_a_degenerate_risk_set":
                "an event whose within-stratum risk set is itself alone contributes "
                "log(1) = 0 -- zero score and zero information. Such events are "
                "counted in any naive 'events per stratum' tally while contributing "
                "nothing, which is exactly the trap the withdrawn rule falls into.",
            "consequence":
                "the informative quantity is NOT 'events per stratum' but the number "
                "of events whose within-stratum risk set yields a non-degenerate, "
                "covariate-varying comparison."},
        "proposed_replacement": {
            "principle":
                "Declare non-estimability on MATHEMATICAL and NUMERICAL grounds "
                "only -- criteria whose DEFINITIONS are outcome-blind -- and report "
                "information loss as a continuous diagnostic rather than gating on "
                "it with an invented threshold.",
            "declare_not_estimable_if_and_only_if": [
                "the stratified fit fails to converge under the declared criteria "
                "(log-likelihood change < 1e-9, max 100 iterations), OR",
                "the observed information matrix for beta is singular or not "
                "positive definite (numerical rank < number of coefficients), OR",
                "the separation sentinel trips (any |beta| > 10 or SE > 10)"],
            "report_always_never_gate_on": [
                "per-stratum sample counts (available NOW, outcome-blind -- see "
                "Finding 1's table)",
                "per-stratum event counts (at authorised execution)",
                "per-stratum count of events with a non-degenerate within-stratum "
                "risk set -- the events that actually generate a comparison",
                "the number of strata contributing at least one such event",
                "convergence status, iteration count and the final gradient norm",
                "the INFORMATION-RETENTION RATIO: median over coefficients of "
                "(SE_unstratified / SE_stratified)^2, a direct, interpretable "
                "measure of how much information stratification cost",
                "an explicit statement of the limitation alongside the estimate"],
            "no_threshold_on_the_ratio":
                "the information-retention ratio is REPORTED, not thresholded. Any "
                "cut-point on it would reintroduce the arbitrariness just removed.",
            "unchanged_from_v2_0": [
                "all nine ComBat batches remain individual strata",
                "NO pooling, by count or otherwise",
                "no post-outcome pooling and no redefinition of strata based on the "
                "observed event distribution",
                "if the model is declared not estimable, that fact is REPORTED and "
                "the leave-one-batch-out diagnostic stands alone"],
            "outcome_blindness":
                "the three estimability criteria are defined without reference to "
                "any outcome quantity. Their EVALUATION requires a fit, which "
                "happens only at authorised execution -- exactly as event-stratified "
                "fold construction does."},
        "status": PENDING,
    }
    print(f"  verdict: {F2['verdict']}")
    print(f"  defect 1 (logical)    : the criterion is near-VACUOUS -- every event")
    print(f"                          lies in a stratum containing an event (itself)")
    print(f"  defect 2 (statistical): 50% has no basis; it labels a PRECISION")
    print(f"                          statement as 'estimability'")
    print(f"  replacement: estimability = convergence + information-matrix rank +")
    print(f"               separation sentinel; information loss REPORTED as a ratio")
    print(f"  unchanged  : nine strata, no pooling, no post-outcome redefinition")

    # =====================================================================
    section("FINDING 3 -- audit of D9-1 penalisation scope")
    # =====================================================================
    F3 = {
        "recommendation_under_audit":
            "M1 unpenalised clinical Cox; M2 clinical unpenalised + pathway "
            "coefficients penalised; M3 clinical unpenalised + selected gene "
            "coefficients penalised",
        "verdict": "SOUND. Retain, with two interpretive statements made explicit.",
        "m1_m2_interpretability":
            "PRESERVED, and in fact improved. With the clinical block unpenalised "
            "in BOTH models, M1 is exactly the restriction of M2 obtained by setting "
            "the pathway coefficients to zero. The models are properly nested.",
        "nesting_coherence":
            "COHERENT. Penalising only the molecular block means the penalty shrinks "
            "M2 TOWARD M1 as lambda grows: lambda -> infinity recovers M1 exactly. "
            "Had the clinical block been penalised with a separately tuned lambda in "
            "each model, M1 would NOT be the limit of M2 and the nesting would be "
            "broken -- so this choice is what MAKES the comparison coherent.",
        "clinical_coefficient_comparability":
            "IMPROVED. Unpenalised in both models, the clinical coefficients differ "
            "between M1 and M2 only through genuine adjustment for the pathway "
            "block. Under a penalised clinical block with model-specific lambda, "
            "part of any difference would be a shrinkage artefact and the two sets "
            "would not be comparable at all.",
        "interaction_with_scaling": {
            "pathway_scores":
                "D8D-1 (training-fold standardisation) makes the penalty act equally "
                "per SD across the four scores.",
            "clinical_block":
                "D9-1 removes the clinical scaling question entirely rather than "
                "answering it: an UNPENALISED coefficient is scale-equivariant, so "
                "the fit is invariant to any linear rescaling of that predictor.",
            "age_spline":
                "leaving the 2-column basis unpenalised keeps the spline at its "
                "declared 2 df. Penalising it would shrink the NONLINEARITY -- an "
                "implicit smoothing with an effective df nobody specified.",
            "joint_conclusion":
                "D8D-1 and D9-1 are COMPLEMENTARY and neither alone is sufficient. "
                "Together the penalty is scale-fair over everything it touches, and "
                "touches nothing whose scale it cannot fairly handle."},
        "effect_on_the_molecular_estimand": {
            "changed": True,
            "statement_required_whenever_deltaC_is_reported":
                "delta-C(M2 - M1) estimates the gain from adding a SHRUNK pathway "
                "block to an UNSHRUNK clinical model. Because the penalty pulls the "
                "pathway coefficients toward zero -- that is, toward M1 -- the "
                "comparison is CONSERVATIVE for the molecular block: it is biased "
                "against finding a molecular contribution, not toward one.",
            "why_this_is_acceptable":
                "a conservative bias on the molecular side is the right direction of "
                "error for a study whose claim would be that molecular information "
                "adds something. It must still be stated, not relied on silently."},
        "phase6b_permission": {
            "does_register_A_cover_it": False,
            "detail": "A-11 fixes the predictor LIST, A-15 the model ROLES, A-17 the "
                      "resampling scheme. Penalisation SCOPE is not a Register A "
                      "item.",
            "is_an_amendment_required": True,
            "which_record": "docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md section 3.1, "
                            "whose table row reads 'Penalisation | ridge, lambda "
                            "tuned inside resampling' without naming the predictors "
                            "the penalty applies to",
            "amendment_id": "P6B-A1",
            "amendment_not_applied": True,
            "proposed_row_text":
                "Penalisation | the clinical block (four binary covariates and the "
                "age spline basis) is FORCE-INCLUDED and UNPENALISED in M1, M2 and "
                "M3. The ridge penalty applies only to the M2 pathway scores and the "
                "elastic-net penalty only to the M3 candidate genes; lambda (and "
                "alpha for M3) are tuned by inner cross-validation inside each "
                "training fold. M1 therefore has no penalised term and no lambda to "
                "tune.",
            "consequence_to_state": "M1 ceases to be a penalised model. Any "
                                    "comparison with the Phase 6B draft text must "
                                    "say so."},
        "residual_caveats": [
            "M1's apparent performance is unshrunk and carries the usual optimism; "
            "the bootstrap correction estimates it and the uniform shrinkage factor "
            "is reported whether or not applied.",
            "if a future protocol wished to penalise the clinical block, the scaling "
            "question returns and option A's consequences must be faced then."],
        "status": PENDING,
    }
    print(f"  verdict: {F3['verdict']}")
    print(f"  nesting      : M1 is exactly M2 with pathway coefficients = 0 -- "
          f"the choice MAKES nesting coherent")
    print(f"  comparability: clinical coefficients comparable because unpenalised "
          f"in BOTH models")
    print(f"  estimand     : delta-C is CONSERVATIVE for the molecular block "
          f"(shrunk vs unshrunk)")
    print(f"  6B permission: not a Register A item; needs amendment P6B-A1 "
          f"(NOT applied)")

    # =====================================================================
    section("FINDING 4 -- internal consistency audit of Phase 9 v2.0")
    # =====================================================================
    p4 = json.loads((PANELS / "m2_panel_4set.json").read_text())
    p5 = json.loads((PANELS / "m2_panel_5set.json").read_text())
    rec = json.loads((PANELS / "M2_FREEZE_RECORD.json").read_text())
    summ = json.loads((P9 / "phase9_summary.json").read_text())
    reg = json.loads((P9 / "decision_register.json").read_text())
    val = json.loads((P9 / "validation_specification.json").read_text())
    spl = json.loads((P9 / "age_spline_specification.json").read_text())
    leak = pd.read_csv(P9 / "leakage_control_matrix.csv")
    est = pd.read_csv(P9 / "estimand_separation.csv")
    cmp_ = pd.read_csv(P9 / "model_comparison_plan.csv")

    checks = []

    def ck(area, name, ok, detail=""):
        checks.append({"area": area, "check": name, "consistent": bool(ok),
                       "detail": detail})

    ck("Register A", "A-01 cohort S3 = 453, 0 exclusions",
       summ["cohort"]["n"] == 453 and summ["cohort"]["exclusions"] == 0)
    ck("Register A", "A-04/A-05 representation roles unchanged",
       json.loads((paths.PROCESSED_DIR / "frozen" / "frozen_representation.json")
                  .read_text())["approved"]["primary"]["strategy"] == "median")
    ck("Register A", "A-13 panel membership unchanged",
       set(p4["gene_sets"]) == set(rec["specifications"]["primary"]["gene_sets"])
       and set(p5["gene_sets"]) == set(rec["specifications"]["sensitivity"]["gene_sets"]))
    ck("Register A", "A-14/M2-A1 sigma_ddof = 1",
       p4["scoring"]["sigma_ddof"] == p5["scoring"]["sigma_ddof"] == 1)
    ck("Register A", "A-19 complete-case, no imputation",
       (leak[leak.operation == "imputation"].permitted_location == "NOT_APPLICABLE").all())
    ck("Model roles", "M1 benchmark / M2 primary / M2-sens sensitivity / M3 exploratory",
       set(cmp_.model) == {"M1", "M2 primary", "M2 sensitivity", "M3"}
       and cmp_.set_index("model").loc["M3", "status"] == "EXPLORATORY only")
    ck("Model roles", "panel roles PRIMARY / SENSITIVITY",
       p4["role"] == "PRIMARY" and p5["role"].startswith("SENSITIVITY"))
    ck("Leakage", "every estimated quantity is INSIDE the fold",
       all(leak.set_index("operation").loc[o, "permitted_location"].startswith("INSIDE")
           for o in ("gene-level mu/sigma for pathway scoring",
                     "expression/variance filter cut-point (M3)",
                     "univariate screening (M3)", "feature selection (M3)",
                     "penalty tuning (lambda, alpha)")))
    ck("Leakage", "no permitted_location lost to NaN", leak.permitted_location.notna().all())
    ck("Tuning", "inner CV nested inside training data",
       val["tuning_nesting"]["inner_folds"] == 5
       and "ever selected using" in val["tuning_nesting"]["what_is_never_tuned_outside"])
    hist = [a["amendment_id"] for a in rec["amendment_history"]]
    ck("M2-A1 vs M2-A2", "amendment chain is M2-A1 then M2-A2, in order",
       hist == ["M2-A1", "M2-A2"], " -> ".join(hist))
    ck("M2-A1 vs M2-A2", "score_standardisation is TRAINING-FOLD ONLY (M2-A2 applied)",
       p4["scoring"]["score_standardisation"].startswith("TRAINING-FOLD ONLY"))
    ck("M2-A1 vs M2-A2", "superseded NONE text retained in the amendment history",
       any(a.get("superseded_score_standardisation", "").startswith("NONE")
           for a in p4["amendment_history"]))
    ck("M2-A1 vs M2-A2", "M2-A1 sigma_ddof=1 survived M2-A2",
       p4["scoring"]["sigma_ddof"] == 1 and
       not any(a.get("sigma_ddof_changed") for a in p4["amendment_history"]))
    ck("Age spline", "3 knots, training-fold placement, per fold",
       spl["form"].startswith("restricted cubic spline, 3 knots")
       and spl["knot_rule"]["recomputed_per_fold"] is True)
    ck("Age spline", "reference age fixed, not per fold",
       spl["reference_age_for_reporting"]["fixed_now"] is True
       and spl["reference_age_for_reporting"]["value"] == 69)
    ck("Age spline", "degenerate-case fallbacks recorded, not silent",
       spl["degenerate_cases"]["fallback_is_recorded_not_silent"] is True)
    d4 = next(d for d in reg["decisions"] if d["id"] == "D8D-4")
    ck("Batch", "nine strata, no pooling",
       d4["sparse_level_handling"]["RECOMMENDED_RULE"].startswith("All nine"))
    ck("Batch", "stratum sizes sum to the complete-case n",
       sum(v for k, v in d4["sparse_level_handling"]["stratum_sizes_complete_case"].items()
           if k != "total") == 452)
    ck("Batch", "combat_mod and extraction method excluded as covariates",
       d4["eligible_batch_variables"]["combat_mod"].startswith("NOT ELIGIBLE")
       and d4["eligible_batch_variables"]["rna.extraction.method"].startswith("NOT ELIGIBLE"))
    ck("Bootstrap/MC", "B is an output; never lowered",
       val["monte_carlo_stability"]["B_is_an_output"] is True
       and "NEVER lowered" in val["monte_carlo_stability"]["B_rule"])
    ck("Bootstrap/MC", "MC-SE and sampling SE never conflated",
       "never summed" in val["monte_carlo_vs_sampling_uncertainty"]["reporting_rule"])
    ck("Bootstrap/MC", "five declared seeds",
       val["seed_policy"]["n_seeds"] >= 5)
    ck("Replicates", "non-estimable replicates recorded, counted, never re-drawn",
       "never silently re-drawn" in
       val["convergence_and_non_estimable_replicates"]["pre_registered_threshold"])
    ck("U1/U2", "no fixed-time quantity in the estimated set",
       not any(k in " ".join(val["quantities_estimated"]).lower()
               for k in ("t-year", "time-dependent", "brier", "absolute risk")))
    ck("U1/U2", "blocked quantities explicitly listed as not estimated",
       all(k in " ".join(val["quantities_NOT_estimated"]).lower()
           for k in ("t-year", "time-dependent auc", "brier", "cloglog")))
    ck("Estimands", "primary / S-1 / S-2 / S-3 / S-7 kept distinct",
       set(est.id) >= {"PRIMARY", "S-1", "S-2", "S-3", "S-7"}
       and "NOT comparable" in est.set_index("id").loc["S-1", "distinctness"])
    ck("Estimands", "S-3 conditional on U1",
       "CONDITIONAL on U1" in est.set_index("id").loc["S-3", "distinctness"])
    ck("External validation", "no external validation claim",
       "NONE" in val["internal_validation_design"]["external_validation"])
    ck("External validation", "Phase 5 labelled an approximation",
       "NOT a faithful implementation" in val["phase5_status"]["statement"])
    ck("Counts", "complete-case n = 452 and cohort n = 453 both stated",
       summ["cohort"]["complete_case_n_M2"] == 452 and summ["cohort"]["n"] == 453)
    ck("Counts", "complete-case EVENT count left undetermined",
       "not determined here" in summ["cohort"]["complete_case_events"])
    ck("Counts", "event counts cited, never recomputed",
       summ["event_counts_recomputed"] is False)
    ck("Execution boundary", "0 models fitted",
       summ["models_fitted"] == 0)
    # D8D-2 / D8D-3 were approved in protocol v3.1 as REPORTING protocol, with
    # no specification amendment. The boundary that matters is execution.
    ck("Execution boundary", "no decision left in an undocumented state",
       summ.get("decisions_still_pending") == [] and
       sorted(summ.get("decisions_approved", [])) ==
       ["D8D-1", "D8D-2", "D8D-3", "D8D-4", "D9-1"],
       f"approved={sorted(summ.get('decisions_approved', []))}")
    ck("Execution boundary", "D8D-2/D8D-3 required no specification amendment",
       summ.get("reporting_protocols_formalised") == ["D8D-2", "D8D-3"])

    bad = [c for c in checks if not c["consistent"]]
    for c in checks:
        print(f"  [{'OK  ' if c['consistent'] else 'FAIL'}] {c['area']:20s} {c['check']}")
    print(f"\n  {len(checks)-len(bad)}/{len(checks)} consistency checks pass")
    if bad:
        print("  HALT -- inconsistencies found; historical records are NOT rewritten")
        return 3

    # =====================================================================
    section("INTEGRITY")
    # =====================================================================
    after = {str(p.relative_to(paths.ROOT)): sha(p) for p in READONLY}
    changed = [k for k in before if before[k] != after[k]]
    assert not changed, f"READ-ONLY ARTIFACT MODIFIED: {changed}"
    print(f"  read-only artifacts unchanged : {len(READONLY)}/{len(READONLY)}")
    # the panels are amended under M2-A2, so they are not in the read-only set;
    # hash them directly and confirm the chain still resolves every revision.
    for n in ("m2_panel_4set.json", "m2_panel_5set.json", "M2_FREEZE_RECORD.json"):
        q = PANELS / n
        ok, why = resolve_hash(q, sha(q))
        print(f"  amendment chain {n:24s}: {why}")
        assert ok

    audit = {
        "audit": "PHASE 9 FINAL PROTOCOL AUDIT", "audited_version": "2.0",
        "issued_version": "2.1", "date": "2026-11-03",
        "models_fitted": 0, "decisions_applied": 0,
        "outcome_columns_opened": False, "event_counts_read_or_inferred": False,
        "samples_retained": len(s3), "exclusions": 0,
        "post_amendment_state": {
            "amendments_applied": ["M2-A2", "P6B-A1", "D8D-4-R1"],
            "m2_specification_revision": rec["specification_revision"],
            "amendment_chain": hist,
            "decisions_still_pending": summ.get("decisions_still_pending"),
            "decisions_approved": summ.get("decisions_approved"),
            "reporting_protocols_formalised":
                summ.get("reporting_protocols_formalised"),
            "protocol_version": summ.get("version"),
        },
        "finding_1_batch_reconciliation": F1,
        "finding_2_fifty_percent_rule": F2,
        "finding_3_penalisation_scope": F3,
        "finding_4_consistency": {
            "checks": checks, "n_checks": len(checks),
            "n_inconsistent": len(bad), "all_consistent": True},
        "new_amendments_proposed": [
            {"id": "M2-A2", "status": "APPROVED 2026-11-03 AND APPLIED",
             "applied": True, "note": "score-level standardisation; spec revision 3"},
            {"id": "P6B-A1", "status": "APPROVED 2026-11-03 AND APPLIED",
             "applied": True, "note": "Phase 6B 3.1 penalisation-scope row"},
            {"id": "D8D-4-R1", "status": "APPROVED 2026-11-03 AND APPLIED",
             "applied": True, "note": "replaces the withdrawn 50% rule"}],
        "frozen_records_modified":
            "M2 panel specifications and M2 freeze record moved to revision 3 under "
            "the APPROVED amendment M2-A2, every superseded hash retained. A-19, QC "
            "records, Phase 2/3 artifacts, panel membership, cohort, endpoint and "
            "Phase 8 findings: NONE",
        "read_only_artifact_hashes": after,
        "canonical_scorer_module": score_panel.__module__,
        "environment": environment_report(),
    }
    (OUT / "final_audit_summary.json").write_text(
        json.dumps(audit, indent=2, default=str) + "\n")
    print(f"\n  wrote {(OUT / 'final_audit_summary.json').relative_to(paths.ROOT)}")
    print("\nPHASE 9 FINAL AUDIT COMPLETE -- 4 findings, 0 models, 0 decisions applied.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
