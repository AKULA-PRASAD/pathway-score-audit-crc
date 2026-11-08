# Phase 9 — Alignment of θ̂ with the corrected nested-bootstrap interval

**Date:** 7 November 2026. **Scope:** read-only audit. No estimate was recomputed or replaced, and no frozen file, configuration, checkpoint or historical result was changed.
**Question:** does the pre-specified point estimate θ̂ come from the same estimator whose sampling distribution the corrected nested bootstrap approximates, so that the two can be reported together?

**Conclusion: yes, at the level of the estimator.** θ̂ and the corrected outer replicates are the same functional, implemented by the same frozen fitting code on the same pinned inputs. They differ only in (i) inner Monte Carlo precision and (ii) how random numbers are assigned to the cross-validation that tunes λ. Neither changes the estimand or the expected value of the estimator. One residual cross-stream link in θ̂ has not been evaluated empirically. Its effect has a strict worst-case bound of 1.6×10⁻³ (about 3 % of the interval width), and a typical-size figure of about 2×10⁻⁵ (section 4.3). Both are analytical calculations, not measurements. No additional analysis is required by the plan. The disclosures in section 6 must accompany any report.

One governance point is corrected here. Earlier working notes, the report tool and the completion inventory attributed the retention of θ̂ to "Decision 3". In the owner's recorded approvals (configuration `approvals` and `owner_decisions_outside_schema`), **Decision 3 is B_in,max = 100 together with the RNG structural argument; it is not θ̂.** The retention-and-disclosure proposal for θ̂ is Decision 3 of the adjudication record, which is still a DRAFT with unticked boxes. θ̂'s authority is therefore the plan itself (section 2). The owner may wish to record an explicit approval of the disclosure wording, but the plan does not require one for reporting.

## 1. Source, calculation and provenance of θ̂

| Item | Evidence |
|---|---|
| Stored value | θ̂ = `monte_carlo.M2_primary.corrected_C_mean − monte_carlo.M1.corrected_C_mean` = **+0.004316342804030704**, in `reports/phase9/execution/phase9_execution_results.json` (pinned by the corrected input manifest) |
| Produced by | `scripts/38_phase9_execute.py` (Phase 9 primary execution, authorised 2026-11-03), B = 500 per seed, seeds 20261102–20261106 |
| Apparent part | Full pipeline fit on D for each model, fresh `default_rng(20261102)` per model (lines 90–92). C_app(M1) = 0.6242431994542509 (λ = 0), C_app(M2) = 0.6382067024814531 (λ = 0.3511). **ΔC_app(D) = +0.013963503027202173** |
| Optimism part | For each seed: 500 event-stratified patient resamples (`boot_indices`, lines 41–47). Each model is refit on the resample and evaluated on D; opt = C(resample) − C(D) (lines 119–137). Corrected C per seed = C_app − mean opt (line 145). θ̂ averages the five seeds. B_effective = 500 in every seed for both models, so this equals one pooled mean over 2,500 replicates: **ΔOpt(D) = +0.00964716022317149** |
| Identity | θ̂ = ΔC_app(D) − ΔOpt(D). Reconstructed: 0.004316342804030683; stored: 0.004316342804030704. The difference of 2×10⁻¹⁷ is floating-point summation order |
| Failures | 0 non-estimable replicates for every model (`non_estimable_replicates.rate` = 0.0 of 2,500 each) |

## 2. The estimand and point-estimate definition in the plan

`docs/PHASE9_DELTAC_INTERVAL_PLAN.md` v2.1:
- Line 6: "Estimand (unchanged): ΔC_corrected = C_corrected(M2 four-set) − C_corrected(M1) = +0.0043."
- Section 3.1, line 139: `theta_hat = DeltaC_corrected(D)  # +0.0043`. Lines 141–159 then apply the same steps to each outer resample D*ₖ: refit both models; compute ΔC_app(k); run B_in inner resamples of D*ₖ, evaluated on D*ₖ; apply pairing-preserving exclusion; ΔC_corrected(k) = ΔC_app(k) − ΔOpt(k).
- Section 3.2, lines 165–167: "The estimator applied to the original data is θ̂(D) = ΔC_app(D) − Δopt̂(D). The nested bootstrap applies **the identical functional** to each D*ₖ."
- Section 1 (lines 26–76) derives the inner noise from the same stored Phase 9 primary aggregates (5 seeds × 500) and sizes B_in = 100 on that basis. The plan therefore treats θ̂'s 2,500-replicate optimism and the 100-replicate inner optimism as the same quantity estimated with different Monte Carlo precision.

The corrected engine (`src/phase9_corrected/engine.py`, frozen) implements exactly these steps. Its docstring and lines 63–155 cite plan lines 138–161. Per outer replicate: `dC_app = C_app(M2) − C_app(M1)`; `dOpt = mean opt(M2) − mean opt(M1)` over paired-complete inner replicates; `dC_corrected = dC_app − dOpt`. The primary interval is the percentile interval of `dC_corrected(k)` (plan line 161), which does not use θ̂. θ̂ enters only the basic interval, z₀, the z₀-only BC diagnostic, A7 and A8 (plan lines 316 and onward).

## 3. Component-by-component comparison

| Component | θ̂ on D (`scripts/38`) | Each outer replicate D*ₖ (corrected engine) | Same functional? |
|---|---|---|---|
| Population and inputs | A-19 complete cases, 452 patients / 136 events; expression, covariates, endpoint, membership, panel | Same selection logic (`data.py:36–47`); the same files, verified by SHA-256 against the execution record's `input_hashes` (all 11 match) | Yes |
| Fitting code | `src.phase9_execution.pipeline_fit`, `c_index`, `linear_predictor` | The same functions, unchanged (`fitting.py`; frozen, hash-verified) | Yes |
| λ grid | `geomspace(1e-3, 10, 12)` | Identical list (checked element by element) | Yes |
| Models | M1 clinical; M2 clinical + frozen four-set panel | Same (`data.py:31–33`) | Yes |
| Resampling | Event-stratified patient bootstrap | Event-stratified patient bootstrap (`resampling.strat_boot`), outer and inner | Yes |
| Apparent ΔC | Full fit on D, C on D | Full fit on D*ₖ, C on D*ₖ | Yes (plan 3.2) |
| Optimism evaluation set | Refit on a resample of D, evaluated on D | Refit on a resample of D*ₖ, evaluated on D*ₖ | Yes (plan 3.2, first bullet) |
| Inner replicates | 2,500 (5 × 500) | B_in = 100 | Same expectation; different Monte Carlo precision (section 4.1) |
| Failure handling | Per-model mean over successful replicates | Pairing-preserving exclusion | Numerically identical on D: there were 0 failures |
| Apparent-fit CV random numbers | Fresh `default_rng(20261102)` per model | Fresh generator per model from a dedicated stream `[BASE, run_id, T]`, k (owner Decision 2) | Same procedure, different random numbers (section 4.2) |
| Inner-fit CV random numbers | The resample's generator, advancing after the resample draws | Dedicated child stream per (k, b) (plan line 364) | Same procedure, different random numbers (section 4.2) |

**Executed reproduction check** (`scripts/47_phase9_theta_alignment_check.py`, read-only, about one minute). It applies the corrected engine's `fit_and_score` to D under θ̂'s historical apparent convention. Result: C_app(M1) and C_app(M2) are **bit-identical** to the stored values, with λ identical (0.0 and 0.3511191734215131). It also confirms the identity θ̂ = ΔC_app(D) − ΔOpt(D) and the 11 input hashes. This shows that the corrected engine computes θ̂'s apparent component with the same code path. The optimism component cannot be compared replicate by replicate, because `scripts/38` stores only per-seed aggregates; recomputing it would be a new 2,500-replicate analysis (section 7).

## 4. The differences, and whether they matter

### 4.1 Inner Monte Carlo precision (2,500 versus 100)
Both are Monte Carlo averages of the same bootstrap optimism, so both target the same quantity. The 100-replicate version adds simulation noise to each ΔC_corrected(k). Plan section 1.2 quantifies this: Var(observed outer) = Var(sampling) + σ_Δ²/B_in, with σ_Δ ≈ 0.00943. That is about 1.2 % SD inflation at B_in = 100, conservative and deliberately not deconvolved (plan item 7, line 508). θ̂'s own inner noise is σ_Δ/√2500 ≈ 0.00019. The plan accepts this asymmetry by design. W3 monitors it (raise B_in if the inflation exceeds 5 %).

### 4.2 CV random-number assignment
The λ-tuning CV draws random fold assignments. Under both conventions the folds are pseudo-random draws, independent of the resample's content, so the distribution of the selected λ given the data, and hence the functional, is the same.

The **defect** that motivated the corrected rerun belongs to the as-run nested driver (`scripts/40`), not to θ̂. In the nested driver, the CV generator was **re-created from the resample's own seed**, so each CV replayed the bits that produced the resample (`CV_SEED_DEVIATION_DRAFT_v2.md` section 8). `scripts/38` never re-creates the generator: its CV draws come from bits **after** the resample (deviation record, section 8, "Different"). θ̂'s optimism replicates therefore do not have the replay link. They carry only a sequential dependence: the resample for replicate b+1 depends on how many draws replicate b's CV consumed. That makes replicates not a pure function of (seed, b), which matters for reproducibility bookkeeping. It does not bias a mean of pseudo-random replicates.

The apparent fit on D is one CV realisation. That is the estimator as defined: the plan applies "the same procedure" to D and D*ₖ, not the same random numbers (adjudication record, Decision 2, as approved by the owner). Each D*ₖ receives its own independent CV realisation, so the bootstrap distribution, and the interval, include λ-selection variability. θ̂'s reported MC-SE does not include it (disclosure 2).

### 4.3 The one unevaluated link, bounded
The apparent-CV generator on D, `default_rng(20261102)`, starts from the same state as the generator of the first resample of seed 20261102 (`scripts/38` lines 91 and 118). The apparent CV folds and that one resample are therefore functions of the same random bits. The link touches 1 of the 2,500 optimism replicates. Its influence on ΔOpt(D), and so on θ̂, is (Δopt_b − ΔOpt)/2500. Two analytical figures follow:
- **Strict worst-case bound.** Each per-model optimism lies in [−1, 1], so Δopt_b and ΔOpt each lie in [−2, 2], and |Δopt_b − ΔOpt| ≤ 4. The link therefore moves θ̂ by at most 4/2500 = **1.6×10⁻³**: about 3 % of the corrected interval widths (0.056–0.061), and smaller than the distance from θ̂ to either interval endpoint.
- **Typical-size figure, not a bound.** If the linked replicate behaves like the others (plan per-replicate SD σ_Δ ≈ 0.00943), a five-SD excursion moves θ̂ by about **1.9×10⁻⁵**, one tenth of θ̂'s MC-SE (0.00019). This assumes the link does not make that replicate atypical, which has not been tested.

Neither figure is an empirical measurement: `scripts/38` keeps no per-replicate values, so the actual effect was not measured. The link cannot affect the percentile interval, which does not use θ̂.

## 5. Comparison with the approved protocol and owner decisions

- **Plan v2.1** pre-specifies θ̂ = ΔC_corrected(D) = +0.0043 and states that the nested bootstrap applies the identical functional (lines 6, 139, 165–167). The corrected campaign runs under this plan (Decision 14 authorisation).
- **Owner Decision 1** (CV-seed deviation adjudicated; plan line 364 binding) and **owner Decision 2** (independent apparent-fit CV stream, T = 20261106) change only random-number assignment inside the outer replicates. The adjudication record's specification of that stream, which the owner approved in Decision 2, states that D is not an outer bootstrap replicate, that no new seed convention for D is introduced, and that the procedure, not the random numbers, must match between D and D*ₖ.
- **Draft Decision 3** (retain θ̂ = +0.004316342804030704, do not re-centre, mandatory MC-SE disclosure) has **not** been formally approved. Its substance is a restatement of the plan plus disclosures, and section 6 adopts those disclosures in full.
- **Deviation record, question 4** ("does θ̂'s third convention affect comparability in the basic interval and z₀?") is answered by sections 3–4. Comparability is affected only through the bounded link in 4.3 and through θ̂'s own Monte Carlo and CV-realisation variability. A7 and A8 are warnings only, and their decomposition separates the tail-asymmetry term, which does not depend on θ̂.

Nothing here changes the estimand, the interval method, a threshold or a seed.

## 6. Mandatory disclosures with the final result

1. θ̂ = +0.004316 is the plan's pre-specified point estimate, computed once on D by the Phase 9 primary execution with 2,500 optimism replicates. It was not re-estimated under the corrected random-number streams.
2. θ̂'s Monte Carlo SE (≈ 0.00019) excludes the randomness of the cross-validation on D that selected λ. The apparent ΔC on D is a single CV realisation, held fixed across the five seeds (`scripts/38` line 145).
3. θ̂ used the historical seeding convention (CV draws from one advancing generator after each resample). The apparent-CV generator on D shares its starting state with one of the 2,500 optimism resamples. Its effect on θ̂ was not measured. Analytically it is at most 1.6×10⁻³ (strict worst case), and about 2×10⁻⁵ if that replicate is typical (section 4.3).
4. The primary result is the five per-run 95 % percentile intervals, which do not depend on θ̂. The basic interval, z₀, the BC diagnostic, A7 and A8 depend on θ̂ and carry disclosures 1–3.
5. The interval is conditional on the observed 136 recurrence events and 316 censored observations in the 452-patient A-19 population.

## 7. Additional analysis

**None is necessary under the plan.** A re-estimate of θ̂'s optimism under dedicated CV streams, or across several CV seeds on D, would be a separately labelled sensitivity analysis. It would never replace θ̂ (draft Decision 3, last bullet). It is not authorised and not requested. It would refine disclosure 2; it is not needed to align the estimate with the interval.
