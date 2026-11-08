# Final report: optimism-corrected value of transcriptomic pathway scores for recurrence risk in stage II–III colon cancer (GSE39582)

**Date:** 8 November 2026. **Status:** final scientific report. It supersedes, for every finding, the provisional historical reports in `reports/phase9/final_project/documentation/`, which are kept unchanged.

**Revision (8 November 2026):** wording and citations only, following `docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md`. No number, model, endpoint or conclusion changed.

**Evidence base.**
- **Corrected ΔC results (sections 7–8):** `reports/phase9/corrected_rerun_report/corrected_interval_results.json` (SHA-256 `1c504b13…eb60`). It was produced by `scripts/46` only after the frozen post-run verification returned PASS with reporting permitted (`reports/phase9/corrected_rerun/verification/L00750.report.json`, `9558aaa5…ff25`).
- **Independent reproduction:** `reports/phase9/independent_audit/INDEPENDENT_AUDIT_REPORT.md` reproduced every number from the raw checkpoints by a separate calculation path (maximum difference 2.2×10⁻¹⁶). Classification: **INDEPENDENTLY REPRODUCED**. Here "independent" means a separate re-implementation within this project, run on the same raw checkpoints; it is not an external replication by another group or in another cohort.
- **Figures:** `docs/SCIENTIFIC_RESULTS_AND_FIGURES.md` presents these results with 9 main and 2 supplementary figures (`reports/figures/final/`), drawn by `scripts/49_make_results_figures.py` from the frozen outputs only. Figure map:

  | Section | Figure |
  |---|---|
  | §2 | Fig. 1, S2 |
  | §3 | Fig. 2 |
  | §4 | Fig. 4 |
  | §7 | Fig. 7 |
  | §8 | Fig. 8, S1 |
  | §9 | Fig. 5, 6 |
  | §10–12 | Fig. 9 |
  | Exploratory structure | Fig. 3 |

  The figures add no new result.
- **Coefficients, calibration and sensitivity results (sections 9–12):** the Phase 9 primary and S-7 execution records (`reports/phase9/execution/phase9_execution_results.json`, `05d13b8b…`; `reports/phase9/execution_s7/phase9_s7_results.json`). They are labelled as such where they appear.

> **Conditional on the observed 136 recurrence events and 316 censored observations in the 452-patient A-19 population.**

## Summary

| Item | Result |
|---|---|
| Question | Does a frozen four-pathway transcriptomic score set improve optimism-corrected discrimination for recurrence beyond clinical covariates? |
| Point estimate (pre-specified) | ΔC_corrected = **+0.0043** (C_corrected 0.6059 for M2 vs 0.6016 for M1) |
| 95 % percentile intervals (primary; five runs differing only in random seed, B_out = 750, B_in = 100) | [−0.0086, +0.0519] · [−0.0090, +0.0469] · [−0.0093, +0.0484] · [−0.0091, +0.0518] · [−0.0091, +0.0494] |
| Reproducibility across runs | endpoint SD 0.00025 (lower), 0.00217 (upper) |
| Monte Carlo precision | EP-1 met: largest endpoint MC-SE is 4.2 % of width |
| Verification | frozen verifier PASS, 60/60; an independent re-implementation within the project reproduced all values (not external replication) |
| Hard failures | none, so the interval is not labelled UNRELIABLE or IMPRECISE |
| Warnings | A7/A8 in every run: the bootstrap distribution sits above θ̂ and the percentile and basic intervals diverge. W4 once |
| Molecular coefficients | primary four-set model: all four pathway hazard-ratio intervals include 1 (historical, penalised fit). Five-set sensitivity model only: E2F HR 0.862 (0.746–0.996), with correlated scores, not a primary finding (§9.1) |
| Batch | molecular coefficients unstable under leave-one-batch-out refits (Wnt and p53 change sign); batch effects not resolved (§12) |
| Calibration | corrected slopes 0.83 (M1) and 0.86 (M2): residual overfitting (historical) |
| **Conclusion** | **No demonstrated improvement.** The data are compatible with no gain and with a gain of up to about 0.05 in C; the result is not evidence of equivalence. Validation is internal only; there is no claim of clinical utility or transportability |

## 1. Research question and rationale

**Question.** Does adding reproducible transcriptomic pathway information improve recurrence-risk prediction beyond clinical characteristics? Are any improvements stable enough to be statistically and biologically credible?

**Rationale.** Expression signatures are often reported to improve prognosis in colorectal cancer. Part of any apparent improvement is optimism: a larger model fits its own training data better than new patients. This project therefore:
- fixed a small, biologically motivated panel **before any outcome was examined**;
- estimated the gain in Harrell's C **after** optimism correction;
- put a confidence interval on that corrected gain with a nested bootstrap, so that the size **and** stability of any gain are quantified.

## 2. Cohort and endpoint

| Item | Specification |
|---|---|
| Data | GEO GSE39582, Affymetrix HG-U133 Plus 2.0 (GPL570); 585 arrays deposited after RMA normalisation and ComBat batch correction by the depositors |
| Cohort | Stage II–III primary colon cancer with documented recurrence follow-up: the S3 cohort of 453 tumours. QC retained all 453, with 29 soft-flagged and none excluded for lack of evidence of technical invalidity (Phase 7C) |
| Analysis population (A-19) | the 452 S3 patients with all five clinical covariates recorded: **136 recurrence events, 316 censored** |
| Endpoint | Register A-07 cause-specific recurrence, full follow-up (`rfs_event_full`, `rfs_time_full`) |
| Time unit (U1) | `rfs.delay` has no documented unit (provisional: months). Only rank-based quantities are reported, so the unit does not affect any result |
| Death without recurrence (U2) | no published rule. The deposited data censor these deaths (46 of 453 in S3, 10.2 %). Phase 6A §2.3 showed this censoring is **informative**: compared with patients alive without recurrence, the 46 were older (median age 75.0 vs 68.0), less often stage III (32.6 % vs 38.7 %) and less often treated with adjuvant chemotherapy (23.9 % vs 42.1 %), the profile of competing non-cancer mortality |
| Estimand implied by the endpoint | Because death is censored, every model estimates the **cause-specific hazard** of recurrence: the recurrence rate among patients still alive and recurrence-free. This differs from the **cumulative incidence** of recurrence, which treats death as a competing event (Fine–Gray / CIF), and from **composite disease-free survival**, which counts death as an event (182 events, a different estimand). Neither alternative was analysed (Phase 6A §3.2.2–3.2.3) |
| Events per coefficient | 136 events for 10 estimated M2 coefficients (≈ 13.6 per coefficient; the protocol's parameter budget gives EPV 15.0–15.1, `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md` §2.1). Adequate by convention, but overfitting is still visible (§10) |

## 3. Expression preprocessing and probe-to-gene mapping

- **No re-normalisation.** The deposited RMA + ComBat values were used as provided. The project applied no further normalisation, batch correction or expression filtering (Phases 2, 7 and 8).
- **Probe classes.** Of the 54,675 GPL570 probe sets:
  - 43,163 map to exactly one gene symbol and are used;
  - 8,893 are unmapped and excluded;
  - 2,619 are ambiguous (two or more symbols) and excluded, never assigned arbitrarily.

  This gives 21,815 genes reachable from uniquely mapped probes (Phase 2).
- **Gene-level representation** (Phase 2B, approved):
  - **median** across a gene's uniquely mapped probes is the **primary** representation;
  - a representative probe chosen by largest IQR on the 443-sample discovery partition (S-7) is the **sensitivity-only** representation.

  The median rule was chosen as uniform and reportable. It is not claimed to be biologically optimal.
- **Outcome-blind.** All preprocessing and mapping steps were run without opening any outcome column.

## 4. Pathway panel and scoring

- **Panel** (Phase 6C, frozen, outcome-blind). Four MSigDB Hallmark sets: E2F targets, epithelial–mesenchymal transition (EMT), p53 pathway and Wnt/β-catenin signalling. Their union is 607 genes (members resolved: 198/195/41/194).
- **Gene identity.** Entrez-primary mapping to GPL570, with symbol fallback.
- **Five-set sensitivity panel.** Adds PI3K/AKT/mTOR signalling (702 genes). It is sensitivity-only and may not be promoted on the basis of performance.
- **Score.** For pathway P, sample i:

  score_P(i) = mean over member genes g of (x_gi − μ_g)/σ_g

  - μ_g and σ_g (sample SD, ddof = 1; amendment M2-A1) are estimated on the **training samples of each fit only**, and genes with zero variance in that fold are dropped.
  - A score is dropped for a fold unless ≥ 80 % of source members resolve **and** ≥ 20 genes remain.
  - The four scores are then standardised again, on the training samples (amendment M2-A2). Each hazard ratio is therefore **per one training-fold SD of the pathway score**.
- **Implementation.** `src/pathway_scoring.py` is the single, hash-pinned implementation.

## 5. Model design

| Model | Covariates | Fitting |
|---|---|---|
| **M1 (clinical)** | age (restricted cubic spline, knots at training percentiles 10/50/90: 2 terms), sex, TNM stage III vs II, tumour location (proximal vs distal), adjuvant chemotherapy | Cox proportional hazards, unpenalised; 6 coefficients |
| **M2 (primary)** | M1 + the four standardised pathway scores | Cox; ridge penalty on the four molecular coefficients only; 10 coefficients |
| M2 sensitivity | M1 + five pathway scores | as M2; 11 coefficients; sensitivity-only |

- **Tuning.** The ridge λ is chosen from 12 values, log-spaced from 10⁻³ to 10, by five-fold **event-stratified** cross-validation, maximising the mean held-out C.
- **No leakage.** Every estimated quantity is re-estimated inside **every** fit from that fit's training samples only (`pipeline_fit`, `src/phase9_execution.py`): spline knots, gene μ/σ, score standardisation and λ.
- **Not fitted.** No model with transcriptome-wide feature selection (M3) was fitted.
- **Penalisation asymmetry.** Only the molecular block is penalised; the clinical block is not. The M2-vs-M1 comparison is therefore conservative by design, and ΔC is not an unbiased estimate of the molecular contribution. A conservative design does not turn a null-consistent result into a positive one (`docs/PHASE9_EXECUTION_REPORT.md` §15).

## 6. Validation design: nested cross-validation and nested bootstrap

**Three nested layers:**
1. **Inner cross-validation:** five-fold λ tuning inside every M2 fit.
2. **Bootstrap optimism correction** (Harrell). For each event-stratified bootstrap resample, refit the full pipeline, then compute optimism = C(resample) − C(original data). The corrected C is the apparent C minus mean optimism.
3. **Outer nested bootstrap for the interval** (plan v2.1 section 3.1).
   - Each outer replicate k resamples the 452 patients, stratified by event status.
   - It recomputes the apparent ΔC on D*ₖ.
   - It estimates optimism from B_in = 100 inner resamples of D*ₖ, evaluated on D*ₖ, with M1 and M2 paired on the same resamples and dropped together on failure.
   - ΔC_corrected(k) = ΔC_app(k) − ΔOpt(k).

**Primary interval:** the 95 % percentile interval in each of **five runs that differ only in random seed** (Monte Carlo replicates of one analysis of one cohort, not independent cohorts or experiments).

**Random numbers and sizing:**
- Every random draw has its own stream: outer resample, inner resample, inner-CV and apparent-CV, all from NumPy `SeedSequence`.
- This corrects a deviation in the earlier as-run driver, where cross-validation reused the resample's random state (owner Decisions 1–2).
- Runs were completed at B_out = 500, then topped up to 750 because the frozen EP-1 sizing rule required it (ρ_max = 0.059 at 500). At 750 EP-1 is met and the sizing decision is terminal.

**Point estimate.** θ̂ = +0.004316342804030704 is the plan's pre-specified estimate, ΔC_corrected on the full data (2,500 optimism replicates).
- `docs/PHASE9_THETA_ALIGNMENT_AUDIT.md` shows it is the same functional the interval bootstraps; its apparent component reproduces bit for bit.
- Each outer replicate uses B_in = 100 inner resamples instead of 2,500. The estimator is the same; only the inner Monte Carlo precision differs.
- Its Monte Carlo SE (≈ 0.00019) excludes the cross-validation randomness on the full data.
- One historical RNG link (the apparent-CV generator on the full data shares its starting state with one of the 2,500 optimism resamples) was not measured. Its strict analytical worst-case effect on θ̂ is **1.6×10⁻³** (each per-replicate optimism difference lies in [−2, 2], so one replicate moves a 2,500-replicate mean by at most 4/2500). The typical-size estimate is about **1.9×10⁻⁵**, conditional on the linked replicate behaving like the others. **Neither figure is an empirical measurement**; `scripts/38` stores no per-replicate values (`docs/PHASE9_THETA_ALIGNMENT_AUDIT.md` §4.3).
- The primary percentile interval does not use θ̂ and cannot be affected by this link. The basic interval, z₀, the BC diagnostic, A7 and A8 do use θ̂.

**Engineering controls:**
- frozen code identity: 27 files, SHA-256 `618f20eb…a790`;
- a pre-execution gate;
- append-only, hash-linked checkpoints;
- a frozen post-run verifier (section 7.1);
- an independent audit (section 7.2).

## 7. Corrected discrimination results

### 7.1 Frozen verification (B_out = 750)
**PASS: 60 of 60 checks**, with 0 FAIL, 0 NOT_APPLICABLE and 0 NOT_COMPARABLE. Reporting is permitted, and H5 is false.

The checks covered:
- invocation binding and run eligibility;
- the sizing decision at 500 and 750;
- resample regeneration (750 outer and 75,000 inner per run, 0 mismatches);
- stream assignment;
- clinical-model invariance;
- the as-run comparison with history (500 per run);
- a full bit-identical recomputation of 20 replicates.

### 7.2 Independent audit
Every interval quantity was recomputed from raw checkpoint bytes with a standard-library-only implementation:
- 3,750 outer replicates, 375,000 inner pairs, 0 fit failures;
- all per-record hashes valid;
- the original 500-record prefixes byte-identical;
- 138 reported quantities reproduced within 2.2×10⁻¹⁶.

Classification: **INDEPENDENTLY REPRODUCED**.

This classification means the calculations were re-implemented separately within this project and reproduced from the same raw checkpoints. It is not an external replication: no other group, cohort or dataset was involved.

### 7.3 Corrected C by model
**Full-data estimates** (Phase 9 primary; these define θ̂):

| | M1 | M2 | M2 five-set |
|---|---:|---:|---:|
| Apparent C | 0.6242 | 0.6382 | 0.6413 |
| Optimism-corrected C | **0.6016** | **0.6059** | 0.6066 |

Both models discriminate modestly. Optimism removes about 0.023 (M1) and 0.032 (M2) of apparent C. M2's larger optimism consumes most of its apparent advantage of +0.0140.

## 8. Corrected ΔC estimate and interval

| Run | 95 % percentile interval (primary) | Width | Median ΔC_corr(k) | Share of replicates < 0 |
|---|---|---:|---:|---:|
| 0 | [−0.00862, +0.05189] | 0.0605 | +0.0120 | 0.225 |
| 1 | [−0.00896, +0.04689] | 0.0559 | +0.0123 | 0.188 |
| 2 | [−0.00929, +0.04837] | 0.0577 | +0.0119 | 0.189 |
| 3 | [−0.00915, +0.05178] | 0.0609 | +0.0114 | 0.231 |
| 4 | [−0.00908, +0.04944] | 0.0585 | +0.0128 | 0.179 |

**Point estimate:** θ̂ = **+0.0043**.

**Between-run reproducibility** (df = 4):

| Endpoint | SD (ddof = 1) | Range |
|---|---:|---:|
| Lower | 0.00025 | 0.00067 |
| Upper | 0.00217 | 0.00499 |

**Endpoint MC-SE** (kernel, maximum over 0.5h/h/2h):
- lower endpoints: 0.00044–0.00054;
- upper endpoints: 0.00150–0.00237;
- EP-1 ratio ≤ 0.042.

**Pooled draws (descriptive only; the plan designates no pooled interval):**
- n = 3,750; percentile range [−0.00906, +0.05023];
- median +0.0121; 20.2 % of draws below zero.

**Hard failures H1–H6:** none.
- 0 of 3,750 outer replicates excluded;
- 0 inner fit failures, so the bounding analysis is vacuous;
- EP-1 met, frozen identity intact, V5 bit-identical, paired counts equal.

**Warnings** (plan section 11.2; investigate and report):

| Warning | Result |
|---|---|
| W6/A7 | 0.525–0.572 in every run (threshold 0.25) |
| W7/A8 = \|z₀\| | 0.431–0.591 in every run (threshold 0.25) |
| W4 | once (run 4 lower endpoint, kernel/between-run SE ratio 2.01) |
| W3 | inner-noise SD inflation 0.33–0.39 %, under the 5 % threshold |
| W1, W5 | not triggered |

**Diagnostic intervals (θ̂-dependent):**

| Run | Basic | BC (z₀ only) | z₀ | A7 (tail asymmetry + median offset) |
|---|---|---|---:|---|
| 0 | [−0.0433, +0.0173] | [−0.0112, +0.0299] | −0.475 | 0.572 (+0.0096, +0.0077) |
| 1 | [−0.0383, +0.0176] | [−0.0116, +0.0249] | −0.548 | 0.525 (+0.0067, +0.0080) |
| 2 | [−0.0397, +0.0179] | [−0.0118, +0.0269] | −0.544 | 0.528 (+0.0077, +0.0075) |
| 3 | [−0.0432, +0.0178] | [−0.0117, +0.0319] | −0.431 | 0.558 (+0.0099, +0.0071) |
| 4 | [−0.0408, +0.0177] | [−0.0117, +0.0265] | −0.591 | 0.542 (+0.0073, +0.0085) |

**What the warnings show:**
- **Location.** The bootstrap distribution sits **above** θ̂: about a third of replicates fall below it, and the median exceeds θ̂ by about 0.008.
- **Mechanism.** Replicates show a larger apparent ΔC (mean +0.022 vs +0.014 on the full data) and a smaller optimism difference (+0.0082 vs +0.0096).
- **Contrast with the plan.** This is the opposite direction to the dilution mechanism anticipated in plan section 7.2, and the analysis does not explain it.
- **Consequence.** Per plan section 7, a material percentile–basic divergence warns that neither interval should be taken at face value. Per section 7.1, it neither validates nor invalidates the interval.

**The qualitative result is robust to this:** the percentile, basic and bias-corrected intervals **all include zero in all five runs** (15 of 15).

### 8.1 A distinct quantity: the paired out-of-original ΔC (historical)
The Phase 9 primary execution also reported a **paired out-of-original ΔC**: the difference in C on the original data between models fitted on each bootstrap resample. It was +0.0121, bootstrap 95 % percentile [+0.0002, +0.0240], under the median representation (`docs/PHASE9_EXECUTION_REPORT.md` §7), and +0.0105 [−0.0001, +0.0224] under S-7 (`docs/PHASE9_S7_SENSITIVITY_REPORT.md` §4).
- It is **not** the optimism-corrected ΔC. It mixes in-sample and out-of-sample information, and its percentiles describe variation of bootstrap-fitted models, not the sampling uncertainty of ΔC_corrected.
- Its interval narrowly excluded zero under one representation and included zero under the other: a knife-edge, not evidence.
- It must not be substituted for, or transferred to, the corrected-ΔC intervals above, which are the only intervals for the estimand.

## 9. Molecular coefficient estimates and uncertainty (historical, full-data M2 fit)

Hazard ratio per one training-fold SD of each score, from the ridge-penalised M2 fit (λ = 0.351):

| Pathway | HR (95 % Wald CI), median representation | HR (95 % Wald CI), S-7 representative-IQR |
|---|---|---|
| E2F targets | 0.906 (0.807–1.016) | 0.872 (0.758–1.003) |
| EMT | 1.070 (0.951–1.202) | 1.084 (0.937–1.254) |
| Wnt/β-catenin | 1.016 (0.906–1.138) | 1.011 (0.879–1.162) |
| p53 pathway | 0.995 (0.882–1.123) | 0.993 (0.852–1.158) |

**Clinical covariates, for reference:** in M1, stage III vs II HR 1.97 (1.30–2.98) and male sex HR 1.38 (0.97–1.96). In M2, stage HR 1.94 (1.28–2.94); with the nine batches as strata it attenuates to 1.68 (1.07–2.65) (coefficient 0.662 → 0.520; `docs/PHASE9_EXECUTION_REPORT.md` §10). The clinical baseline's strongest term is therefore itself sensitive to batch structure.

**Interpretation:**
- **No independent association shown (primary model).** Every pathway interval in the primary four-set model includes 1, under both representations, so **no pathway shows an independent association** with recurrence after adjustment.
- **Ridge shrinkage.** The penalty shrinks the coefficients towards zero, so the Wald intervals are approximate and descriptive. They are not valid post-selection inference.
- **Direction.** E2F below 1 and EMT above 1 keep their direction under S-7, under batch stratification and in every leave-one-batch-out refit (§12). Their magnitudes do not stay stable, and neither is distinguishable from zero.
- **Association, not discrimination.** Hazard ratios describe adjusted association. Added predictive discrimination is assessed by ΔC (§8), and the two need not agree.

### 9.1 Five-set sensitivity model (historical; sensitivity-only)
Ridge-penalised, λ = 0.152 (`docs/PHASE9_EXECUTION_REPORT.md` §5):

| Pathway | HR (95 % Wald CI) |
|---|---|
| E2F targets | **0.862 (0.746–0.996)** |
| EMT | 1.094 (0.948–1.263) |
| Wnt/β-catenin | 1.025 (0.894–1.175) |
| p53 pathway | 0.978 (0.832–1.151) |
| PI3K/AKT/mTOR | 1.043 (0.891–1.221) |

- **The E2F interval in the five-set model excludes 1, only just.** This is a **sensitivity-model result, not a primary finding**.
- Adding PI3K introduced correlated scores:
  - p53 ↔ PI3K r = 0.636;
  - condition number 2.9 → 7.1;
  - every molecular SE inflated ×1.20–1.34.
- The E2F coefficient moved from −0.099 to −0.149 in step with that collinearity. The narrowly excluding interval coincides with this collinearity-driven shift rather than reflecting a stable effect.
- The panel was fixed as sensitivity-only before any fit (Phase 6C), and its corrected C (0.6066) does not differ from the primary model (0.6059) beyond Monte Carlo error. It is not promoted.

### 9.2 Biological interpretation: plausibility versus evidence
A Hallmark score here is the mean training-standardised **mRNA** abundance of a fixed gene set in bulk tumour tissue (§4). It is **not** a measure of pathway activation, protein activity, signalling flux or mutation status, and no causal mechanism can be inferred from it. The background column below summarises published colorectal-cancer literature (references R1–R8 at the end of this report) as **plausibility only**. It was not tested in this project.

| Pathway | Background (plausibility, not tested here) | Demonstrated by this dataset |
|---|---|---|
| EMT | Mesenchymal and stromal expression programmes characterise the CMS4 consensus subtype, which had worse relapse-free survival (R1). GSE39582 was one of the data sets used to derive the CMS, so this is not independent evidence. Mesenchymal-subtype transcripts in bulk colorectal tumours are mostly of stromal origin (R2, R3); whether this holds for the Hallmark EMT score was not tested here | HR 1.07 (0.95–1.20), consistently > 1; not distinguishable from 1. The score cannot separate stroma from tumour-cell EMT |
| E2F targets | Proliferation and cell-cycle genes. Evidence on proliferation and colorectal-cancer outcome is mixed: one gene-signature study found lower proliferation in relapsed primary tumours (R4), whereas a meta-analysis of Ki-67 staining associated high proliferation with poorer survival (R5). Neither measured E2F-target mRNA | HR 0.91 (0.81–1.02), consistently < 1; includes 1 in the primary model. The single excluding interval is in the collinear five-set model (§9.1) |
| Wnt/β-catenin | Wnt-pathway genes are genomically altered in about 93 % of colorectal cancers, with APC or CTNNB1 alterations in about 80 % (R6). Genomic alteration is not the same as pathway activation, and Wnt-related expression varies between tumours (marked WNT/MYC activation characterises CMS2; R1). Target-gene expression is an indirect proxy | HR 1.02 (0.91–1.14); changes sign under leave-one-batch-out; least representation-stable score (r = 0.857). No signal |
| p53 pathway | The gene set captures p53 **transcriptional targets**. TP53 is mutated in about 60 % of non-hypermutated colorectal cancers (R6), but TP53 mutation is not equivalent to low target expression, and TP53 status was not analysed here | HR 1.00 (0.88–1.12); changes sign twice under leave-one-batch-out. No signal, and nothing can be said about TP53 status |
| PI3K/AKT/mTOR (sensitivity only) | PIK3CA is mutated in about 15 % of non-hypermutated colorectal cancers (R6); gene expression is not pathway activity | HR 1.04 (0.89–1.22); collinear with p53; mainly destabilised other coefficients |

The consistent E2F (< 1) and EMT (> 1) directions are **hypotheses consistent with part of the prior literature** (EMT/stroma: R1–R3; proliferation: mixed, R4 vs R5). This dataset demonstrates neither an association nor a mechanism.

## 10. Calibration and overfitting (historical, single-layer bootstrap)

| Model | Apparent slope | Optimism-corrected slope |
|---|---:|---:|
| M1 | 1.000 | **0.833** |
| M2 | 1.083 | **0.862** |
| M2 five-set | 1.075 | 0.836 |

- **Residual overfitting.** Corrected slopes below 1 indicate residual overfitting: predictions applied to new patients would be too extreme and would need shrinkage of about 14–17 %.
- **M1 included.** This applies to the clinical model too, because with 136 events even six coefficients overfit measurably.
- **Log-hazard slope only.** This is a slope on the log-hazard scale. Calibration-in-the-large and absolute-risk calibration cannot be assessed while the time unit is undocumented (U1).
- **Not re-estimated.** Calibration was not re-estimated under the corrected random-number streams.

## 11. Median versus representative-IQR sensitivity (S-7, historical)

| | Median (primary) | Representative-IQR (S-7) |
|---|---:|---:|
| Apparent C, M2 | 0.6382 | 0.6387 |
| Corrected C, M2 / M1 | 0.6059 / 0.6016 | 0.6059 / 0.6016 |
| Corrected ΔC | +0.0043 | +0.0043 |
| Selected λ | 0.351 | 0.152 |
| Corrected slope, M2 | 0.862 | 0.832 |
| Corrected C, five-set (sensitivity) | 0.6066 | 0.6041 (−0.0025, exceeds MC-SE 0.00041) |
| p53 ↔ PI3K score correlation (five-set model) | 0.636 | 0.688 |

- **Discrimination.** The corrected gain is the same under both representations.
- **Coefficients.** Coefficient directions agree (section 9). Magnitudes differ modestly, mainly for E2F and through the different λ.
- **Wnt.** Wnt is the score whose two representations correlate least (r = 0.858, Phase 8C), but its coefficient is near zero under both.
- **Five-set model.** It is the only model whose corrected C moves beyond Monte Carlo error under S-7. Its collinearity is worse under representative-IQR: condition number 7.08 → 9.55, and a second correlated pair (E2F ↔ PI3K, r = 0.545). The PI3K coefficient falls from +0.0417 to +0.0103 (`docs/PHASE9_S7_SENSITIVITY_REPORT.md` §4, §8). This reinforces its sensitivity-only role.
- **Not replication.** S-7 re-summarises the same arrays for the same patients. Agreement shows robustness of the pipeline to probe summarisation, not independent confirmation (`docs/PHASE9_S7_SENSITIVITY_REPORT.md` §10).
- **Scope.** S-7 has no nested interval; its interval is not applicable under the plan.

## 12. Batch-related instability

- **Confounding.** Processing batch (nine ComBat batches) is confounded with the discovery/validation partition (Cramér's V = 0.49; Phase 4/5). Residual batch-associated structure remains in the deposited, ComBat-corrected matrix (Phase 8B: ComBat batch associated with PC3, ε² = 0.061, BH q = 2.8×10⁻⁵).
- **Leave-one-batch-out refits of M2** (historical; nine refits, 0 failures). Apparent C ranges from 0.614 (batch03 omitted, n = 102) to 0.659 (batch04 omitted, n = 103), against 0.638 on all data. Apparent performance therefore depends visibly on which batches are included. Molecular coefficients (log-hazard per training-fold SD) are markedly unstable:

  | Term | All data | Leave-one-batch-out range | Sign changes |
  |---|---:|---|---:|
  | EMT | +0.067 | +0.009 to +0.204 | 0 |
  | E2F targets | −0.099 | −0.249 to −0.014 | 0 |
  | Wnt/β-catenin | +0.015 | −0.002 to +0.067 | 1 |
  | p53 pathway | −0.005 | −0.113 to +0.014 | 2 |
  | Stage III (clinical) | +0.662 | +0.498 to +0.826 | 0 |

  Sources: `reports/phase9/execution/phase9_execution_results.json` (`batch_sensitivity.leave_one_batch_out`); `docs/PHASE9_EXECUTION_REPORT.md` §10. Batches differ about seven-fold in size (14–103 patients, 4–52 events), and omitting batch03 or batch04 removes about 23 % of the cohort. The spread therefore conflates batch identity with sample size, as the protocol caveat anticipated.
- **Batch-stratified fit** (nine strata). It converged with a positive-definite covariance. The molecular coefficients barely moved (Δ ≤ 0.027), but the clinical terms moved more: stage 0.662 → 0.520 (HR 1.94 → 1.68), chemotherapy 0.137 → 0.029 (`docs/PHASE9_EXECUTION_REPORT.md` §10).
- **Implication.** Coefficient values, and in two cases their signs, depend on cohort composition. Together with the partition confounding and residual batch structure, part of the pathway-score variation may be technical. Coefficient-level biological interpretation is therefore not credible. **Batch effects have not been resolved**, only characterised.
- **Not performed.** The batch-perturbation sensitivity analysis (Phase 8D candidate) was not performed.
- **Consequence.** Residual batch effects could bias the pathway scores in either direction.

## 13. Internal-validation limitations

1. **Internal only.** Optimism correction estimates performance in new samples **from the same population and processing**. It is not external validation. **No external validation was performed.**
2. **CIT partition is not a test set.** It is confounded with batch and has 29 events, so it is at most an internal holdout and was not used (Phase 5).
3. **No feasible external cohort.** Phase 5 found that no available external cohort (GSE14333: 48 events; GSE33113: 18 events) could support confirmatory external validation.
4. **Conditional interval.** The interval is conditional on the observed event and censoring counts.
5. **Five runs.** They bound Monte Carlo noise, not systematic error common to all runs.
6. **No absolute risk.** The cause-specific Cox model gives relative hazards and rank-based discrimination. It does **not** give absolute recurrence probability or cumulative incidence (U1, U2).
7. **Partly historical.** θ̂, calibration, coefficients, S-7, batch analyses and the five-set results come from the historical single-layer Phase 9 primary execution, without nested intervals, and were not re-estimated under corrected streams. Only ΔC_corrected has a verified nested interval. The as-run nested intervals from the seed-defective driver are superseded and excluded from the findings.
8. **Competing risks (U2).** The cause-specific hazard is a valid aetiological quantity, but it cannot be converted into an absolute probability of recurrence. 1 − Kaplan–Meier would overstate recurrence risk here, and calibrating cause-specific predictions against observed cumulative incidence would be an estimand mismatch. No Fine–Gray or composite-endpoint analysis was performed (Phase 6A §3.2.3).
9. **Limits of the C-index.** C measures ranking only:
   - it is insensitive to added predictors when the baseline model is moderate;
   - it carries no information on calibration or on the consequences of decisions.
   - No net-benefit (decision-curve) or reclassification analysis was performed, and time-dependent discrimination cannot be computed while U1 is unresolved.
   - A ΔC near zero therefore says little about clinical usefulness either way.
10. **Limits of calibration.** Only log-hazard calibration slopes (single-layer, no intervals) were estimated. Calibration-in-the-large and absolute-risk calibration are not possible (U1, U2).
11. **Events.** 136 events (about 13.6 per M2 coefficient) limit precision and drive the width of the ΔC interval. Overfitting is visible even in the six-coefficient clinical model.

## 14. Interpretation and conclusion

**Observed:**
- The pre-specified optimism-corrected gain in C from adding four pathway scores is +0.0043.
- Five verified 95 % percentile intervals, from runs that differ only in random seed and reproduced by an independent re-implementation within the project, span about −0.009 to +0.05, and all include zero.
- No pathway coefficient in the primary model is distinguishable from no effect. The only excluding interval (E2F, 0.746–0.996) is in the collinear five-set sensitivity model.
- Molecular coefficients are unstable across batches (§12).
- Both models show residual overfitting.

**Statistical interpretation:**
- The study **does not demonstrate** that the pathway scores improve discrimination beyond clinical covariates in this cohort.
- A loss of up to about 0.009 and a gain of up to about 0.05 in C are both compatible with the data. An interval that contains zero does **not** establish equivalence or absence of a pathway contribution.
- ΔC is a measure of **added discrimination**, not of causal molecular contribution. A marginal gain would not establish a robust pathway contribution in any case.
- The interval is precise as a Monte Carlo estimate. Its width reflects sampling uncertainty with 136 events.
- Bootstrap bias diagnostics flag the interval's location relative to θ̂. All interval constructions agree on the qualitative result.

**Stability:**
- **Stable** across Monte Carlo runs, across the two expression representations, and under the seed correction: the corrected intervals are close to the superseded as-run ones (pooled range [−0.0091, +0.0510]).
- **Visibly sensitive** to batch composition.
- Taken together, the evidence does not meet the standard of a stable, credible pathway contribution.

**Biological hypothesis (not a finding):**
- EMT/stromal programmes are plausible recurrence correlates (R1–R3); evidence on proliferation (E2F) is mixed (R4, R5) (§9.2).
- In this dataset their estimated effects are small, imprecise and not distinguishable from zero, and their magnitudes vary with batch composition.
- Testing them would need larger, externally independent cohorts with documented endpoints, ideally with orthogonal measures such as stromal content, mutation data or protein assays.

**Status of the open items listed in the historical execution report** (`docs/PHASE9_EXECUTION_REPORT.md` §15, "Ready for final interpretation? No"):
1. **S-7 representation sensitivity:** since executed (`docs/PHASE9_S7_SENSITIVITY_REPORT.md`); the corrected ΔC is unchanged (§11).
2. **No interval for the corrected ΔC:** resolved by the verified nested bootstrap, reproduced by an independent re-implementation within the project (§8).
3. **Coefficient instability under batch perturbation:** characterised (§12), but **not explained or resolved**. It limits coefficient interpretation, not the ΔC interval.

**Not claimed:**
- statistical significance of a gain;
- clinical utility;
- improved patient stratification;
- absolute risk prediction;
- external validity or transportability;
- any causal mechanism.

## 15. Reproducibility and audit trail

| Artifact | Path | SHA-256 |
|---|---|---|
| Freeze record (27 files) | `config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json` | `618f20eb…a790` |
| Run configuration | `config/analysis/PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json` | config hash `bf939a4e…6a56` |
| Input manifest | `reports/phase9/corrected_rerun/manifest/input_manifest.json` | `1ff24387…3f27` |
| Interval plan | `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` | v2.1 |
| Modelling protocol | `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md` | v3.1 |
| Sizing decision at 500 | `reports/phase9/corrected_rerun/sizing/level00500.json` | `28c689bb…cf03` |
| Checkpoints (750 records each) | `reports/phase9/corrected_rerun/runs/run{0..4}.jsonl` | `0e27494f…`, `9dfe640b…`, `0204da98…`, `05274f0e…`, `af753007…` |
| Verification at 500 (FAIL on sizing; superseded by the top-up) | `reports/phase9/corrected_rerun/verification/L00500.report.{json,txt}` | — |
| Verification at 750 (PASS) | `reports/phase9/corrected_rerun/verification/L00750.report.{json,txt}` | `9558aaa5…ff25` |
| Interval results and tables | `reports/phase9/corrected_rerun_report/` | `1c504b13…eb60` (JSON) |
| Independent audit | `reports/phase9/independent_audit/` | `b67560a8…e58d` (JSON) |
| θ̂ alignment audit | `docs/PHASE9_THETA_ALIGNMENT_AUDIT.md` | — |
| Archive manifest | `reports/phase9/ARCHIVE_MANIFEST_SHA256.txt` | see `docs/LARGE_ARTIFACT_ARCHIVE.md` |

**Note on decision numbers.** The results file labels its pooled-draws note "adjudication Decision 9, item 5". That number follows the draft adjudication record, not the owner's recorded approvals, where Decision 9 is the verification gate. The reporting format here follows plan sections 2.3, 7 and 9.1 step 5.

**Compute.**
- Each run: 3.3–5.1 h at B_out = 500, plus 2.1–2.4 h for the top-up to 750, on 8 worker processes.
- 153,000 model fits per run (750 × [2 apparent + 2 as-run apparent + 200 inner]); each M2 fit adds 60 cross-validation fits for λ.
- Verification at 750: about 41 min.

## References

- R1. Guinney J, Dienstmann R, Wang X, et al. The consensus molecular subtypes of colorectal cancer. *Nat Med* 2015;21(11):1350–1356. doi:[10.1038/nm.3967](https://doi.org/10.1038/nm.3967)
- R2. Calon A, Lonardo E, Berenguer-Llergo A, et al. Stromal gene expression defines poor-prognosis subtypes in colorectal cancer. *Nat Genet* 2015;47(4):320–329. doi:[10.1038/ng.3225](https://doi.org/10.1038/ng.3225)
- R3. Isella C, Terrasi A, Bellomo SE, et al. Stromal contribution to the colorectal cancer transcriptome. *Nat Genet* 2015;47(4):312–319. doi:[10.1038/ng.3224](https://doi.org/10.1038/ng.3224)
- R4. Anjomshoaa A, Nasri S, Humar B, et al. Slow proliferation as a biological feature of colorectal cancer metastasis. *Br J Cancer* 2009;101(5):822–828. doi:[10.1038/sj.bjc.6605229](https://doi.org/10.1038/sj.bjc.6605229)
- R5. Luo ZW, Zhu MG, Zhang ZQ, et al. Increased expression of Ki-67 is a poor prognostic marker for colorectal cancer patients: a meta analysis. *BMC Cancer* 2019;19:123. doi:[10.1186/s12885-019-5324-y](https://doi.org/10.1186/s12885-019-5324-y)
- R6. Cancer Genome Atlas Network. Comprehensive molecular characterization of human colon and rectal cancer. *Nature* 2012;487(7407):330–337. doi:[10.1038/nature11252](https://doi.org/10.1038/nature11252)
- R7. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst* 2015;1(6):417–425. doi:[10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)
- R8. Marisa L, de Reyniès A, Duval A, et al. Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. *PLoS Med* 2013;10(5):e1001453. doi:[10.1371/journal.pmed.1001453](https://doi.org/10.1371/journal.pmed.1001453)

All references were verified on 8 November 2026 through Europe PMC (`docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md` §1).
