# Pre-specified pathway-expression scores and recurrence discrimination in stage II–III colon cancer: an optimism-corrected evaluation in GSE39582

**Prasad Akula** · November 2026

Audience: scientific reviewers and interviewers who want the complete argument in one document. Every number is taken from frozen project outputs. Source files are named in brackets; the full audit trail is in `docs/PHASE9_CORRECTED_FINAL_REPORT.md`.

---

## Abstract

**Background.**
- Gene-expression signatures are often reported to improve recurrence prediction in colorectal cancer.
- Apparent improvements are inflated by optimism: a larger model fits its own training data better than new patients.

**Question.** Do four biologically selected tumour-expression pathway scores add reproducible recurrence-discrimination information beyond established clinical variables in GSE39582?

**Methods.**
- **Cohort:** 452 patients with stage II–III colon cancer (136 recurrences, 316 censored) from GEO GSE39582, profiled on Affymetrix HG-U133 Plus 2.0 arrays.
- **Pathway scores:**
  - four MSigDB Hallmark sets (EMT, E2F targets, Wnt/β-catenin, p53 pathway; 607 genes), chosen without outcome data and frozen before any association analysis;
  - each scored as the mean of member-gene z-scores standardised inside each training fold.
- **Models:**
  - a clinical Cox model (age spline, sex, stage, tumour location, adjuvant chemotherapy);
  - the same model plus ridge-penalised pathway scores.
- **Validation:** an optimism-corrected nested bootstrap with five Monte Carlo runs differing only in random seed, each with 750 outer × 100 inner resamples.

**Results.**
- **Optimism:** the apparent gain in Harrell's C (+0.0140) was mostly optimism; the corrected gain (ΔC) was +0.004316.
- **Uncertainty:** all five 95 % percentile intervals included zero (about −0.009 to +0.05).
- **Pathways:** no pathway hazard-ratio interval excluded 1 in the primary model.
- **Calibration:** corrected slopes were 0.833 (clinical) and 0.862 (pathway model), indicating residual overfitting.
- **Batch and subtype:** pathway coefficients changed sign when single processing batches were omitted, and all four scores were associated with mismatch-repair (MMR) status in exploratory tests.

**Conclusion.**
- The study did not demonstrate a reproducible improvement in recurrence discrimination from adding the four pathway scores.
- This does not establish that these pathways have no role in recurrence.
- Interpretation is limited by:
  - internal validation only;
  - an undocumented recurrence-time unit;
  - informative censoring of deaths;
  - unmeasured cell composition;
  - unresolved batch structure.

---

## 1. Introduction

Recurrence after curative surgery for stage II–III colon cancer varies widely, and stage alone sorts that risk imperfectly. Tumour gene-expression profiles summarise biology that clinical variables cannot capture, and many signatures have been reported to add prognostic value.

Two problems recur in such reports:
- **Optimism.** Performance measured on the data used to build a model overstates performance in new patients, and the overstatement grows with model size.
- **Selection.** Selecting genes from tens of thousands of candidates after looking at outcomes produces signatures that fit noise. In this cohort, a simulation of transcriptome-wide selection under no true signal gave an apparent C of 0.92 when the true value was 0.50 (`docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md` §6.2).

This study tests a deliberately narrow question:
- four biologically motivated pathway programmes;
- fixed before outcomes were examined;
- added to a fixed clinical model;
- with the gain measured after optimism correction and given a valid uncertainty interval.

## 2. Research question

Do four biologically selected tumour-expression pathway scores add reproducible recurrence-discrimination information beyond established clinical variables in GSE39582?

The question was originally broader, a test of transportability to independent cohorts. It was narrowed after a feasibility assessment showed that no public cohort had enough recurrence events (section 3.7).

## 3. Methods

### 3.1 Data and cohort
- **Source.** GEO GSE39582 (Marisa et al., 2013, R8): 585 Affymetrix HG-U133 Plus 2.0 arrays, deposited as RMA-normalised and ComBat-batch-corrected expression. Expression was used as deposited, with no re-normalisation or batch correction.
- **Cohort.** Stage II–III primary tumours with documented recurrence data.
  - The primary cohort contained 453 patients.
  - The analysis population (complete clinical covariates) contained **452 patients, 136 recurrences and 316 censored**.
  - One patient was excluded for missing adjuvant-chemotherapy status.
  [`reports/figures/final/figure_values.json`, `fig1`]
- **Partitions.** The deposit includes a discovery set (353 patients, 107 events) and a CIT partition (99 patients, 29 events). The partitions were analysed together. The CIT partition shares processing history with the discovery set and is confounded with batch, so it was treated as **internal, never external**, validation.

### 3.2 Endpoint
- **Estimand.** The models estimate the **cause-specific hazard of recurrence**.
- **Censoring of deaths.** 46 deaths without recurrence are censored in the deposit. Those patients were older (median 75.0 vs 68.0 years), less often stage III (32.6 % vs 38.7 %) and less often treated with chemotherapy (23.9 % vs 42.1 %), so the censoring is informative.
- **Consequences.**
  - The cause-specific hazard remains estimable but cannot be converted into an absolute recurrence probability.
  - The recurrence-time unit is undocumented (a provisional assumption of months was recorded but not relied on), so only rank-based quantities are reported.
  - No competing-risk or composite-endpoint analysis was performed.

### 3.3 Probe-to-gene mapping
- 54,675 probe sets.
- 43,163 map to exactly one gene symbol (8,893 unmapped and 2,619 ambiguous probe sets excluded).
- These give 21,815 unique genes.
- Genes with several probes (10,751) are summarised by the median.
- A representative-probe rule (highest-IQR probe per gene, selected on discovery samples) was frozen as a sensitivity analysis, called S-7.

[`figure_values.json`, `fig2.probes`]

### 3.4 Pathway panel and scoring
- **Selection.**
  - Ten candidate MSigDB Hallmark sets (v2026.2.Hs, R7) were grouped by score correlation.
  - Six of the ten sat on one mesenchymal axis (r = 0.60–0.91 with EMT).
  - One representative per cluster was chosen from biology and literature, without outcome data.
  - The panel was frozen before any association analysis. [`docs/PHASE6C_PATHWAY_PANEL_REVIEW.md`]
- **Primary panel (four sets).**

  | Set | Member genes |
  |---|---:|
  | EMT | 198 |
  | E2F targets | 195 |
  | Wnt/β-catenin | 41 |
  | p53 pathway | 194 |
  | **Union** | **607** |

  The descriptive full-cohort score-matrix condition number is 2.95, against 77 for all ten candidates.
- **Sensitivity panel (five sets).** Adds PI3K/AKT/mTOR (condition number 6.95). It was fixed in advance as sensitivity-only.
- **Scoring.**
  - Each member gene is standardised with the mean and sample SD of the **training fold only**.
  - Each score is the mean of its member z-scores.
  - Every estimated constant is re-fitted inside each resample, so held-out patients never influence the features.
- **Biological rationale** (R1–R6):

  | Set | Axis |
  |---|---|
  | EMT | Mesenchymal/stromal state; CMS4 poor-prognosis subtype |
  | E2F targets | Proliferation |
  | Wnt/β-catenin | Initiating driver pathway |
  | p53 pathway | DNA-damage checkpoint |

  Each score summarises bulk mRNA. It is not a measurement of pathway activation, protein activity or mutation status.

### 3.5 Models
- **M1 (clinical):** cause-specific Cox model with age (restricted cubic spline), sex, stage III vs II, tumour location and adjuvant chemotherapy. Unpenalised.
- **M2 (primary):** M1 plus the four pathway scores.
  - A ridge penalty is applied to the pathway block only.
  - λ is tuned by five-fold event-stratified cross-validation inside each training set.
  - The full-data value is λ = 0.351.
- **Design choices.**
  - Leaving the clinical block unpenalised keeps the comparison scale-equivariant and makes M1 exactly nested in M2.
  - Shrinking the molecular block makes the gain in C a conservative measure of molecular contribution.
- **Sensitivity models.**
  - the five-set panel (λ = 0.152);
  - S-7, the representative-probe representation;
  - a batch-stratified fit;
  - leave-one-batch-out refits.

### 3.6 Validation and uncertainty
- **Discrimination.** Harrell's C. The estimand is the optimism-corrected difference ΔC = C(M2) − C(M1).
- **Point estimate.** θ̂ = +0.004316, the full-data optimism-corrected ΔC.
- **Nested bootstrap.**
  - For each of 750 event-stratified outer resamples, both models are refitted and optimism is re-estimated from 100 inner resamples.
  - The two models are paired on the same resamples.
  - Five runs differ only in random seed. They bound Monte Carlo error and are not independent experiments.
  - Each random draw has its own stream.
  - Runs were extended from 500 to 750 outer resamples under a pre-specified Monte Carlo precision rule (largest endpoint Monte Carlo SE ≤ 5 % of interval width). The final value is 4.2 %.
- **Verification.**
  - A frozen post-run verifier passed 60 of 60 checks before any result was reported.
  - A separate standard-library re-implementation recomputed 138 reported quantities from the recorded checkpoints, agreeing within 2.2×10⁻¹⁶. This re-computation did not refit models and is not an external replication.
- **Calibration.** Optimism-corrected calibration slope of the linear predictor.

### 3.7 External-validation feasibility
- Candidate GPL570 colorectal cohorts were screened for events, endpoint compatibility and patient overlap.
  - GSE14333 had 48 recurrence events and GSE33113 had 18. Both are exploratory only.
  - A precision-based assessment indicated roughly 338 events would be needed at an assumed C of 0.65. This is an indicative simulation approximation.
  - GSE17536 shares 129 unambiguous patient pairs with GSE14333 and was excluded.
- Confirmatory external validation was therefore not attempted. [`docs/PHASE4_PREACQUISITION.md`, `docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md`]

## 4. Results

### 4.1 Clinical model
- Stage III vs II: HR 1.97 (95 % CI 1.30–2.98).
- No other clinical term's interval excluded 1.
- The chemotherapy coefficient is confounded by indication and has no causal reading.

### 4.2 Discrimination

| Model | Apparent C | Optimism-corrected C |
|---|---:|---:|
| M1 clinical | 0.6242 | 0.6016 |
| M2 four-set (primary) | 0.6382 | 0.6059 |
| M2 five-set (sensitivity) | 0.6413 | 0.6066 |

- Optimism removed 0.023 (M1) and 0.032 (M2).
- The apparent gain of +0.0140 fell to a corrected **ΔC of +0.004316**.

**Corrected ΔC, 95 % percentile intervals (five seed-replicate runs):**

| Run | Interval |
|---|---|
| 0 | [−0.0086, +0.0519] |
| 1 | [−0.0090, +0.0469] |
| 2 | [−0.0093, +0.0484] |
| 3 | [−0.0091, +0.0518] |
| 4 | [−0.0091, +0.0494] |

- Every interval includes zero. Including zero does not establish that the true effect is zero.
- Basic and bias-corrected diagnostic intervals also include zero in all runs.
- The bootstrap distribution sits above θ̂ in every run, which is flagged and not explained.
- A different statistic, the paired out-of-original ΔC interval ([+0.0002, +0.0240] under the median representation; [−0.0001, +0.0224] under S-7), mixes in-sample and out-of-sample fits and is **not** an interval for the corrected ΔC.

[`reports/phase9/corrected_rerun_report/corrected_interval_results.json`; Figures 7 and 8]

### 4.3 Calibration
- Corrected calibration slopes: 0.833 (M1), 0.862 (M2), 0.836 (five-set).
- Slopes below 1 indicate residual overfitting.
- Predictions applied to new patients would need about 14–17 % shrinkage.

### 4.4 Pathway estimates
Hazard ratios are per one training-fold SD of score.

| Pathway | Primary HR (95 % Wald CI) | S-7 HR (95 % CI) | Batch-stratified HR (point) | Leave-one-batch-out log-HR range (sign changes) |
|---|---|---|---:|---|
| EMT | 1.070 (0.951–1.202) | 1.084 (0.937–1.254) | 1.073 | +0.009 to +0.204 (0) |
| E2F targets | 0.906 (0.807–1.016) | 0.872 (0.758–1.003) | 0.912 | −0.249 to −0.014 (0) |
| Wnt/β-catenin | 1.016 (0.906–1.138) | 1.011 (0.879–1.162) | 1.004 | −0.002 to +0.067 (1) |
| p53 pathway | 0.995 (0.882–1.123) | 0.993 (0.852–1.158) | 1.022 | −0.113 to +0.014 (2) |

- No primary-model interval excludes 1.
- Wald intervals are approximate under ridge penalisation.
- **S-7** re-summarises the same arrays, so it is a sensitivity analysis, not replication.
- **Five-set model:** E2F HR 0.862 (0.746–0.996). This interval excluding 1 coincided with collinearity after PI3K was added:
  - p53 ↔ PI3K r = 0.636;
  - model condition number 2.9 → 7.1;
  - SEs inflated by a factor of 1.20–1.34.

  It was not promoted.
- **Leave-one-batch-out:** omitted batches contain 14–103 patients, so the spread mixes batch identity with sample size. Apparent C across these refits ranged 0.614–0.659.

### 4.5 Subtype and batch structure (exploratory, outcome-free)
- **Dominant expression axis.** The cohort's dominant expression axis is MMR/CIMP/BRAF/tumour location, not stage.
- **MMR association.** All four scores were associated with MMR status (Kruskal–Wallis, n = 410, BH-adjusted):

  | Score | ε² | q |
  |---|---:|---:|
  | EMT | 0.0394 | 1.3×10⁻⁴ |
  | E2F | 0.0251 | 0.0021 |
  | Wnt | 0.0187 | 0.0070 |
  | p53 | 0.1716 | 4.7×10⁻¹⁶ |

  p53's is the largest among non-circular score–annotation tests.
- **Batch.**
  - Batch and extraction-method location effects on scores did not survive a permutation test stratified by the ComBat covariate `combat_mod` (0 of 10).
  - That covariate matches recorded MMR in all 410 samples, but its provenance is unresolved.
  - Four QC-flag dispersion effects survived.
  - Residual batch-associated structure remains (PC3, ε² 0.061).

[`reports/phase8c/score_vs_annotation.csv`; `docs/PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md`; `docs/PHASE8B_PCA_BATCH_INVESTIGATION.md`]

## 5. Discussion

**Main finding.**
- Most of the apparent molecular gain was optimism.
- After correction, the four pathway scores did not demonstrably improve recurrence discrimination beyond clinical variables.
- The interval is wide relative to the estimate because the cohort has 136 events. A small loss and a gain of about 0.05 in C are both compatible with the data, so the result is not evidence of equivalence.

**Biological interpretation.**
- **EMT and E2F directions.** EMT stayed above 1 and E2F below 1 in every reported fit, without any interval excluding 1.
  - EMT's direction is compatible with reports that stromal, mesenchymal expression marks poor-prognosis colorectal cancer (R1–R3). GSE39582 contributed to the CMS classification, so R1 is not independent of this cohort.
  - Evidence on proliferation and outcome is mixed (R4 vs R5).
- **Wnt and p53.** Both sat near zero and changed sign under leave-one-batch-out.
  - Their gene sets are also loosely co-expressed in this tissue: implied mean inter-gene correlation 0.041 and 0.033, against 0.304 for E2F and 0.241 for EMT.
  - That co-occurrence is a hypothesis, not a demonstrated cause. Near-zero coefficients change sign easily.
- **Bulk-tissue ambiguity.**
  - Bulk expression mixes tumour cells with stromal, immune and endothelial cells.
  - Cell composition was not measured, so an elevated EMT score cannot be attributed to tumour-cell EMT rather than stromal content.
  - The strong association of the p53 score with MMR status is an association with molecular subtype, not with recurrence, and establishes no mechanism.

**Methodological lessons.**
- Optimism correction changed the size of the answer about threefold.
- Defining the estimand in advance prevented a different, near-significant statistic from being reported as the result.
- A collinear sensitivity model produced the only "significant" pathway coefficient.

## 6. Limitations

1. **Internal validation only.** No external cohort with enough events was available, and the CIT partition is internal.
2. **Endpoint.** The recurrence-time unit is undocumented. Deaths without recurrence (46) are censored informatively. No competing-risk or absolute-risk analysis was possible.
3. **Bulk tissue.** Cell composition was not measured, and tumour-cell expression was not separated from other cells.
4. **Batch.** Coefficients are batch-sensitive. Batch effects were characterised, not resolved.
5. **Power.** 136 events limit precision. Calibration slopes below 1 show residual overfitting.
6. **Single-layer estimates.** Coefficients, calibration, S-7 and batch analyses have no nested intervals. Only ΔC was given a verified nested interval.
7. **Scope of verification.** The independent re-computation checked the reported statistics against recorded model outputs; it did not refit models.

## 7. Conclusion

The study did not demonstrate a reproducible improvement in recurrence discrimination from adding the four pathway scores. This does not establish that the biological pathways have no role in colorectal-cancer recurrence.

Three studies would advance the question:
- separating tumour-cell and stromal signal by deconvolution or single-cell data;
- testing pathway value within molecular subtypes;
- validating the frozen models in an independent cohort with about 340 or more recurrence events and a documented time unit.

## References

- R1. Guinney J, et al. The consensus molecular subtypes of colorectal cancer. *Nat Med* 2015;21:1350–1356. doi:10.1038/nm.3967
- R2. Calon A, et al. Stromal gene expression defines poor-prognosis subtypes in colorectal cancer. *Nat Genet* 2015;47:320–329. doi:10.1038/ng.3225
- R3. Isella C, et al. Stromal contribution to the colorectal cancer transcriptome. *Nat Genet* 2015;47:312–319. doi:10.1038/ng.3224
- R4. Anjomshoaa A, et al. Slow proliferation as a biological feature of colorectal cancer metastasis. *Br J Cancer* 2009;101:822–828. doi:10.1038/sj.bjc.6605229
- R5. Luo ZW, et al. Increased expression of Ki-67 is a poor prognostic marker for colorectal cancer patients: a meta analysis. *BMC Cancer* 2019;19:123. doi:10.1186/s12885-019-5324-y
- R6. Cancer Genome Atlas Network. Comprehensive molecular characterization of human colon and rectal cancer. *Nature* 2012;487:330–337. doi:10.1038/nature11252
- R7. Liberzon A, et al. The Molecular Signatures Database hallmark gene set collection. *Cell Syst* 2015;1:417–425. doi:10.1016/j.cels.2015.12.004
- R8. Marisa L, et al. Gene expression classification of colon cancer into molecular subtypes. *PLoS Med* 2013;10:e1001453. doi:10.1371/journal.pmed.1001453

References were verified against published records (`docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md`).
