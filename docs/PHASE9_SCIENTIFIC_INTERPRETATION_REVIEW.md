# Scientific interpretation and discussion audit — GSE39582 recurrence-risk project

**Date:** 8 November 2026. **Type:** review only. No analysis was run, and no report or frozen artifact was edited. Proposed corrections to the final report are listed in section 10 and await owner approval.

**Source abbreviations used throughout:**

| Code | File |
|---|---|
| FR | `docs/PHASE9_CORRECTED_FINAL_REPORT.md` |
| CIR | `reports/phase9/corrected_rerun_report/corrected_interval_results.json` |
| IA | `reports/phase9/independent_audit/INDEPENDENT_AUDIT_REPORT.md` |
| VR | `reports/phase9/corrected_rerun/verification/L00750.report.txt` |
| EX | `reports/phase9/execution/phase9_execution_results.json` (keys given) |
| ER | `docs/PHASE9_EXECUTION_REPORT.md` |
| S7 | `docs/PHASE9_S7_SENSITIVITY_REPORT.md` |
| TA | `docs/PHASE9_THETA_ALIGNMENT_AUDIT.md` |
| P6A | `docs/PHASE6A_ENDPOINT_RESOLUTION.md` |
| P5 | `docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md` |
| P6C | `docs/PHASE6C_PATHWAY_PANEL_REVIEW.md` |
| P7C | `docs/PHASE7C_QC_REPORT.md` |
| P8B | `docs/PHASE8B_PCA_BATCH_INVESTIGATION.md` |
| PLAN | `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` v2.1 |
| PROT | `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md` v3.1 |

---

## 1. Results interpretation

### 1.1 Cohort and endpoint
- **Measured:**
  - 452 patients (A-19 complete-case subset of the 453-tumour S3 cohort), 136 recurrences, 316 censored (EX `population`; FR §2).
  - QC excluded no sample (P7C §0 D2).
  - Deaths without recurrence: 46 of 453, censored in the deposited data (P6A §2.2).
- **Estimated:** nothing. These are counts.
- **Supports:** an internally consistent cause-specific recurrence analysis with about 13.6 events per M2 coefficient (136/10). PROT §2.1 gives EPV 15.0–15.1, computed on its parameter budget.
- **Does not establish:**
  - absolute risk at any time point, because the time unit is undocumented (U1; P6A §1, row U1 in the summary table);
  - any statement about cumulative incidence of recurrence (U2; P6A §3.2).

### 1.2 Clinical baseline model (M1)
- **Measured:** apparent C = 0.6242 (EX `apparent.M1.C`).
- **Estimated:**
  - optimism-corrected C = 0.6016, MC-SE 0.00058 (EX `monte_carlo.M1`; FR §7.3);
  - stage III vs II HR 1.97 (1.30–2.98) (EX `apparent.M1.hazard_ratios`; FR §9);
  - male sex HR 1.38 (0.97–1.96);
  - other clinical HRs with intervals including 1.
- **Supports:** modest clinical discrimination, driven mainly by stage.
- **Does not establish:**
  - adequate discrimination for individual decisions;
  - stability of the stage effect under batch stratification. In the M2 fit, the stage coefficient falls from 0.662 (HR 1.94) to 0.520, HR 1.68 (1.07–2.65) (ER §10; EX `batch_sensitivity.stratified_fit`).

### 1.3 Four-pathway primary model (M2)
- **Measured:** apparent C = 0.6382 (EX `apparent.M2_primary.C`).
- **Estimated:**
  - corrected C = 0.6059, MC-SE 0.00047 (FR §7.3);
  - λ = 0.351 (EX).
- **Supports:** a small apparent gain of +0.0140, most of which is optimism, because M2's optimism (about 0.032) exceeds M1's (about 0.023) (FR §7.3).
- **Does not establish:** that the molecular block adds discrimination (section 1.7).

### 1.4 Five-pathway sensitivity model
- **Estimated:**
  - corrected C = 0.6066 (EX `monte_carlo.M2_sensitivity`; FR §8.1, which records it as historical);
  - corrected slope 0.836.
  - **E2F HR 0.862 (0.746–0.996)** is the only pathway interval in any model that excludes 1, and only just. PI3K HR 1.043 (0.891–1.221) (ER §5).
  - Adding PI3K raised collinearity: P53↔PI3K r = 0.636; condition number 2.9 → 7.1; molecular SEs inflated ×1.20–1.34; the E2F estimate moved −0.099 → −0.149 (ER §5).
- **Supports:** no difference from the primary model. The C difference of 0.0007 is within Monte Carlo error (ER §15).
- **Does not establish:**
  - a five-set advantage;
  - a robust E2F association. The borderline interval coincides with collinearity-driven coefficient inflation, in a sensitivity model with no interval for ΔC.

### 1.5 Individual pathway coefficients (four-set, per training-fold SD)
| Pathway | HR (95 % Wald CI), median | S-7 | LOBO coefficient range (sign changes) | Batch-stratified HR |
|---|---|---|---|---|
| E2F targets | 0.906 (0.807–1.016) | 0.872 (0.758–1.003) | −0.249 to −0.014 (0) | 0.912 (0.811–1.025) |
| EMT | 1.070 (0.951–1.202) | 1.084 (0.937–1.254) | +0.009 to +0.204 (0) | 1.073 (0.952–1.208) |
| Wnt/β-catenin | 1.016 (0.906–1.138) | 1.011 (0.879–1.162) | −0.002 to +0.067 (1) | 1.004 (0.895–1.127) |
| p53 | 0.995 (0.882–1.123) | 0.993 (0.852–1.158) | −0.113 to +0.014 (2) | 1.022 (0.904–1.155) |

Sources: EX `apparent.M2_primary`; S7 §3 and §7; EX `batch_sensitivity.leave_one_batch_out` and `.stratified_fit`; ER §10; FR §9.

- **Supports:** no individual pathway association distinguishable from HR = 1 in the primary model.
- **Direction:** E2F < 1 and EMT > 1, consistently across representation, batch stratification and every leave-one-batch-out refit.
- **Does not establish:** independent pathway effects. The intervals come from a ridge-penalised fit and are approximate (FR §9). The Wnt and p53 coefficients change sign when one batch is omitted.

### 1.6 Apparent vs optimism-corrected discrimination
| | Apparent C | Corrected C | Optimism |
|---|---:|---:|---:|
| M1 | 0.6242 | 0.6016 | ≈ 0.023 |
| M2 | 0.6382 | 0.6059 | ≈ 0.032 |

Sources: EX; FR §7.3.

- **Supports:** the apparent advantage (+0.0140) is mostly optimism.
- **Does not establish:** performance outside this population.

### 1.7 Corrected ΔC and the five run-specific intervals
- **Estimated:** θ̂ = +0.004316. Its alignment with the interval is shown in TA, and its value is reproduced in IA §3.
- **Five 95 % percentile intervals** (CIR `runs.*.percentile_primary`; FR §8; IA §3):

  | Run | Interval |
  |---|---|
  | 0 | [−0.0086, +0.0519] |
  | 1 | [−0.0090, +0.0469] |
  | 2 | [−0.0093, +0.0484] |
  | 3 | [−0.0091, +0.0518] |
  | 4 | [−0.0091, +0.0494] |

- **Run-to-run variation:** between-run SD 0.00025 (lower endpoint) and 0.00217 (upper).
- **Diagnostics:** EP-1 met; no hard failure. A7/A8 are flagged in every run (z₀ −0.43 to −0.59). All 15 percentile, basic and bias-corrected intervals include zero (FR §8).
- **Supports:** the data are compatible with roughly −0.009 to +0.05 in C, and the Monte Carlo precision is adequate.
- **Does not establish:**
  - a positive gain;
  - equivalence (FR §14).
  - Because the bootstrap-bias warning is unexplained, the interval's location relative to the true ΔC is less certain than its width suggests (FR §8 "What the warnings show").

### 1.8 Calibration slopes
- **Estimated:** corrected log-hazard slopes of 0.833 (M1), 0.862 (M2) and 0.836 (five-set) (EX `monte_carlo.*.per_seed.corrected_slope`; FR §10). Under S-7, M2's slope is 0.832 (S7 §4).
- **Supports:** residual overfitting in every model under both representations.
- **Does not establish:**
  - calibration-in-the-large;
  - absolute-risk calibration (U1).
  - These values are historical, single-layer estimates without intervals.

### 1.9 S-7 expression-representation sensitivity
- **Estimated:** corrected C and ΔC are identical to four decimals under the median and representative-IQR representations (+0.0043) (S7 §4). The four-set coefficients do not change sign (S7 §7).
- **Paired out-of-original ΔC** (a different quantity; section 3.3):
  - median: +0.0121 [+0.0002, +0.0240];
  - S-7: +0.0105 [−0.0001, +0.0224].
- **Five-set model:** under S-7 its corrected C moves beyond Monte Carlo error (−0.0025), with worse P53–PI3K collinearity (S7 §4 and §8).
- **Supports:** the four-set result is stable to how probes are summarised.
- **Does not establish:** replication. It uses the same patients, arrays and genes (S7 §10).

### 1.10 Batch-related instability
- **Measured:**
  - nine ComBat batches, confounded with the data partition (Cramér's V = 0.488; P5 §4.2);
  - residual batch association with PC3 (ε² = 0.061; P8B §4);
  - batches range 14–103 patients and 4–52 events (ER §10).
- **Estimated:**
  - leave-one-batch-out apparent C ranges 0.614–0.659;
  - molecular coefficients vary by up to about an order of magnitude, with sign changes for Wnt and p53 (ER §10);
  - the batch-stratified fit leaves the molecular coefficients nearly unchanged (Δ ≤ 0.027) but attenuates stage (ER §10).
- **Supports:** the molecular coefficient estimates are fragile to which batches are included.
- **Does not establish:**
  - a batch effect on outcome;
  - which batch is "right".
  - The leave-one-batch-out spread conflates batch identity with sample size, since omitting batch03 or batch04 removes about 23 % of the cohort (PROT, Diagnostic 1 caveat; ER §10).

### 1.11 Independent numerical reproducibility
- **Measured:**
  - 3,750 outer replicates and 375,000 inner pairs;
  - all record hashes valid;
  - 138 quantities match the official results within 2.2×10⁻¹⁶ (IA §2–4);
  - frozen verifier PASS, 60/60, including a bit-identical recomputation of 20 replicates (VR).
- **Supports:** the reported statistics follow exactly from the recorded fits, and the recorded fits are reproducible.
- **Does not establish:** scientific validity, external validity or the absence of shared systematic error (IA §5; PLAN §2.3).

---

## 2. Biological interpretation

Background biology below is general knowledge from the colorectal-cancer literature. It describes **plausibility**, not evidence from this study. A Hallmark score here is the mean standardised **mRNA** abundance of a fixed gene set (FR §4; P6C). It is not a measure of pathway activation, protein activity, signalling flux or mutation status.

| Pathway | Biological background (plausibility) | Observed here (evidence) | What can be said |
|---|---|---|---|
| **EMT** | Mesenchymal and stromal programmes characterise the poor-prognosis CMS4 subtype (Guinney et al., *Nat Med* 2015). Bulk-tumour EMT expression often reflects stromal and fibroblast content as much as tumour-cell EMT | HR 1.07 (0.95–1.20); direction > 1 in every leave-one-batch-out refit, under S-7 and under stratification | The direction is consistent with the hypothesis that mesenchymal or stromal-rich tumours recur more. The effect is not distinguishable from zero, and the score cannot separate stroma from tumour-cell EMT |
| **E2F targets** | Cell-cycle and proliferation genes. In colorectal cancer, higher proliferation signatures have been associated with *better* outcome in several reports, unlike many other cancers | HR 0.91 (0.81–1.02), < 1 in every refit; borderline 0.86 (0.75–0.996) only in the collinear five-set model | The direction matches the reported "proliferation paradox" in colorectal cancer. Unproven here: the interval includes 1 in the primary model, and the one excluding interval is in a collinearity-inflated sensitivity model |
| **Wnt/β-catenin** | Wnt pathway activation (usually through APC loss) is nearly universal in colorectal cancer, which limits between-tumour contrast. Expression of target genes is an indirect proxy | HR 1.02 (0.91–1.14); sign changes once under leave-one-batch-out; lowest cross-representation score correlation (r = 0.857; S7 §7) | No signal. The near-ubiquity of activation and the representation sensitivity make a null expected, not informative |
| **p53 pathway** | The Hallmark set captures p53 **transcriptional targets**. TP53 mutation (common in colorectal cancer) is not equivalent to low target expression | HR 1.00 (0.88–1.12); sign changes twice under leave-one-batch-out | No signal. A target-gene score is not a TP53 mutation proxy, so nothing can be inferred about TP53 status |
| **PI3K/AKT/mTOR** (sensitivity only) | PIK3CA mutations occur in a minority of colorectal cancers. Expression of pathway genes does not equal pathway activity | HR 1.04 (0.89–1.22); collinear with p53 (r = 0.64; ER §5) | No signal. Its addition mainly destabilised other coefficients, as anticipated (P6C; ER §15) |

**Overall biological reading.** The two consistently signed directions, E2F < 1 and EMT > 1, are **hypothesis-generating** and biologically coherent. The study demonstrates neither an association nor a mechanism. Testing them would need independent cohorts and ideally orthogonal measures: stromal content, mutation data, protein-level assays.

---

## 3. Statistical interpretation

1. **Why +0.0043 is marginal.**
   - It is less than 1 % of the C scale (about 0.60).
   - It is small relative to the between-model differences that typically matter for prognostic models, and smaller than the leave-one-batch-out spread of apparent C (0.045).
   - It is about one-third of the apparent gain, because the extra four coefficients carry more optimism (FR §7.3).
2. **Why intervals including zero prevent an improvement claim.**
   - All five percentile intervals span zero, and their lower limits (about −0.009) admit a small loss (FR §8).
   - With 136 events the data cannot separate a null gain from a gain of about 0.05.
   - Absence of a demonstrated gain is not evidence of no gain (FR §14).
3. **Why paired out-of-original intervals are not corrected-ΔC intervals.**
   - The paired quantity, +0.0121 [+0.0002, +0.0240], is ΔC of bootstrap-fitted models evaluated on the original data (ER §7). It mixes in-sample and out-of-sample information and is **not** optimism-corrected.
   - Its interval reflects variation in the bootstrap fits, not sampling uncertainty of ΔC_corrected.
   - It narrowly excluded zero under the median representation and included zero under S-7 (S7 §4). That is a knife-edge, and it must not be transferred to the corrected ΔC.
   - **The final report does not mention this quantity at all** (gap G4).
4. **Why coefficient intervals including HR = 1 matter.** No single pathway has an adjusted association distinguishable from none. So even a positive ΔC could not be attributed to any specific pathway. Under ridge shrinkage the Wald intervals are approximate (FR §9).
5. **Why slopes below 1 matter.**
   - Slopes of 0.83–0.86 mean the linear predictors are too spread out for new patients and need about 14–17 % shrinkage (FR §10).
   - Overfitting affects even the six-coefficient clinical model, which underlines the limited information in 136 events.
6. **Why the five-set model must not be promoted.**
   - It was pre-specified as sensitivity-only, outcome-blind, before any fit (P6C; PROT D8D-2).
   - Its C advantage (0.0007) is within Monte Carlo error.
   - It introduces collinearity that inflates SEs and moves E2F.
   - Under S-7 its C changes beyond Monte Carlo error.
   - Promoting it on performance would be outcome-driven model selection (ER §5, §15; S7 §8).
7. **Why S-7 is sensitivity, not replication.** It uses the same patients, arrays and normalisation, with only the probe summarisation differing. Agreement shows pipeline robustness, not an independent confirmation (S7 §10).
8. **Why batch instability limits coefficient interpretation.**
   - When omitting single batches flips the sign of two coefficients and moves the other two several-fold (ER §10), the coefficient values depend on cohort composition.
   - With batch confounded with partition and residual batch structure present (P5; P8B), part of a pathway score's variation may be technical. Coefficient-level interpretation is therefore not credible.
9. **Penalisation asymmetry** (ER §15). Only the molecular block is penalised, so the M2-vs-M1 comparison is conservative by design. That cannot convert a null-consistent result into a positive one, but it should be stated (gap G5).
10. **Limits of the C-index.** C measures ranking only:
    - it is insensitive to added predictors when the baseline is moderate;
    - it says nothing about calibration or decision consequences;
    - net benefit, reclassification and time-dependent measures were not assessed, and time-dependent measures cannot be until U1 is resolved.

---

## 4. Robustness and sensitivity interpretation

| Dimension | Result | Source |
|---|---|---|
| Monte Carlo (five runs) | Robust: endpoint SD 0.00025 / 0.00217; EP-1 ratio ≤ 0.042 | CIR; FR §8 |
| Seed-correction | Robust: corrected and as-run pooled ranges nearly identical | FR §8.2 (as historical), §14 |
| Expression representation | Robust for the four-set C and ΔC; not for the five-set model | S7 §4, §10 |
| Panel choice (four- vs five-set) | No difference in C; five-set coefficients destabilised | ER §5, §15 |
| Batch composition | **Not robust** for coefficients or apparent C | ER §10 |
| Bootstrap bias diagnostics | Warning in all runs; direction unexplained | CIR `warnings`; FR §8 |
| Numerical reproducibility | Independently reproduced | IA; VR |

Robustness here is mostly **agreement about a null-consistent result** (S7 §10). It shows the pipeline is stable, not that an effect is real.

---

## 5. Limitations

1. **Time unit (U1):** undocumented. No absolute risk, t-year calibration or time-dependent accuracy is possible (P6A §1; PROT, U1 row).
2. **Death without recurrence (U2):**
   - The deposited data censor 46 deaths (10.2 %) in S3, and that censoring is informative: those patients are about 7 years older, less often stage III and about half as likely to have had chemotherapy (P6A §2.3).
   - The cause-specific hazard remains estimable, but its meaning is conditional, and it does not give cumulative incidence (P6A §3.2).
   - A Fine–Gray or composite (recurrence-or-death) endpoint would answer different questions and was not analysed.
3. **Events:** 136 events, about 13.6 per M2 coefficient, which limits precision. Overfitting is visible in every model (FR §10).
4. **Validation:**
   - Internal only.
   - The CIT partition is confounded with batch and has 29 events, so it is at most an internal holdout (P5 §4.2).
   - No external cohort met precision requirements: GSE14333 has 48 events, GSE33113 has 18 (P5).
5. **Batch:** partition confounding; residual batch structure; leave-one-batch-out coefficient instability (section 1.10).
6. **Overfitting:** corrected slopes of 0.83–0.86 under both representations.
7. **C-index:** a ranking measure only, with no clinical-utility assessment.
8. **Partly historical:** θ̂, calibration, S-7, the five-set model, coefficients and batch analyses come from the single-layer Phase 9 primary execution, without intervals and without corrected streams. Only ΔC has a verified nested interval.
9. **Unexplained bootstrap shift:** the outer distribution sits above θ̂, opposite to PLAN §7.2's anticipated direction (FR §8).

---

## 6. Clinical implications

- **None can be drawn for patient care.** The study does not show:
  - that the four-pathway scores improve recurrence-risk discrimination;
  - that they reclassify patients usefully;
  - that they could guide adjuvant therapy.
- A C of about 0.60 for both models is modest, so neither model is suitable for individual risk communication.
- No absolute recurrence probabilities can be produced (U1/U2).
- Generalisation is limited to GPL570 data processed like GSE39582; no transportability is shown.

---

## 7. Overall discussion

The project asked a narrow, pre-specified question and answered it with unusually strong computational controls: frozen code, verified and independently reproduced nested-bootstrap intervals. The answer is **inconclusive-to-negative**. The optimism-corrected gain is small (+0.0043), and its intervals span zero, from a small loss to a moderate gain.

Supporting analyses reinforce that reading:
- no individual pathway coefficient is distinguishable from no effect;
- both models are overfitted;
- coefficients are unstable across batches;
- the only "significant-looking" results do not survive scrutiny:
  - the paired out-of-original interval is a different quantity and fails under S-7;
  - the five-set E2F interval sits in a collinear sensitivity model.

The consistent directions for E2F (protective) and EMT (adverse) are biologically coherent hypotheses, not findings. The main methodological lesson is how much of an apparent omics gain (+0.014) is optimism (about two-thirds here), and that computational stability is not statistical stability (ER §15).

---

## 8. Final conclusion

**Does adding reproducible transcriptomic pathway information improve recurrence-risk prediction beyond clinical characteristics?**

**Not demonstrably, in this cohort.**
- The pre-specified optimism-corrected improvement in Harrell's C is +0.0043.
- Its 95 % percentile intervals, in five independent and independently reproduced runs, span about −0.009 to +0.05.
- Neither a small loss nor a modest gain can be excluded.
- No pathway shows an independent association, and the molecular coefficients are unstable across processing batches.

The result is **inconclusive with respect to a modest benefit and offers no support for a robust one**. It does not establish that the pathways are uninformative. Resolving this would need larger, externally independent cohorts with documented endpoints.

---

## 9. Research question → evidence → conclusion

| Question | Evidence (source) | Conclusion |
|---|---|---|
| Do pathways improve corrected discrimination? | θ̂ +0.0043; five intervals about [−0.009, +0.05] (CIR; FR §8; IA §3) | Not demonstrated; inconclusive for a modest gain |
| Is the estimate numerically trustworthy? | VR PASS 60/60; IA reproduced within 2.2×10⁻¹⁶; EP-1 met (CIR) | Yes, as computation |
| Is the interval location reliable? | A7/A8 in all runs, unexplained direction (FR §8) | Warning; qualitative conclusion holds for all 15 intervals |
| Does any pathway have an independent association? | All four-set HR CIs include 1 (EX; FR §9) | No |
| Are coefficient directions stable? | E2F < 1 and EMT > 1 across S-7, stratified and leave-one-batch-out fits; Wnt and p53 flip (ER §10; S7 §7) | Directions stable for E2F and EMT only; magnitudes unstable |
| Are models well calibrated? | Corrected slopes 0.83–0.86 (EX; FR §10; S7 §4) | Residual overfitting |
| Does the five-set panel help? | ΔC_corrected vs four-set 0.0007; collinearity; S-7 shift (ER §5, §15; S7 §8) | No; remains sensitivity-only |
| Is the result robust to representation? | Corrected ΔC +0.0043 under both (S7 §4) | Yes, for the four-set model; not replication |
| Is the result robust to batch? | Leave-one-batch-out C 0.614–0.659; coefficient sign flips (ER §10) | No, for coefficients and apparent C |
| Is there external validity? | None performed; no adequate cohort (P5) | Not established |
| Is there clinical utility or absolute risk? | U1/U2; C only (P6A) | Cannot be assessed or claimed |

---

## 10. Gaps in the current final report (proposed corrections; not yet applied)

| # | Gap | Where | Proposed correction |
|---|---|---|---|
| G1 | **Leave-one-batch-out coefficient instability is not reported**: only the apparent-C range is. EMT ranges +0.009 to +0.204, E2F −0.249 to −0.014, and Wnt and p53 change sign (ER §10) | FR §12 | Add a coefficient-range table and the sample-size caveat |
| G2 | **Batch-stratified fit results are not reported**: molecular Δ ≤ 0.027, stage 0.662 → 0.520 (ER §10) | FR §12 | Add one sentence |
| G3 | **Five-set E2F interval excludes 1 (0.746–0.996)** but is not mentioned. The README and FR state "every pathway CI includes 1" without restricting it to the primary model | FR §8.1, §9; README | Restrict the statement to the primary model; report the five-set E2F result with the collinearity explanation |
| G4 | **Paired out-of-original ΔC** (+0.0121 [+0.0002, +0.0240]; S-7 +0.0105 [−0.0001, +0.0224]) is omitted, so the report does not explain why it is not a corrected-ΔC interval | FR §7–8 | Add a short "distinct quantities" paragraph (ER §7; S7 §4) |
| G5 | **Penalisation asymmetry** (molecular penalised, clinical not; conservative comparison) is not stated | FR §5, §14 | Add one sentence (ER §15) |
| G6 | **S-7 five-set shift** beyond Monte Carlo error (−0.0025) and **P53–PI3K collinearity** are not reported | FR §11 | Add a row and a sentence (S7 §4, §8) |
| G7 | **Informative-censoring profile** (older, less stage III, less chemotherapy) and the Fine–Gray / composite distinction are missing | FR §2, §13 | Expand the U2 text (P6A §2.3, §3.2) |
| G8 | **C-index limitations** for clinical usefulness (no decision-curve, reclassification or time-dependent analysis) are not explicit | FR §13–14 | Add a limitation item |
| G9 | **Biological section** cites "some published reports" without references and lacks the plausibility-vs-evidence framing and per-pathway caveats | FR §9 | Replace with section 2 of this review, in condensed form |
| G10 | **Section 12 compares an apparent-C range with a ΔC interval distance**, which is not like-for-like | FR §12 | Remove the comparison and state the C range on its own |
| G11 | **ER §15 "not ready" blockers are not closed explicitly.** S-7 was executed; the corrected-ΔC interval now exists; the batch instability remains unexplained | FR §14 or §15 | Add a three-line status note |
| G12 | **Events per coefficient** (about 13.6 for M2) is not stated | FR §2 / §13 | Add to the events limitation |
| G13 | **Stage attenuation under stratification** (M2 stage HR 1.94 → 1.68) qualifies the clinical baseline | FR §9 | Add a note |

No gap changes a reported number or the conclusion. All of them add interpretation or restrict wording. **Approved by the owner and applied on 8 November 2026**; see `docs/release/G1_G13_IMPLEMENTATION.md`.
