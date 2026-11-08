# Phase 5 — Precision-Based Feasibility Assessment (rev. 2, audited)

**Date:** 2026-11-02 · **Revision 2** — methodological audit applied; a claim withdrawn.
**Status:** Assessment complete. **No model built, no expression data touched, no packages installed, Phase 3 bundle unmodified.**

**Reproduce:**
```
.venv/bin/python scripts/25_validation_precision_simulation.py --reps 400
.venv/bin/python scripts/26_phase5_methodological_audit.py --reps 30
```
→ `reports/phase5/`

> **All simulation data is synthetic.** The linear predictor is drawn from an assumed distribution; it is **not** derived from GSE39582 expression. No predictive model was fitted to any real data.

---

## 0. METHODOLOGICAL AUDIT — what this implementation is and is not

### 0.1 Correction and withdrawal

> **WITHDRAWN.** Revision 1 stated that the ~338-event figure *"closely matches Riley's published ~365 — an independent cross-check."* **That claim is withdrawn.** It is not a cross-check. The two numbers come from different estimators, different clinical settings, different linear-predictor distributions, different censoring mechanisms and different absolute risks. Their proximity is **coincidental and carries no validation weight.**

> **This simulation is an APPROXIMATION of the Riley framework, not a faithful implementation of it.** It is labelled as such throughout.

### 0.2 Step-by-step comparison with Riley et al. 2022

| # | Riley 2022 (Step 1 set-up, §4.1) | This implementation | Match |
|---|---|---|---|
| 1 | **Specify a time point of interest.** VTE: 3 years | **No time point.** Calibration assessed over the whole follow-up | ❌ **DIVERGES** |
| 2 | **LP distribution in the validation population.** VTE: **skew-normal** (mean 4.60, var 0.65, skewness 0.5, kurtosis 5), taken from the development article's histogram | **Normal(0, σ²)**, σ back-solved from an *assumed* C-index | ⚠️ **approximation** |
| 3 | **Specify F(t)** — absolute risk at the horizon. VTE: F(3) = 0.17 | **No F(t).** An "event fraction" over the whole follow-up is used instead | ❌ **DIVERGES** |
| 4 | **Survival distribution conditional on LP.** VTE: exponential, baseline rate ≈ 0.00050 chosen so F(3)=0.17 when the log HR for LP is 1 | exponential, rate exp(LP); baseline folded into the censoring quantile | ⚠️ **similar in spirit** |
| 5 | **Censoring distribution + max follow-up.** VTE: **exponential censoring**, ~72% censored before 3 years | **administrative censoring only** | ❌ **DIVERGES** — random censoring removes more information, so this is **optimistic** |
| 6 | **Target SE for the calibration slope** (0.051 / 0.102) | same targets | ✅ |
| 7 | Steps 2–9: iterate n to hit the target SE | grid over event counts; required n obtained by fitting SE ≈ k/√E | ⚠️ equivalent in effect |
| 8 | **Calibration slope from Equation (3):** pseudo-observations F_i(t), **complementary log-log link**, `ln(−ln(1−E[F_i(t)])) = α + β·ln(−ln(1−F_i(t)))`, fitted by **GEE with a sandwich variance estimator** (needed because pseudo-observations are correlated) | **Cox regression of outcome on the linear predictor**; β is the slope; SE from the **model-based information matrix** | ❌ **DIFFERENT ESTIMATOR, DIFFERENT SCALE, DIFFERENT VARIANCE ESTIMATOR** |
| 9 | Also report **mean within-dataset CI width** for the C-index | **1.96 × empirical SD across replications** (Monte-Carlo sampling SD) | ❌ **different quantity** (§0.5) |

**Verbatim, Riley's Equation (3):**

> "Royston suggests using a complementary log-log link function, such that `ln(−ln(1−E[F_i(t)])) = α + β ln(−ln(1−F_i(t)))` … where β is the calibration slope, and F_i(t) is the predicted risk calculated from the existing prediction model for individual i at time t."

> "Model fitting is discussed by Andersen and Perme, who suggest the use of generalized estimating equations followed by a sandwich estimator for the variance of parameter estimates (which is needed to account for the correlation of the pseudo-observations themselves)."

### 0.3 What the implementation actually estimates

A **Cox-on-linear-predictor calibration slope** — a standard, widely used calibration measure, but **not** Riley's pseudo-observation cloglog slope. It:

- assesses calibration **on the log-hazard scale over the whole follow-up**, not absolute risk at a horizon;
- uses the **model-based** information matrix, not a sandwich estimator;
- requires **no predicted risks**, which is why it could be run at all (§0.4).

**Consequence for the conclusions.** The two estimators are not interchangeable, so the *absolute* event requirements in §5.2 should be read as **indicative, not as Riley-derived numbers**. The *ordering and magnitude* of the shortfall — 48 and 18 events against requirements in the hundreds — is driven by the √E scaling of any calibration-slope SE and is not an artifact of the estimator choice. **The feasibility conclusion survives; the precision of the numbers behind it does not.**

### 0.4 Inputs unavailable because no model exists yet

Riley's framework requires inputs that **cannot be supplied until a model and a horizon are specified**:

| Required input | Status | Why unavailable |
|---|---|---|
| **Prediction horizon t** | **UNKNOWN** | No horizon has been chosen. Compounded by the unresolved time unit (U1) — a "3-year" horizon cannot even be expressed in `rfs.delay` units |
| **Predicted risks F_i(t)** | **UNAVAILABLE** | Requires a fitted model with a baseline hazard. None exists |
| **LP distribution (shape, not just spread)** | **ASSUMED** | Riley obtains it from the development article's histogram. Ours is assumed normal |
| **F(t), absolute risk at the horizon** | **UNKNOWN** | Requires the horizon and a baseline survival estimate |
| **C-index or Royston's D of the actual model** | **UNKNOWN** | Assumed 0.60–0.70 from literature, not measured |
| **Censoring distribution in the validation population** | **UNKNOWN** | Not characterised for GSE14333 or GSE33113 |

**A faithful Riley calculation is therefore not currently possible.** It becomes possible once a model, a horizon and a time unit exist — and should be redone then.

### 0.5 C-index confidence-interval audit

| Item | This implementation |
|---|---|
| **Estimator** | **Harrell's C** — concordant/comparable pairs, ties credited 0.5 |
| **Censoring assumption** | pairs comparable only when the earlier time is an observed event. **Harrell's C is known to depend on the censoring distribution.** The censoring-robust alternative (**Uno's IPCW C**) is **not implemented** |
| **Replicates** | 400 (script 25); 300 (audit) |
| **CI method** | **1.96 × empirical SD of C across replications** |
| **What that is** | a **Monte-Carlo estimate of the sampling SD of the estimator** |
| **What Riley reports** | the **mean within-dataset 95% CI width** |
| **Equivalent?** | **No** — related, not the same |

**Empirical check** (`reports/phase5/c_index_ci_audit.csv`): comparing the Monte-Carlo width against a within-dataset bootstrap CI (300 resamples):

| Events | Monte-Carlo width | Bootstrap within-dataset width | ratio |
|---:|---:|---:|---:|
| 48 | 0.1512 | 0.1389 | 0.92 |
| 136 | 0.0945 | 0.0941 | **1.00** |
| 365 | 0.0566 | 0.0547 | 0.97 |

The approximation **tracks a within-dataset bootstrap CI closely** (ratio 0.92–1.00). That is reassuring for the §5.4 figures, but it is a documented approximation, not Riley's measure, and it inherits Harrell's C's censoring dependence.

---

## 1. Headline

> **No available cohort reaches even the relaxed precision target for a calibration slope.**
> **Confirmatory external validation is not feasible with the data identified.**

| Cohort | Events | Slope 95% CI width (C=0.65) | Riley relaxed target (0.4) |
|---|---:|---:|---|
| GSE39582 S3, full follow-up | 136 | **0.63** | ❌ not met |
| GSE39582 S3, 60-unit (conditional) | 130 | **0.64** | ❌ not met |
| **GSE14333** (external) | **48** | **1.07** | ❌ not met |
| **GSE33113** (external, stage II, composite) | **18** | **1.79** | ❌ not met |

---

## 2. Method

### 2.1 Framework

Riley RD, Debray TPA, Collins GS, Archer L, Ensor J, van Smeden M, Snell KIE.
*Minimum sample size calculations for external validation of a clinical prediction model with a time-to-event outcome.*
**Statistics in Medicine** 2022;41(7):1280–1295. DOI [10.1002/sim.9275](https://doi.org/10.1002/sim.9275) · PMID 34915593

Quoted from that paper:

> "A simulation-based framework is proposed, which calculates the sample size required to target a particular confidence interval width for the calibration slope… The process requires assumptions about the validation population in terms of the (i) distribution of the model's linear predictor and (ii) event and censoring distributions."

> "We focus predominantly on the sample size required to precisely estimate the calibration slope as our previous articles, and those of others, have shown… that the calibration slope is usually the measure that requires the largest sample size to estimate precisely."

**Their stated targets:**

| Target | SE of calibration slope | CI width | Slope interval |
|---|---:|---:|---|
| **Primary** | **0.051** | **0.2** | 0.9 – 1.1 |
| **Relaxed** | **0.102** | **0.4** | 0.8 – 1.2 |

In their worked VTE example, the primary target required **~14,250 individuals (~1,430 events)** and the relaxed target **~3,600 (~365 events)**. They note plainly: *"Such numbers may not be possible in practice."*

**Supporting references:**
- Collins GS, Ogundimu EO, Altman DG. *Sample size considerations for the external validation of a multivariable prognostic model: a resampling study.* **Stat Med** 2016;35:214–226. DOI [10.1002/sim.6787](https://doi.org/10.1002/sim.6787)
- Snell KIE et al. *External validation of clinical prediction models: simulation-based sample size calculations were more reliable than rules-of-thumb.* **J Clin Epidemiol** 2021. PMC8352630

**These are precision criteria, not significance thresholds**, and **100/200 events are not used here as universal cutoffs.** The screening bands adopted in Phase 4 (D3) were explicitly labelled screening thresholds; this assessment replaces them with the simulation they called for.

### 2.2 Implementation

`src/precision_sim.py` — synthetic cohort generator, univariate Cox partial-likelihood fit (Breslow ties, Newton–Raphson), Harrell's concordance.

**Validation of the implementation:** the custom Cox fit was checked against `lifelines.CoxPHFitter` on a 3,000-sample synthetic cohort — β agreed to **1.4 × 10⁻⁷** and SE to **5.7 × 10⁻¹¹**; concordance agreed to 3 × 10⁻⁴.

**Calibration slope** is obtained by fitting the outcome on the linear predictor: the coefficient is 1.0 under perfect calibration. Its SE is the quantity targeted.

---

## 3. Assumptions — every one is a choice, not a fact

| Assumption | Value used | Direction of error if wrong |
|---|---|---|
| **Linear-predictor distribution** | Normal(0, σ²) | Riley's framework requires this; skewed LPs generally need **more** events |
| **Model strength** | σ calibrated to **C = 0.60, 0.65, 0.70** | spans what CRC recurrence signatures typically report; **weaker models need far more events** |
| **Baseline hazard** | exponential | Weibull/flexible baselines would change absolute n, not the ordering |
| **Censoring** | administrative only, at the quantile giving the target event fraction | **random/informative censoring reduces effective information → more events needed** |
| **Event fraction** | 0.30 (sensitivity 0.20–0.40, §5.3) | conclusion insensitive (§5.3) |
| **True calibration slope** | **1.0 — perfect calibration** | **OPTIMISTIC.** A genuinely miscalibrated model needs **more** events, not fewer |
| **Follow-up** | complete to the administrative horizon | incomplete follow-up reduces effective events |

**Every assumption here is favourable to feasibility.** The results below are therefore a **best case**.

### 3.1 Assumptions deliberately NOT made

- **No time unit is assumed.** The simulation is on an arbitrary time scale and all results are reported in **events**, which is unit-free. The unresolved `rfs.delay` unit (Phase 1 U1) does not affect any number here.
- **No endpoint harmonisation.** GSE33113's composite endpoint is assessed **separately** and its events are **never added** to GSE14333's.
- **No death-without-recurrence rule assumed.** That rule (Phase 1 U2) affects *which* records are events in the real cohorts; it does not affect the precision arithmetic given an event count. Its unresolved status is carried forward in §6.

---

## 4. Internal vs external validation feasibility — different questions

These are **not** the same assessment and must not be conflated.

### 4.1 Internal model-development and optimism-correction feasibility

**This is a different question from external-validation precision and is assessed separately.**

GSE39582 S3 has **136 events** (130 under the conditional 60-unit endpoint). Internal validation asks whether *development-sample optimism* can be corrected — addressed by **bootstrap optimism correction or repeated cross-validation** over the full sample, not by a held-out calibration-slope CI. The companion development paper is Riley et al., *Stat Med* 2019;38:1276–1296.

**Feasibility depends entirely on model design.** Three designs, assessed separately:

#### (a) Restricted clinical model — **potentially feasible**

Candidate predictors from the Phase 1 audit: age, sex, stage, T, N, location, adjuvant chemotherapy. Roughly **8–10 parameters against 136 events**.

Requirements: penalisation (ridge or elastic net) with the penalty tuned **inside** the resampling; bootstrap optimism correction for C-index, calibration slope and calibration-in-the-large; a formal **development** sample-size calculation (Riley 2019) once candidate predictors are fixed. **Treatment is a covariate, never a causal term** (Phase 1 §4.5).

#### (b) Pre-specified pathway-level model — **potentially feasible, with conditions**

A **small, pre-specified** set of pathway or signature scores (order 10–50 summary features, fixed *before* looking at outcome) keeps the effective dimension comparable to (a).

**The pre-specification must be genuine and documented**: the pathway list, its source, and the scoring rule written down before any outcome is touched. A pathway set chosen after inspecting associations is transcriptome-wide selection wearing a smaller label, and inherits the behaviour in (c).

#### (c) Transcriptome-wide feature selection — **NOT feasible at 136 events**

Quantified by simulation (`scripts/26_phase5_methodological_audit.py`, `reports/phase5/selection_optimism_null.csv`). Design: n = 453, ~136 events, **21,815 predictors, every one pure noise**. The true C-index is **0.500** for any model built from them.

| Genes selected by univariate screening | Apparent C | **Optimism** |
|---:|---:|---:|
| 5 | 0.7011 | **+0.201** |
| 10 | 0.7580 | **+0.258** |
| 25 | 0.8360 | **+0.336** |
| 50 | 0.8862 | **+0.386** |
| **100** | **0.9173** | **+0.417** |

**Selecting 100 pure-noise genes produces an apparent C-index of 0.92 when the truth is 0.50.** Every point above 0.500 is selection-induced optimism.

The mechanism, from the same simulation (one cohort, 21,815 null genes):

| Statistic | Value |
|---|---:|
| median \|z\| | 0.667 |
| 95th percentile | 1.970 |
| genes with \|z\| > 1.96 | **1,121 (5.1%** — exactly the nominal rate) |
| genes with \|z\| > 3.00 | **54** |
| **maximum \|z\|** | **4.267** (nominal two-sided p = 2.0 × 10⁻⁵) |

**With ~21,800 null genes, several always reach \|z\| > 4 by chance.** At 136 events, univariate screening cannot distinguish them from real signal — the noise maximum sits where a genuine effect would.

**Consequences, stated plainly:**

1. An **apparent** C-index from transcriptome-wide selection at this event count is **uninterpretable** without correction.
2. Correction must be **nested** — the selection step repeated inside every bootstrap/CV fold. Selecting once on the full data and cross-validating afterwards does not correct this; it measures nothing.
3. Even nested correction yields a **very imprecise** corrected estimate at 136 events, since the C-index CI width is ~0.095 (§5.4) before accounting for selection variability.
4. **Recommendation: do not pursue transcriptome-wide selection as the primary design.** If pursued at all, it is exploratory hypothesis generation, must use nested resampling, and must report the null-simulation optimism above alongside any result.

**Important caveat for all three designs:** internal validation estimates **reproducibility**, not **transportability**. It says nothing about performance in a new population — which is the project's actual question.

### 4.2 CIT partition — not an external test

Already established (Phase 4 §5): partition and processing batch are confounded (Cramér's V = **0.488**, p = 1.0 × 10⁻¹⁹; two batches discovery-only; 56% of S3 in ≥90%-skewed batches). The CIT split is an internal holdout at most, and with **29 events** it would fail precision regardless.

### 4.3 External validation — not feasible for confirmatory claims

§5 below. **GSE14333 (48 events) and GSE33113 (18 events) are assessed separately and never pooled.**

---

## 5. Results

### 5.1 Calibration-slope precision by event count

Mean SE and 95% CI width of the calibration slope, 400 replications per cell:

| Events | C = 0.60 | C = 0.65 | C = 0.70 |
|---:|---:|---:|---:|
| **18** | 2.601 | **1.791** | 1.346 |
| **48** | 1.558 | **1.071** | 0.795 |
| 100 | — | 0.731 | 0.545 |
| **130** | 0.934 | **0.643** | 0.477 |
| **136** | 0.917 | **0.625** | 0.465 |
| 200 | 0.751 | 0.517 | 0.383 ✅ |
| 365 | 0.557 | 0.383 ✅ | 0.283 ✅ |
| 500 | 0.475 | 0.326 ✅ | 0.241 ✅ |
| 1000 | 0.335 | 0.231 ✅ | 0.171 ✅✅ |

✅ = meets relaxed target (0.4) · ✅✅ = meets primary target (0.2)

### 5.2 Events required to meet each target

Fitted as SE ≈ k/√E across the simulated grid:

| Model strength | SE relationship | Events for **width 0.4** (relaxed) | Events for **width 0.2** (primary) |
|---|---|---:|---:|
| C = 0.60 | SE ≈ 2.732/√E | **~717** | ~2,870 |
| C = 0.65 | SE ≈ 1.875/√E | **~338** | ~1,351 |
| C = 0.70 | SE ≈ 1.393/√E | **~186** | ~746 |

> **No cross-check is claimed.** Revision 1 compared ~338 here with Riley's ~365 in their VTE example. **That comparison is withdrawn** (§0.1): different estimator, different disease, different LP distribution, different censoring, different absolute risk. The proximity is coincidental.
>
> **These numbers are indicative, not Riley-derived.** Read them as order-of-magnitude requirements from a Cox-on-LP approximation (§0.3).

### 5.3 Sensitivity to the event-fraction assumption

| Event fraction | GSE14333 (48 ev) slope SE | GSE39582 S3 (136 ev) slope SE |
|---:|---:|---:|
| 0.20 | 0.2700 (n = 240) | 0.1589 (n = 680) |
| 0.30 | 0.2730 (n = 160) | 0.1605 (n = 453) |
| 0.40 | 0.2759 (n = 120) | 0.1623 (n = 340) |

**The conclusion is insensitive.** Precision is governed by the **event count**, not the sample size — doubling n at a fixed event count barely changes the SE. This is the central reason "GSE14333 has 290 samples" is not a relevant number.

### 5.4 Discrimination precision

95% CI width of Harrell's C (empirical SD across replications), C = 0.65:

| Events | C 95% CI width | Interpretation |
|---:|---:|---|
| 18 | **0.269** | a true C of 0.65 could read 0.52–0.78 |
| 48 | **0.144** | a true C of 0.65 could read 0.58–0.72 |
| 136 | **0.097** | 0.60–0.70 |
| 365 | 0.057 | 0.62–0.68 |

At 48 events, a C-index estimate **cannot distinguish a useful signature (0.70) from a near-useless one (0.58)**.

### 5.5 What a 48-event calibration slope actually admits

At 48 events and C = 0.65, the expected 95% CI is **[0.48, 1.55]**. That single interval simultaneously admits:

| Slope | Meaning |
|---|---|
| ~0.5 | **severe over-fitting** — predictions far too extreme |
| 1.0 | **perfect calibration** |
| ~1.5 | **under-fitting** — predictions too conservative |

**An interval compatible with all three distinguishes nothing.** This is the concrete reason 48 events cannot support a validation claim — not a threshold, but the width of the answer.

---

## 6. Decision table — four categories

**❌ unsupported** · **⚠️ exploratory** (may be done, must be labelled as such) · **🔶 potentially feasible under a restricted design** · **❓ unknown until a model and horizon are specified**

### 6.1 External validation

| Claim | GSE14333 (48 ev) | GSE33113 (18 ev) |
|---|---|---|
| Confirmatory external validation | ❌ **unsupported** | ❌ **unsupported** |
| Calibration slope estimate | ❌ **unsupported** (width ≈ 1.07) | ❌ **unsupported** (≈ 1.79) |
| Calibration curve across risk range | ❌ **unsupported** | ❌ **unsupported** |
| C-index with CI | ⚠️ **exploratory**, wide (±0.07) | ❌ **unsupported** (±0.13) |
| Time-dependent performance | ❌ **unsupported** — splits events further | ❌ **unsupported** |
| Qualitative direction-of-effect check | ⚠️ **exploratory** | ⚠️ **exploratory** |
| Pooling 48 + 18 events | ❌ **prohibited** — different endpoint constructs, time units, stage coverage | ❌ **prohibited** |
| **Faithful Riley sample-size calculation** | ❓ **unknown** — needs horizon, predicted risks, LP shape, F(t), censoring distribution (§0.4) | ❓ **unknown** |

### 6.2 Internal development and optimism correction (GSE39582 S3, 136 events)

| Design / claim | Status |
|---|---|
| **Restricted clinical model** (~8–10 params, penalised, nested tuning) | 🔶 **potentially feasible** |
| **Pre-specified pathway-level model** (~10–50 fixed features) | 🔶 **potentially feasible**, conditional on genuine documented pre-specification |
| **Transcriptome-wide selection** (21,815 genes) | ❌ **not feasible** — null simulation gives apparent C 0.92 at true C 0.50 (§4.1c) |
| Transcriptome-wide as **exploratory hypothesis generation** | ⚠️ **exploratory only**, nested resampling mandatory, null-optimism reported alongside |
| Bootstrap / repeated-CV optimism correction | 🔶 **feasible** for designs (a) and (b) |
| Apparent (uncorrected) performance | ❌ **must never be reported as performance** |
| Internal C-index with CI | ⚠️ **exploratory** (±0.05 before selection variability) |
| Internal calibration slope | ⚠️ **exploratory** (width ≈ 0.63) |
| Development sample-size calculation (Riley 2019) | ❓ **unknown** until candidate predictors are fixed |
| **Transportability claim of any kind** | ❌ **unsupported** — internal validation cannot establish it |

### 6.3 Blocked on unresolved endpoint questions

| Item | Blocked by | Affects |
|---|---|---|
| Any **horizon-specific** analysis (t-year risk, time-dependent AUC/Brier, faithful Riley calculation) | **U1 — time unit unresolved** | a horizon cannot be expressed in `rfs.delay` units |
| The **60-unit conditional endpoint** | **U1** | "60 units" is only "5 years" if the unit is months |
| **Which records count as events** | **U2 — death-without-recurrence rule** | event counts, and therefore every precision figure |
| Event-count arithmetic **in events** | **neither** | §5 is unit-free and unaffected by U1 |

---

## 7. Limitations

1. **The simulation is optimistic by construction.** Perfect calibration, normal LP, administrative censoring only, complete follow-up. Every departure increases the events required.
2. **Model strength is assumed, not known.** No model has been built, so C = 0.60–0.70 is a plausible range from the literature, not a measurement. If the real signature achieves C = 0.60, requirements roughly **double** relative to C = 0.65.
3. **Exponential baseline hazard.** A different baseline would shift absolute numbers; it would not change which cohorts fall short by factors of 7–40×.
4. **Harrell's C under censoring** is known to depend on the censoring distribution; the CI widths in §5.4 are indicative.
5. **This is an approximation of Riley's framework, not an implementation of it** (§0.2). The calibration-slope estimator, the scale, the variance estimator, the censoring mechanism and the absence of a prediction horizon all diverge. Absolute event requirements are **indicative**; the order-of-magnitude shortfall is robust because it follows from √E scaling.
6. **Harrell's C is not censoring-robust.** Uno's IPCW C is not implemented. The C-index CI is a Monte-Carlo sampling SD, not Riley's mean within-dataset width — though it tracks a bootstrap CI to within 0–8% (§0.5).
7. **A faithful Riley calculation is currently impossible** — six required inputs are unavailable until a model and horizon exist (§0.4).
8. **The unresolved endpoint issues are carried, not solved.** The time unit (U1) and death-without-recurrence rule (U2) remain unresolved. They do not affect the precision arithmetic — which is in events — but they do affect *which records count as events*, and therefore remain prerequisites for any analysis.
9. **Shortfall magnitude, stated plainly:** at C = 0.65, GSE14333 has **48 of ~338** required events (14%); GSE33113 has **18 of ~338** (5%).

---

## 8. Recommendations

**8.1 Do not pursue confirmatory external validation with the identified cohorts.** The shortfall is 7–40×, not a marginal miss. No preprocessing choice, cohort combination or analytic refinement closes it.

**8.2 Proceed with a single-cohort, internally-validated study.** GSE39582 S3 (453 patients, 136 events) supports development with bootstrap optimism correction, mandatory shrinkage, and feature selection nested inside the resampling. Report it as **internally validated**, with transportability explicitly unestablished.

**8.3 If external data is used at all, label it exploratory.** GSE14333 and GSE33113 may support a **qualitative direction-of-effect check** — reported per cohort, with each cohort's endpoint construct, time unit and stage coverage attached, and never described as validation. **Their events must not be combined.**

**8.4 Settle the endpoint prerequisites first** (U1 time unit, U2 death handling). They are cheap relative to everything else and gate any analysis.

**8.5 Revisit only if a materially larger cohort appears.** A GPL570 CRC cohort with **~340+ stage II–III recurrence events** and a compatible endpoint would change this conclusion. GSE92921, GSE31595 and GSE29621 remain unassessed but are 37–65 samples — they cannot supply it.

**8.6 The fRMA pilot is now lower priority.** Uniform preprocessing was needed to make external validation defensible. With external validation not feasible, it is only required if the exploratory assessment in 8.3 is pursued. **The hard gate (`PHASE4_FRMA_PILOT.md` §3) still costs nothing and can be run whenever convenient.**

---

## 9. Integrity

| Check | Result |
|---|---|
| Phase 0 + Phase 3 test suite | **57/57 passed** |
| Raw data SHA-256 | both **OK** |
| Phase 3 bundle hash | `a1ca41807c534584` — **unchanged** |
| Expression data read | **none** |
| Model fitted to real data | **none** |
| Packages installed | **none** |

**Files created:** `src/precision_sim.py` · `scripts/25_validation_precision_simulation.py` · `reports/phase5/` (4 files) · this document.
**Files modified:** none.
