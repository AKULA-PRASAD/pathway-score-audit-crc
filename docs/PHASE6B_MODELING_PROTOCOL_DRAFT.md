# Phase 6B — Modeling Protocol (DRAFT, rev. 2)

**Date:** 2026-11-02 · **Status:** **DRAFT STRUCTURE — approved as structure, NOT execution-ready.** Not executed. No model trained, no expression normalized, no CEL files, no packages installed, Phase 3 bundle untouched.

---

## 1. Scope and estimand

**Question:** do transcriptomic features add prognostic information beyond routinely available clinical variables for recurrence in stage II–III colon cancer?

**Estimand (primary):** the **cause-specific hazard ratio** for recurrence, and **discrimination** (C-index), in the S3 cohort.

**Explicitly NOT the estimand:** absolute t-year recurrence risk. Blocked by U1 (unit) and complicated by U2 (competing risks). See `PHASE6A_ENDPOINT_RESOLUTION.md` §3.

**Validation status:** internal only. Phase 5 established that confirmatory external validation is not feasible with available cohorts. **No transportability claim will be made.**

---

## 2. Fixed inputs — all frozen before any modelling

| Component | Specification | Frozen in |
|---|---|---|
| **Cohort (primary)** | **S3** — 453 samples, 136 events, 354 discovery / 99 validation | Phase 1 |
| Cohort sensitivity | S0 (461/139), S1 (458/136), S2 (456/139) | Phase 1 |
| **Representation (primary)** | **median aggregation**, 21,815 genes × 43,163 probes | Phase 2B |
| **Representation (sensitivity)** | **representative-IQR selected on discovery_443**, 21,815 genes | Phase 2B |
| Endpoint (primary) | `rfs_event_full` / `rfs_time_full`, cause-specific | Phase 6A |
| Endpoint sensitivity | S-1 competing risks · S-2 composite DFS · S-3 60-unit (conditional) | Phase 6A |
| Analysis bundle | `data/processed/analysis/`, hashes recorded | Phase 3 |

**None of these may be revisited after seeing model results.**

---

## 3. Model families

### 3.1 M1 — Parsimonious clinical baseline *(the benchmark)*

| | |
|---|---|
| **Predictors** | age, sex, stage (II/III), tumour location, adjuvant chemotherapy. **Pre-specified; no data-driven selection.** |
| Optional (pre-specified, only if missingness permits) | T category, N category, MMR status |
| Parameters | ~6–9 against 136 events |
| Model | Cox proportional hazards, cause-specific |
| Penalisation | **AMENDED by P6B-A1 (2026-11-03)** — the clinical block (the four binary covariates and the age spline basis) is **force-included and UNPENALISED** in M1, M2 and M3. The ridge penalty applies **only** to the M2 pathway scores; the elastic-net penalty **only** to the M3 candidate genes. λ (and α for M3) are tuned by inner cross-validation **inside** each training fold. **M1 therefore has no penalised term and no λ to tune.** *Superseded text:* `ridge, λ tuned inside resampling`. See §3.1.1 |
| **Chemotherapy** | **prognostic covariate only. Its coefficient is never interpreted as a treatment effect** (Phase 1 §4.5 — confounding by indication, with the driving confounders absent from the metadata) |
| Stage | included; stratification by stage assessed as a pre-specified alternative |
| Missing data | pre-specified: multiple imputation **inside** resampling, or complete-case with the drop reported. **Decided before execution** |

#### 3.1.1 Amendment P6B-A1 — penalty scope (2026-11-03)

> **The original wording is preserved verbatim here and is not deleted from the record.**
>
> **Superseded text (original Phase 6B draft, M1 table, *Penalisation* row):**
> `ridge, λ tuned **inside** resampling`
>
> **Replacement:** the clinical block (four binary covariates + the age spline basis) is
> **force-included and unpenalised** in M1, M2 and M3. The ridge penalty applies only to the M2
> pathway scores and the elastic net only to the M3 candidate genes; λ (and α for M3) are tuned by
> inner cross-validation inside each training fold. **M1 therefore has no penalised term and no λ
> to tune.**

| | |
|---|---|
| **Amendment** | **P6B-A1** · approved by the project owner, 2026-11-03 |
| **Basis** | Phase 8D decision D9-1; Phase 9 protocol v2.1 §3.5; Phase 9 Final Audit Finding 3 |
| **Register A item?** | **No.** A-11 fixes the predictor *list*, A-15 the model *roles*, A-17 the resampling. Penalisation **scope** was protocol text only |
| **Rationale** | Ridge is not scale-invariant. The clinical block is four binaries plus a **year-scale** age spline basis; penalising them makes shrinkage depend on *units*. An **unpenalised** coefficient is **scale-equivariant**, so the question disappears rather than being answered by convention |
| **Outcome data used** | **No** |

**Model definitions after P6B-A1 — the clinical block is identical in all three:**

| Model | Clinical block | Molecular block |
|---|---|---|
| **M1** | the five pre-specified clinical predictors, **unpenalised** | *none* |
| **M2** | **the same** clinical predictors, **unpenalised** | four (primary) / five (sensitivity) pathway scores, **ridge-penalised**, λ tuned by inner CV inside the training fold |
| **M3** | **the same** clinical predictors, **unpenalised** | selected transcriptomic genes, **elastic-net-penalised**, α and λ tuned by inner CV inside the training fold |

**Nesting.** M1 is exactly M2 with all molecular coefficients set to zero, on the same data and the
same risk sets. λ → ∞ recovers M1. The clinical block is **never separately penalised or retuned**
between models, which is what makes the nesting hold.

> ### ⚠ Required interpretation of ΔC — binding
>
> **ΔC = C(M2) − C(M1)** compares an **unpenalised clinical baseline** against a model with an
> **unpenalised clinical block plus a shrunken molecular block**. Molecular shrinkage can make the
> comparison **conservative** with respect to detecting an incremental molecular contribution.
>
> **ΔC must not be described as an unbiased estimate of the molecular block's contribution.**

---

**M1 is built first and never re-tuned afterwards.** Every genomic model is measured against this fixed benchmark.

### 3.2 M2 — Pre-specified pathway-level model

**This section is a SPECIFICATION OF REQUIREMENTS, not a completed specification.** M2 cannot be executed until every field below is filled in and frozen to a file. The fields are deliberately left blank: filling them is a separate, approvable step.

#### 3.2.1 ✅ SPECIFICATION COMPLETE AND FROZEN — 2026-11-02T21:13:12.824057+00:00

All seven required fields are filled. **M2 is frozen** (Phase 6C rev. 3).

| # | Field | Value |
|---|---|---|
| **1** | Database and version | **MSigDB Human Hallmark (H), v2026.2.Hs**, released **2026-10-15**. GMT SHA-256 `56bba1300555fb6302844fc8…` (symbols), `a26578a72aed3171fffce4d2…` (entrez) |
| **2** | Pathway IDs and gene membership | **PRIMARY (4-set):** EMT 198 · E2F_TARGETS 195 · WNT_BETA_CATENIN 41 · P53_PATHWAY 194 — union **607**. **SENSITIVITY (5-set):** + PI3K_AKT_MTOR 104 — union **702**. Full membership and true-unmatched lists in the spec files |
| **3** | Biological rationale | Phase 6C §4 — four distinct axes: stromal/mesenchymal (CMS4 representative), proliferation, canonical Wnt driver, DNA-damage/checkpoint |
| **4** | Overlap handling | Retained as-is, overlap quantified. Gene-level max Jaccard 0.054; **score-level** max \|r\| **0.318** (4-set), condition number **2.95**, 0 pairs > 0.5 |
| **5** | Scoring method | **mean of training-fold-standardised member genes**; numpy 2.5.4 / pandas 2.3.3; **no new packages** |
| **6** | Standardization | per-gene μ̂/σ̂ from **training fold only**; σ̂ is the **sample** SD, `ddof=1` (**pinned by amendment M2-A1**, 2026-11-02 — read from `scoring.sigma_ddof`, never from a library default); no score-level standardisation before the Cox model |
| **7** | Fold-specific calculation | **confirmed** — membership is frozen metadata (outside the loop); every estimated constant is re-fitted per training fold |

**Spec files (revision 2, after amendment M2-A1):** `config/pathways/m2_panel_4set.json` (`cb4de071e2f64a6b…`, superseding `bab99da9f139db34…`) · `m2_panel_5set.json` (`10dae4874a918160…`, superseding `6058f733b2acda1c…`) · `M2_FREEZE_RECORD.json`

> **PI3K/AKT/mTOR is sensitivity-only and must NOT be promoted to primary or removed from sensitivity on the basis of observed model performance.** Either change reclassifies M2 as M3-exploratory (§3.2.4).

#### 3.2.2 Fold-specific computation — the binding rule

> **Pathway scores are NOT computed once on the full cohort.** Any scoring step that involves an estimated quantity — a mean, an SD, a PCA rotation, an ssGSEA normalisation constant — is a model-fitting step and **must be re-estimated inside every bootstrap replicate / CV fold, from the training portion only**, then applied unchanged to the held-out portion.

| Scoring component | Where computed |
|---|---|
| **pathway membership** (which genes belong) | **outside** — it is frozen metadata, not an estimate |
| gene-level centring/scaling constants | **inside**, from training fold |
| PCA rotation / loadings (if used) | **inside**, from training fold |
| ssGSEA normalisation constants (if used) | **inside**, from training fold |
| score-level standardisation | **inside**, from training fold |
| penalty tuning over the scores | **inside**, nested within training fold |

A score computed once on all 453 samples and then cross-validated **leaks held-out information into the feature definition** and makes the resulting performance uninterpretable.

#### 3.2.3 Model

M1 clinical terms **+** pathway scores, penalised Cox (cause-specific). Incremental value over M1 reported as ΔC-index with CI, likelihood-ratio test, and calibration — all optimism-corrected.

#### 3.2.4 Pre-specification is binding — and cannot be retrofitted

> **Pathways must NOT be selected using outcome associations and subsequently labelled "pre-specified."** This is the single most important constraint on M2.

Operationally:

1. The seven fields above are written to a **version-controlled file with a timestamp**, before any outcome column is joined to any expression data.
2. **Any** of the following reclassifies M2 to **M3-exploratory**, permanently and in the written record:
   - adding, removing or substituting a pathway after any outcome contact;
   - filtering the list by association, p-value, variance-with-outcome, or model fit;
   - choosing among candidate lists by which performs better;
   - re-running with a different scoring method and reporting the better result.
3. Reclassification is **not reversible** by reverting the change — once outcome information has informed the list, the list carries that information.
4. If the pathway list is revised for a **non-outcome** reason (e.g. a database version correction, or genes absent from our index), the revision, its reason and its timestamp are recorded, and the prior version is retained.

### 3.3 M3 — Transcriptome-wide *(EXPLORATORY ONLY)*

| | |
|---|---|
| **Status** | **exploratory hypothesis generation. Not a candidate for the primary result.** |
| Features | all 21,815 genes |
| Method | penalised Cox (elastic net), **entire pipeline nested inside resampling** |
| **Mandatory reporting** | the Phase 5 null simulation **alongside every result**: at 136 events with 21,815 null predictors, selecting the top 100 by univariate screening gives apparent C = **0.917** against a true C of **0.500** |

**Why it cannot be primary** (Phase 5 §4.1c): with ~21,800 null genes the maximum \|z\| reaches **4.267** by chance alone, with 54 genes above \|z\| = 3. At 136 events, univariate screening cannot distinguish noise maxima from real signal.

**Any apparent M3 performance that is not optimism-corrected is uninterpretable and must never be reported as performance.**

---

## 4. Leakage prevention — what happens inside the resampling loop

> **Rule: anything estimated from data must be re-estimated inside every fold/bootstrap replicate, using only that replicate's training portion.**

### 4.1 Inside the loop (re-estimated every replicate)

| Step | Why it must be inside |
|---|---|
| **Gene filtering** (expression/variance thresholds) | thresholds computed on all data see the held-out samples |
| **Standardisation / centring / scaling** | means and SDs are estimates |
| **Pathway score computation** | if scoring standardises member genes, that is an estimate |
| **Univariate screening** (M3) | **the dominant leakage source** — this is what produces C = 0.92 from noise |
| **Penalty tuning** (λ, α) | inner cross-validation within the training portion only |
| **Imputation models** | fitted on training only, applied to held-out |
| **Coefficient estimation** | obviously |
| **Any threshold or cut-point** | including risk-group cut-offs |

### 4.2 Outside the loop (fixed inputs, never re-estimated)

| Item | Why it is legitimately outside |
|---|---|
| Cohort membership (S3) | frozen in Phase 1 from eligibility criteria, **outcome-blind by construction** |
| Gene→probe representation | frozen in Phase 2B; selection used **discovery_443 only** and is verified validation-free (test T-29) |
| The pathway **list** (M2) | pre-specified before any outcome contact |
| The clinical predictor list (M1) | pre-specified |
| Endpoint definition | frozen in Phase 6A |

### 4.3 The specific trap to avoid

**Selecting features once on the full dataset and then cross-validating the fitted model does not correct selection optimism — it measures nothing.** The held-out samples already influenced which features exist. Phase 5's null simulation quantifies exactly this failure.

### 4.4 A note on the frozen representation

The representative-IQR map was selected on `discovery_443` and is therefore **not** outcome-blind with respect to *dataset membership* — but it is outcome-blind with respect to **outcome values**, which is what matters for leakage. The median representation involves no selection at all. Verified by tests T-18/T-19/T-20/T-29.

---

## 5. Validation and optimism correction

| | |
|---|---|
| **Method** | bootstrap optimism correction (Harrell's approach), with repeated 10-fold CV as a cross-check |
| **B** | **B = 500 is an INITIAL specification, not a final one.** It is a starting value, subject to the Monte-Carlo stability checks in §5.1. The final B is whatever those checks require |
| **Scope** | the **entire** pipeline of §4.1 is repeated inside every replicate |
| **Metrics** | C-index (optimism-corrected) · calibration slope (log-hazard scale — unit-free, §6A) · calibration-in-the-large · ΔC-index vs M1 with CI |
| **Shrinkage** | uniform shrinkage factor estimated from the bootstrap and **reported**, whether or not applied; penalised models report effective degrees of freedom |
| **Reported together** | apparent **and** corrected performance, with optimism shown explicitly. Apparent performance is **never** reported alone |
| **Precision caveat** | at 136 events the corrected C-index CI is ~±0.05 **before** selection variability (Phase 5 §5.4). Every estimate carries its CI |

### 5.1 Monte-Carlo stability of the bootstrap — required, not optional

B = 500 is a **starting point**. The resampling itself introduces Monte-Carlo error, which is separate from — and additional to — the sampling uncertainty already quantified in Phase 5. It must be measured, not assumed negligible.

**Required checks, run before any result is reported:**

| # | Check | Criterion |
|---|---|---|
| **MC-1** | **Running-estimate trace.** Plot the optimism-corrected C-index and calibration slope as a function of B (B = 50, 100, 200, 300, 400, 500, …) | the estimate must be visibly stable over the final ~40% of replicates |
| **MC-2** | **Repeat-with-different-seeds.** Run the whole B-replicate procedure ≥ 5 times with different seeds | the **SD of the corrected estimate across seeds** is the Monte-Carlo SE |
| **MC-3** | **Relative-magnitude rule** | Monte-Carlo SE must be **≤ 10% of the sampling SE** of the same quantity. If not, **increase B and repeat** |
| **MC-4** | **Report it** | the Monte-Carlo SE is reported alongside every corrected estimate, explicitly distinguished from sampling uncertainty |

**If MC-3 fails at B = 500**, B is increased (1000, 2000, …) until it passes. **B is therefore an output of the stability analysis, not an input fixed by this protocol.** The *rule* (MC-1…MC-4 and the 10% criterion) is frozen now; the *value of B* is not.

M3's nested selection is the most computationally expensive and most Monte-Carlo-unstable case, and will likely require the largest B.

---

**Not performed:** time-dependent AUC/Brier (blocked by U1) · external validation (not feasible) · a faithful Riley cloglog calibration (blocked by U1).

---

## 6. Two explicit decision registers

### 6.1 REGISTER A — frozen before any outcome-driven analysis

> **These must be fixed before any outcome column is joined to any expression data, and must not change afterwards.** A change to any of them after outcome contact invalidates the pre-specified status of the affected analysis.

| # | Decision | Value | Frozen in |
|---|---|---|---|
| A-01 | Primary cohort | **S3** (453 samples, 136 events) | Phase 1 |
| A-02 | Sensitivity cohorts | S0 (461/139), S1 (458/136), S2 (456/139) | Phase 1 |
| A-03 | Cohort eligibility criteria | Criterion 1 + Criterion 2, kept separate | Phase 1 |
| A-04 | Primary representation | **median aggregation**, 21,815 genes | Phase 2B |
| A-05 | Sensitivity representation | **representative-IQR @ discovery_443** | Phase 2B |
| A-06 | Gene→probe mapping | frozen specs, hash-recorded | Phase 2B |
| A-07 | Primary endpoint | cause-specific recurrence, full follow-up | Phase 6A |
| A-08 | Endpoint sensitivity set | S-1 CIF/Fine–Gray · S-2 composite DFS (182) · S-3 60-unit (conditional) · S-4…S-7 | Phase 6A |
| A-09 | Estimand of the primary analysis | cause-specific hazard ratio + discrimination | Phase 6A |
| A-10 | Prohibition on absolute-risk / t-year claims | binding while U1 unresolved | Phase 6A |
| A-11 | M1 clinical predictor list | age, sex, stage, location, adjuvant chemo (+ pre-specified optionals) | §3.1 |
| A-12 | Chemotherapy is a covariate, never a causal term | binding | Phase 1 §4.5 |
| A-13 | **M2 pathway list — all 7 fields of §3.2.1** | ✅ **FROZEN** — 4-set primary / 5-set sensitivity | Phase 6C, 2026-11-02 |
| A-14 | M2 scoring method and parameters | ✅ **FROZEN** — training-fold z-mean | Phase 6C, 2026-11-02 |
| A-15 | Model families and their status (M1 benchmark · M2 candidate · M3 exploratory) | §3 | §3 |
| A-16 | Leakage rules — what is inside vs outside resampling | §4 | §4 |
| A-17 | Resampling scheme + Monte-Carlo stability rule (MC-1…MC-4, 10% criterion) | §5, §5.1 | §5 |
| A-18 | Primary metric | optimism-corrected C-index; ΔC vs M1 | §5 |
| A-19 | **Missing-data strategy** | ✅ **FROZEN** — complete-case for M1/M2/M2-sens/M3; MI not required | Phase 6D, 2026-11-02 |
| A-20 | Three-framing comparison rules | Phase 6A §4.4 | Phase 6A |

**A-13, A-14 (Phase 6C) and A-19 (Phase 6D) are all FROZEN. Every Register A item is now settled.**

### 6.2 REGISTER B — may legitimately depend on QC findings

> **These may be informed by data, but only by data that is blind to outcome.** Each carries a constraint that prevents it becoming a back-door for outcome-driven choice.

| # | Decision | Permitted basis | Binding constraint |
|---|---|---|---|
| B-01 | **Array exclusion** | QC metrics **on the deposited matrix** — **NOT GNUSE**, which requires fRMA and is inapplicable under Rank 1 | **Specified in Phase 7B** (`PHASE7B_QC_PROTOCOL.md`). Thresholds pre-specified before any QC distribution was computed. Outcome must never inform exclusion. Excluded arrays enumerated |
| B-02 | Expression/variance filter **cut-point** | distribution of probe intensities **in the training fold** | the *rule* is frozen now; only the *value* is data-derived, and it is recomputed **inside** each fold |
| B-03 | ~~fRMA vs conventional RMA~~ → **RESOLVED** | ✅ **Rank 1 APPROVED (2026-11-02): deposited matrix as-is.** No CEL acquisition, no fRMA, no re-normalisation, no re-ComBat (Phase 7A §0) | frozen artifacts unchanged; cross-cohort claims foreclosed |
| B-04 | Batch-structure handling | QC of the 9 recovered labels | **descriptive only.** Any correction is a separate decision requiring approval |
| B-05 | Proportional-hazards remedy | Schoenfeld residual diagnostics | the *diagnostic* is pre-specified; any *remedy* is reported as a protocol deviation |
| B-06 | Final **B** for the bootstrap | Monte-Carlo stability checks (§5.1) | B increases until MC-3 passes; it never decreases to suit a result |
| B-07 | Penalty values (λ, α) | inner cross-validation | **inside** the training fold only |
| B-08 | Imputation model specification | observed missingness pattern | fitted **inside** folds; the *choice* of strategy is A-19 and frozen first |

### 6.3 NEVER PERMITTED

| |
|---|
| Changing any Register A item after seeing outcome-related results |
| **Selecting pathways by outcome association and labelling them pre-specified** (§3.2.4) |
| Reporting apparent (uncorrected) performance as performance |
| Any absolute-risk or t-year claim while U1 is unresolved |
| Pooling external cohort events |
| Describing any analysis here as external validation |
| Claiming that agreement among the three endpoint framings proves death handling is immaterial (Phase 6A §4.4) |
| Reducing B to make an estimate look more stable |

## 6.4 A-19 — Missing-data strategy ✅ **FROZEN**

> **FROZEN 2026-11-02T21:47:11.883173+00:00** · record: `config/analysis/A19_FREEZE_RECORD.json`
>
> ## **COMPLETE-CASE ANALYSIS** is the primary missing-data strategy for **M1, M2 primary (4-set), M2 sensitivity (5-set) and M3.**
> **Multiple imputation is NOT required for the current model specifications.**

### 6.4.1 The frozen decision and its evidence

| Model | complete cases | % | excluded |
|---|---:|---:|---:|
| **M1 primary (clinical only)** | **452 / 453** | **99.8** | **1** |
| **M2 primary (frozen 4-set panel)** | **452 / 453** | **99.8** | **1** |
| **M2 sensitivity (frozen 5-set panel)** | **452 / 453** | **99.8** | **1** |
| **M3 exploratory (transcriptome-wide)** | **452 / 453** | **99.8** | **1** |

Basis: the **Phase 6D definitive missingness audit** (`docs/PHASE6D_MISSINGNESS_AUDIT.md`). Expression inputs carry **0 NaN and 0 zero-variance genes**, so pathway and transcriptome predictors are complete.

#### The single excluded observation — identified, never silently discarded

| | |
|---|---|
| **sample_id** | **GSM972262** |
| title | CIT306 |
| dataset | discovery |
| combat_batch | batch05 |
| **missing variable** | **`chemotherapy.adjuvant`** (raw value `N/A`) |
| observed values | age 88, Female, stage 2, distal |

Recorded in `reports/phase6d/complete_case_exclusions.{json,csv}` and in the freeze record. **Any analysis using the complete-case cohort must report this exclusion.**

#### Prohibited imputation — unchanged and binding

**No imputation of:** expression values · pathway members · endpoint fields · sample identifiers · cohort membership · provenance fields.

#### Frozen M2 coverage rules — unchanged

The M2 specification is untouched: missing gene members are **excluded from the score and enumerated, never substituted**, and the frozen rule (**≥80% of resolved members AND ≥20 genes** surviving per-fold zero-variance exclusion) applies as frozen. Verified byte-identical: `m2_panel_4set.json` `bab99da9…`, `m2_panel_5set.json` `6058f733…`.

### 6.4.1b Contingency — optional covariates remain excluded from M1

`tnm.t`, `tnm.n`, `tnm.m`, `mmr.status`, `cimp.status`, `cin.status`, `kras.mutation`, `braf.mutation`, `tp53.mutation` are **excluded from M1** and remain so.

**If a future model proposes adding any of them, all five are required before the model specification changes:**

1. reassess the scientific rationale **and** the missingness;
2. recalculate complete-case retention (adding TNM/MMR alone would drop it to **390/453, 86.1%**);
3. reconsider whether multiple imputation is justified;
4. **account for assay/batch structure** — molecular missingness is near-deterministic in processing batch;
5. **obtain explicit approval.**

> **Multiple imputation must NOT be activated automatically.** Nothing in this protocol causes MI to begin without a new approval.

### 6.4.1c Deferred — pooling across imputations

> **Status: DEFERRED — NOT APPLICABLE.**
>
> Pooling across imputations (Rubin's rules on coefficients vs averaging individual predictions) has **no role under a complete-case strategy**. **No pooling strategy is specified, and none is invented.** The question reopens only if multiple imputation is proposed **and** approved in a future analysis.

### 6.4.1d Preserved scientific distinctions (Phase 6D)

| Finding | Status |
|---|---|
| The deposit documents `ND: not determined` and `NA: Not applicable`, **but `ND` is used 0 times** — both meanings collapse to `N/A` | preserved |
| **Structural `N/A` values must NOT be treated as ordinary missing observations** (e.g. chemotherapy type when no chemotherapy was given; mutation detail columns when WT) | preserved |
| Molecular assay availability is **strongly associated with batch** — `tp53.mutation` missingness is 0% in four batches and 100% in three | preserved |
| **MCAR is inconsistent** with the evidence; **MAR (conditional on batch) and MNAR/structural remain QUALIFIED ASSESSMENTS, not proven** | preserved |
| The one missing chemotherapy observation (GSM972262) is **explicitly documented** | preserved |

### 6.4.2 Definitive audit specification — COMPLETED (Phase 6D)

Completed 2026-11-02; artifacts in `reports/phase6d/`. All six required outputs delivered:

| # | Required output |
|---|---|
| 1 | Variable-level count and percentage missing |
| 2 | **Missingness patterns** (co-occurrence across variables) |
| 3 | Missingness **by discovery/validation subset**, with the differential |
| 4 | Three categories reported **separately**: missing **clinical covariates** · missing **pathway-score inputs** (gene-level) · missing **endpoint** data |
| 5 | **Complete-case variables listed separately** from those with substantial or structurally missing data |
| 6 | Explicit statement of which variables, if any, actually require imputation |

**No sample or variable is excluded on missingness grounds without a documented, pre-specified rule.** Under the recommendation no such rule is needed: the planned models exclude the high-missingness variables **by model design**, not by a missingness threshold.

### 6.4.3 What may and may not be imputed

| May be imputed | **MUST NOT be imputed** |
|---|---|
| Candidate **clinical predictors**, where the audit shows it appropriate | **Recurrence event status** (`rfs_event_source`) |
| | **Survival/follow-up time** (`rfs_time_source`, `os_*`) |
| | **Sample identifiers** |
| | **Cohort membership** (`in_S0/S1/S2/S3`) |
| | **Provenance fields** (batch, dataset, exclusion criterion, flags) |
| | **Expression values** — see below |

> **Imputation must NOT be used to resolve the documented endpoint ambiguities.** U1 (time unit) and U2 (death-without-recurrence) are **definitional** questions. Imputing around them would manufacture an answer to a question the data cannot answer, and is prohibited.

**M2 pathway scores retain the frozen scoring specification exactly:**

- missing gene members are **excluded from the score**, not replaced;
- **coverage is recorded per fold**;
- the frozen **minimum-coverage rule applies**: ≥80% of resolved members **and** ≥20 genes surviving per-fold zero-variance exclusion;
- **no silent substitution** of missing expression values with zero, mean, or any other arbitrary value.

**Any expression-level imputation requires separate scientific justification and approval.** It is not authorised by this proposal.

### 6.4.4 Fold isolation — binding

| Rule |
|---|
| For **every** training/validation split, the imputation procedure is **fitted on training data only** |
| The **training-derived** procedure is applied to the held-out fold |
| **Never** impute the complete dataset before cross-validation or bootstrap |
| Held-out samples must never influence: imputation parameters · feature filtering · standardisation · pathway-scoring parameters · feature selection · penalty tuning |

This is the same rule as §4.1, extended to imputation. The leakage test of §7.3 applies: *if one more held-out sample were added, would any fitted quantity change?* It must not.

### 6.4.5 Variable-appropriate imputation models

| Requirement |
|---|
| Method chosen **by variable type and distribution** — e.g. predictive mean matching for continuous, logistic/multinomial for categorical |
| **Documented:** predictor set per imputation model · transformations · number of imputations *m* · iterations · convergence diagnostics · random seeds |
| **Deterministic mean/mode substitution is NOT the primary strategy** |
| **No method is specified for a variable whose missingness and structure have not been inspected.** Concretely: the per-variable models are left **unspecified** until §6.4.2 completes |

### 6.4.6 Survival-aware specification

| Rule |
|---|
| For **training-fold** covariate imputation, a survival-aware specification is assessed. If event status and a cumulative baseline-hazard summary (e.g. a Nelson–Aalen estimate) are used as auxiliary variables, **they must be derived exclusively from the training fold** |
| **Held-out predictor imputation must NOT condition on held-out event status or survival time** — these are unavailable at prediction time, and using them would make the held-out performance unachievable in practice |
| The held-out imputation therefore uses **only** held-out predictors plus **training-derived** parameters |

### 6.4.7 Pooling across imputations

> **Pooling coefficients and averaging individual survival predictions are NOT interchangeable**, and the protocol must pre-specify which is used.

| Approach | What it does | When it differs |
|---|---|---|
| **(a) Rubin's rules on coefficients**, then predict once | pools β and its variance across *m* imputations, then forms one linear predictor | |
| **(b) Predict in each imputed dataset**, then average predictions | averages S(t) or the risk estimate across *m* | **These diverge** because S(t) = exp(−H₀(t)·e^LP) is a **non-linear** function of β. Averaging on the coefficient scale ≠ averaging on the probability scale |

**To be decided and recorded at freeze time:** which is primary, and how the Monte-Carlo and between-imputation components of uncertainty are reported. **Not decided here.**

### 6.4.8 Sensitivity analysis

| Comparison | Reported |
|---|---|
| Primary imputation strategy **vs complete-case** | where sample size and missingness permit — for the M1 list this costs only 1 sample (§6.4.1), so it is clearly feasible |
| Changes in **sample size** | yes |
| Changes in **coefficient estimates** | yes |
| Changes in **uncertainty** (CI width) | yes |
| Changes in **discrimination and calibration** | where estimable |

> **Differences are interpreted as sensitivity to missing-data handling — never as proof that one method is correct.** Neither complete-case nor MI is designated "true"; both rest on untestable assumptions about the missingness mechanism, and §6.4.1 gives positive evidence that the mechanism is **not** MCAR for several variables.

### 6.4.9 Status and unresolved decisions

| Item | Status |
|---|---|
| Strategy | ✅ **FROZEN: complete-case**. MI not required |
| Definitive missingness audit | ✅ **DONE** — `PHASE6D_MISSINGNESS_AUDIT.md` |
| Per-variable imputation models | **MOOT** — no imputation performed |
| Number of imputations *m*, iterations, seeds | **MOOT** — no imputation performed |
| Pooling approach (a) vs (b) | **DEFERRED — NOT APPLICABLE** (§6.4.1c) |
| Exclusion rule for high-missingness variables | **NOT NEEDED** — excluded by model design |
| Whether MI is needed **at all** for the planned models | ✅ **RESOLVED: NO** — 452/453 complete for all four |
| **A-19 frozen?** | ✅ **YES — FROZEN 2026-11-02** (complete-case) |

### 6.4.10 Tests

Implemented now (`tests/test_phase6d_a19.py`): policy assertions — A-19 is recorded as PROPOSED not FROZEN; the never-impute list contains the endpoint, identifier, membership and provenance fields; the M2 scoring rule is unchanged; expression imputation is not authorised.

**Not implementable yet:** behavioural fold-isolation and outcome-leakage tests require the imputation code to exist. They are specified here and will be added with the implementation:

| Planned test | Assertion |
|---|---|
| T-A19-1 | imputation parameters fitted on a training fold are **bit-identical** when held-out rows are deleted from the input |
| T-A19-2 | held-out imputation produces identical values whether or not held-out outcome columns are present |
| T-A19-3 | no endpoint, identifier, membership or provenance column is ever written by the imputation step |
| T-A19-4 | pathway scores use exclusion + coverage recording, never substitution |

---

## 7. Reporting commitments

1. **M1 before M2 before M3**, with M3 labelled exploratory throughout.
2. **Apparent and corrected performance side by side**, with optimism and shrinkage shown.
3. **All sensitivity analyses reported**, including those that disagree with the primary.
4. **U1 and U2 stated in the limitations**, with their consequences (§6A §3).
5. **No transportability claim.** Internal validation establishes reproducibility, not transportability.
6. **The Phase 5 null simulation reported alongside any M3 result.**
7. **Protocol deviations recorded explicitly**, with the reason and the stage at which they occurred.

---

## 8. Prerequisites before execution

| # | Prerequisite | Status |
|---|---|---|
| 1 | This protocol approved | **pending** |
| 2 | Pathway list + scoring rule written to a file and frozen | ✅ **DONE** — 2026-11-02 |
| 3 | fRMA identifier hard gate | **MOOT under Rank 1** — fRMA not used. Residual loaded-object check deferred with Rank 2 |
| 4 | Preprocessing strategy decided | ✅ **DONE** — Rank 1, deposited matrix as-is (2026-11-02) |
| 5 | QC thresholds pre-specified with sources | ⚠️ **PROPOSED** — `PHASE7B_QC_PROTOCOL.md`; awaiting approval |
| 6 | Missing-data strategy decided | ✅ **FROZEN** — complete-case (2026-11-02) |
| 7 | U1 correspondence attempted | **not done** — user action |

**Nothing in §3–§5 may begin until 1–6 are complete.**
