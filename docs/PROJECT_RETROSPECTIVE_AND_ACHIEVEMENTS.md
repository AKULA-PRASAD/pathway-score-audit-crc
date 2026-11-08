# Project retrospective: original intent, achievements, findings and remaining gaps

**Date:** 8 November 2026. **Type:** retrospective, documentation only. It uses existing records only: no model was fitted, no analysis was run, and no frozen artifact, configuration, figure or existing report was changed.

**Central question.** Did the project achieve what it set out to do? What did it set out to discover, what did it actually find, what did it build, and what is the work worth scientifically and technically?

**Short answer.**
- The project began as a **cross-cohort transportability** study.
- That aim was formally abandoned at Phases 4–5, when no external cohort with enough events could be found.
- It became an **internally validated, pre-specified test** of whether four Hallmark pathway scores add recurrence discrimination beyond clinical covariates in GSE39582.
- That narrower question was answered. The answer is null: **no demonstrated improvement**, with an interval that cannot exclude a small loss or a modest gain.
- The strongest achievement is the rigour of the answer, not a positive result.

**Conventions used throughout.**
- Section and table references point to files in this repository.
- **Observation** = a frozen numerical result. **Hypothesis** = a literature-based or speculative interpretation. **Unsupported** = a claim the records do not support.
- "Five runs" always means five Monte Carlo runs of one analysis that **differ only in random seed**, never five cohorts or experiments.
- "Independent audit" always means a **separate re-implementation within this project** that recomputed reported quantities from recorded checkpoint values. It did **not** refit models, and it is **not** external replication.

---

## Contents
1. Original intentions, reconstructed
2. Objective-by-objective achievement assessment
3. What we actually discovered
4. Biological interpretation
5. Unexpected findings and lessons
6. What we built
7. What we did not achieve
8. Overall assessment
9. Portfolio translation
10. Evidence sources and unrecoverable information

---

## 1. Original intentions, reconstructed

### 1.1 Earliest records

| Date | Record | What it shows |
|---|---|---|
| 2 Nov 2026 | Commit `6fda652`, "Baseline: existing work prior to Phase 0" (local `main` only; the published branch's root commit `1cd2bba` contains the same baseline scripts) | 19 files: GEO acquisition (scripts 01–12), clinical-metadata and RFS audits, expression completeness, GPL570 annotation and probe–gene mapping audits. No README, plan or stated objective |
| 2 Nov 2026 | Repository name `coad_transportability`; `config/analysis_config.yaml` line 2: "COAD transportability — INFRASTRUCTURE configuration" | The project's working title names **transportability** as its aim |
| 2 Nov 2026 | `docs/PHASE0_CHANGELOG.md` | "Whether those values support **cross-cohort transport** is a Phase 4 question." Phase 0 separated computational fidelity from biological or normalisation validity: "A pipeline can pass every Phase 0 test and still be scientifically invalid." |
| 2 Nov 2026 | `docs/PHASE1_COHORT_ENDPOINT_INVESTIGATION.md` (rev. 4) | First approved scientific decisions: RFS definition (A1), stage II–III cohort (A2), S3 primary cohort (A3), "CIT subset = internal validation only, never external" (A6), chemotherapy as a prognostic covariate only (A7). Two endpoint assumptions were flagged, never verified: deaths censored (P1) and months as the time unit (P2) |
| 2 Nov 2026 | `docs/PHASE2_GENE_REPRESENTATION.md` §4 | The representative-probe representation was valued partly because "it **transfers**" to any other GPL570 cohort, which shows external use was still the intent |
| 2 Nov 2026 | `docs/PHASE4_PLAN.md`, `PHASE4_PREACQUISITION.md` | "External Validation Feasibility". Candidate GPL570 cohorts were screened, a GSE14333 ↔ GSE17536 patient overlap was found, and acquisition was ruled **NO-GO** |
| 2 Nov 2026 | `docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md` §8 | "Do not pursue confirmatory external validation… Proceed with a single-cohort, internally-validated study… transportability explicitly unestablished." |
| 2–3 Nov 2026 | `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md` §1 | The **revised** question: "do transcriptomic features add prognostic information beyond routinely available clinical variables for recurrence in stage II–III colon cancer?" "No transportability claim will be made." |

### 1.2 Reconstructed original objectives (before Phase 5)

| # | Original objective | Source |
|---|---|---|
| O1 | Build a computationally faithful, reproducible foundation for GSE39582 | Phase 0 |
| O2 | Define a defensible recurrence endpoint and stage II–III cohort without altering source data | Phase 1 |
| O3 | Produce a gene-level expression representation that could transfer to other GPL570 cohorts | Phase 2, 2B |
| O4 | **Assess whether a GSE39582-derived recurrence model transports to independent cohorts** (external validation) | repository name; Phase 0; Phase 2 §4; Phase 4 |
| O5 | Place GSE39582 and external cohorts on a common preprocessing scale (fRMA) | `PHASE4_FRMA_PILOT.md` |

### 1.3 Objectives as revised at Phases 5–6, and later additions

| # | Objective | When | Source |
|---|---|---|---|
| R1 | Test whether pre-specified transcriptomic pathway scores add discrimination beyond a clinical baseline (cause-specific hazard and C-index; internal validation only) | Phase 6B, 2 Nov | `PHASE6B` §1 |
| R2 | Fix the pathway panel **outcome-blind**, before any outcome was examined | Phase 6C, frozen 2 Nov 21:13 UTC | `PHASE6C` header |
| R3 | Prevent leakage: fold-internal standardisation and scoring, nested tuning | Phase 6B §3.2.2; Phase 9 protocol | — |
| R4 | Fit an exploratory transcriptome-wide model, M3 | Phase 6B §3, Phase 9 protocol §4.4 | **never fitted** |
| R5 | Run pre-specified endpoint and cohort sensitivity analyses S-1 to S-7 | Phase 6A | only **S-7** executed |
| R6 | Report optimism-corrected discrimination and calibration | Phase 6B §3.2.3 | done, single layer |
| L1 | **Added later:** a nested-bootstrap interval for the corrected ΔC | added after the execution report found "no interval was computed for the corrected ΔC" (`PHASE9_EXECUTION_REPORT.md` §15) | `PHASE9_DELTAC_INTERVAL_PLAN.md` |
| L2 | **Added later:** a corrected re-run after the as-run nested driver was found to reuse the resample's random state in cross-validation | owner Decisions 1–2 | FR §6 |
| L3 | **Added later:** a frozen post-run verifier, an independent numerical audit, a θ̂ alignment audit, a checkpoint archive release, public-release and citation audits, and results figures | 6–8 Nov | FR §15; commit history |

### 1.4 Answers to the specific questions in the brief

| Question | What the records show |
|---|---|
| Why GSE39582? | **Not recorded.** The dataset choice predates every planning document (it is in the baseline commit). The records do document its relevant properties: 585 GPL570 arrays (566 tumours), deposited RMA + ComBat, stage annotation, RFS fields and a built-in discovery/CIT partition (`PHASE0`, `PHASE1`, `PHASE3_DATASET.md` §1). Whether it was chosen for size or availability cannot be recovered |
| Biological and clinical problem | Recurrence after surgery in stage II–III colon cancer, where stage alone stratifies risk imperfectly. Whether expression adds information beyond routine clinical variables (`PHASE6B` §1; FR §1) |
| What we expected to learn | No formal hypothesis with an expected effect size was recorded. The Phase 6C literature grading **expected EMT to be the strongest candidate** ("the only candidate with strong progression AND recurrence evidence") and **expected Wnt to discriminate poorly** (`PHASE6C` §4, §6). The 8 November citation audit later withdrew the Wnt reasoning as unsupported |
| Why pathways rather than genes or clinical only | Phase 5 simulation: transcriptome-wide selection at 136 events gives an **apparent C of 0.92 when the true C is 0.50**, so it is "not feasible". A pre-specified pathway-level model with about 10–50 fixed features is "potentially feasible" (`PHASE5` §6.2). Pathways were the largest molecular feature set the event count could support |
| Why these Hallmark sets | MSigDB Hallmark v2026.2.Hs. Ten candidates were reduced to one representative per correlation cluster, on biology and literature only, never outcome (`PHASE6C` §3–§6). Four-set primary (EMT, E2F, Wnt, p53; condition number 2.95, no pair with \|r\| > 0.5); five-set (+PI3K) as pre-specified sensitivity (`PHASE6C` §0.5) |
| What M1 and M2 were to demonstrate | M1: a fixed clinical benchmark, "built first and never re-tuned". M2: whether a shrunken pathway block improves on it. Binding caveat: ΔC "must not be described as an unbiased estimate of the molecular block's contribution" because shrinkage makes it conservative (`PHASE6B` §3.1) |
| Rigour and reproducibility intent | Present from Phase 0 (hash tracking, tests, fidelity vs validity) and explicit throughout: outcome-blind panel, fold-internal scoring, frozen inputs "none of these may be revisited after seeing model results" (`PHASE6B` §2) |
| Original deliverables and success criteria | **No project-level success criteria were recorded at the start.** Criteria appear phase by phase: the Phase 4 acquisition conditions (D5) and event thresholds (D3: < 100 events exploratory only), Phase 5 precision targets, and later the plan's EP-1 Monte Carlo criterion and hard failures H1–H6 (`PHASE9_DELTAC_INTERVAL_PLAN.md`) |

---

## 2. Objective-by-objective achievement assessment

| Original objective | Why it mattered | What we implemented | Evidence of completion | Actual outcome | Status |
|---|---|---|---|---|---|
| O1 Reproducible, faithful data foundation | Every later result depends on reading the deposit correctly | Shared SOFT parser, cache verified bit-identical, SHA-256 tracking | `PHASE0_CHANGELOG.md`: 21 tests; final suite 1,097 passed | Data read faithfully and reproducibly | **Achieved** |
| O2 Defensible cohort and endpoint | Endpoint errors silently invalidate survival models | S0–S3 scenarios, verbatim source columns, A-19 complete-case population (452 patients, 136 events, 316 censored) | `PHASE1` rev. 4; `PHASE3_DATASET.md` (57/57 tests); `PHASE6A` | Cohort defined. The time unit (P2/U1) and death handling (P1/U2) remain **undocumented** | **Partially achieved**: a defensible cohort, but an endpoint with two unresolved assumptions |
| O3 Transferable gene representation | Needed for external use | Median aggregation (primary) and representative-IQR (sensitivity) | `PHASE2`, `PHASE2B`; S-7 | Both built. The choice between them is a convention: "A real answer would require an independent reference measurement… None exists" (`PHASE2B`) | **Achieved** as engineering; transfer never used |
| O4 External transportability | The project's original purpose | Cohort screen, overlap test, event-count policy, precision simulation | `PHASE4_PREACQUISITION.md` D1–D5; `PHASE5` §8 | **Not feasible**: GSE14333 has 48 events and GSE33113 18, against about 338 required at C = 0.65 (7–40× short). GSE17536 excluded for patient overlap | **Not assessable with available data** |
| O5 Common preprocessing scale (fRMA) | Prerequisite for O4 | Pilot protocol with a hard identifier gate | `PHASE4_FRMA_PILOT.md` | Prepared, **never executed**; deprioritised once O4 failed | **Not achieved** (deliberately not pursued) |
| R1 Do pathway scores add discrimination beyond clinical covariates? | The revised scientific question | M1 vs M2 penalised Cox, optimism-corrected C, nested-bootstrap interval | FR §7–§8; verifier PASS 60/60 | ΔC_corrected **+0.0043**; all five 95 % intervals include zero (about −0.009 to +0.05) | **Achieved**: answered, and the answer is **no demonstrated improvement** |
| R2 Outcome-blind pre-specification | Prevents retrofitted "pre-specified" panels | Panel frozen before any outcome analysis; 22 freeze tests | `PHASE6C` header; `tests/test_phase6c_m2.py` | Held: the five-set E2F result was not promoted | **Achieved** |
| R3 Leakage control | Leakage inflates performance | Fold-internal gene and score standardisation; nested λ tuning | `PHASE6B` §3.2.2; `PHASE8D`: production scorer agrees with an independent re-implementation to 3×10⁻¹⁵; fold isolation holds under deletion and poisoning | Verified | **Achieved** |
| R4 Exploratory transcriptome-wide model (M3) | Hypothesis generation | Specified only | `PHASE9_EXECUTION_REPORT.md` §1: "M3 was **not** fitted" | Not run | **Not achieved** (never authorised) |
| R5 Endpoint and cohort sensitivity analyses S-1 to S-7 | Tests robustness to endpoint construct and exclusions | S-7 (representative-IQR) only | `PHASE9_S7_SENSITIVITY_REPORT.md`; FR §13.8 "No Fine–Gray or composite-endpoint analysis was performed" | S-7: ΔC unchanged. The records contain no execution of S-1 to S-6 | **Partially achieved** |
| R6 Optimism-corrected discrimination and calibration | Apparent performance is optimistic | Single-layer bootstrap, B = 500 × 5 seeds | FR §7.3, §10 | Corrected C and slopes reported; no calibration intervals; no absolute calibration (U1) | **Partially achieved** |
| L1 Interval for the corrected ΔC | The execution report left the headline quantity without uncertainty | Nested bootstrap, 5 runs × 750 outer × 100 inner | FR §8; `corrected_interval_results.json` | Five verified percentile intervals, all including zero | **Achieved** |
| Demonstrate pathway-specific associations | Biological credibility | Primary-model HRs, LOBO, stratified fit | FR §9, §12 | No primary-model interval excludes 1; coefficients unstable across batches | **Not achieved** (the scientific answer is null) |
| Establish clinical usefulness | The ultimate translational aim | Not attempted | FR §13.9 | No decision-curve, reclassification or absolute-risk analysis | **Not assessable with available data** (U1 blocks horizon-based analysis) |
| Demonstrate external generalisability | Same as O4 | — | `PHASE5` §8 | No external cohort | **Not assessable with available data** |

**Key distinctions.**
- **Building a predictive pipeline:** achieved.
- **Showing molecular information improves prediction:** not shown.
- **Showing biological associations:** not shown.
- **Establishing clinical usefulness:** not assessable.
- **Showing external generalisability:** not assessable.

---

## 3. What we actually discovered

All values come from the frozen primary execution (`reports/phase9/execution/phase9_execution_results.json`), the S-7 results (`reports/phase9/execution_s7/phase9_s7_results.json`) and the corrected interval results (`reports/phase9/corrected_rerun_report/corrected_interval_results.json`).

**Analysis population (A-19):**
- 452 patients, 136 recurrences, 316 censored.
- Discovery partition: 353 patients, 107 events. CIT partition: 99 patients, 29 events.
- One patient (GSM972262) was excluded for missing chemotherapy status.

### 3A. Clinical findings (M1, unpenalised Cox)

| Term | HR (95 % Wald CI) | Observation |
|---|---|---|
| Stage III vs II | **1.97 (1.30–2.98)** | The only clinical term whose interval excludes 1 |
| Male sex | 1.38 (0.97–1.96) | Includes 1 |
| Proximal location | 0.80 (0.56–1.16) | Includes 1 |
| Adjuvant chemotherapy | 1.17 (0.76–1.79) | Includes 1. **Never a treatment effect:** confounding by indication (`PHASE1` A7) |
| Age (restricted cubic spline, two terms) | 1.005 (0.98–1.03); 1.010 (0.98–1.04) | Includes 1 |

- **Discrimination:** apparent C 0.6242; optimism-corrected C **0.6016** (Monte Carlo SE 0.0006).
- **Calibration:** corrected slope **0.833**. Even six coefficients overfit measurably at 136 events.
- **Batch sensitivity of stage:** in M2 the stage HR is 1.94 (1.28–2.94). With the nine batches as strata it falls to **1.68 (1.07–2.65)**. Under leave-one-batch-out its coefficient ranges from +0.498 to +0.826 (FR §9, §12).

**Interpretation:**
- Stage is the dominant measured clinical prognostic factor. Discrimination is modest: the clinical model orders about 60 % of comparable patient pairs correctly, against 50 % by chance.
- This does **not** show that age, sex, location or chemotherapy are unimportant clinically. Their intervals are wide at 136 events.
- The chemotherapy coefficient cannot be read causally.

### 3B. Molecular findings (M2, ridge-penalised pathway block, λ = 0.351)

| Pathway | HR per training-fold SD (95 % Wald CI), primary median | S-7 representative-IQR | Leave-one-batch-out log-HR range (sign changes) |
|---|---|---|---|
| E2F targets | 0.906 (0.807–1.016) | 0.872 (0.758–1.003) | −0.249 to −0.014 (0) |
| EMT | 1.070 (0.951–1.202) | 1.084 (0.937–1.254) | +0.009 to +0.204 (0) |
| Wnt/β-catenin | 1.016 (0.906–1.138) | 1.011 (0.879–1.162) | −0.002 to +0.067 (1) |
| p53 pathway | 0.995 (0.882–1.123) | 0.993 (0.852–1.158) | −0.113 to +0.014 (2) |

- **Observation:** no pathway interval excludes 1 in the primary model, under either representation. **No pathway shows an independent (adjusted) association** with recurrence.
- **Observation:** E2F stays below 1 and EMT above 1 under S-7, batch stratification and every leave-one-batch-out refit. Their magnitudes vary several-fold.
- **Caveat:** ridge shrinkage makes the Wald intervals approximate and descriptive, not post-selection inference (FR §9).
- **Five-set sensitivity model** (λ = 0.152): E2F HR **0.862 (0.746–0.996)**, narrowly excluding 1.
- **Why it was not promoted:**
  1. The panel was fixed as sensitivity-only before any fit (`PHASE6C` §0.5).
  2. Adding PI3K introduced collinearity: p53 ↔ PI3K r = 0.636, condition number 2.9 → 7.1, and every molecular SE inflated by a factor of 1.20–1.34.
  3. The E2F coefficient moved from −0.099 to −0.149 in step with that collinearity.
  4. The model's corrected C (0.6066) does not differ from the primary model's (0.6059) beyond Monte Carlo error.
  5. Under S-7 its corrected C fell by 0.0025 and its collinearity worsened (condition number 9.55) (FR §9.1, §11).

### 3C. Predictive performance

| | M1 | M2 four-set (primary) | M2 five-set (sensitivity) |
|---|---:|---:|---:|
| Apparent C | 0.6242 | 0.6382 | 0.6413 |
| Optimism-corrected C | 0.6016 | 0.6059 | 0.6066 |
| Optimism | 0.023 | 0.032 | 0.035 |

- **Apparent vs corrected:** M2's apparent advantage of +0.0140 is mostly optimism. **ΔC_corrected = θ̂ = +0.004316** (plan v2.1 §3.1).
- **Corrected nested bootstrap:** five runs differing only in random seed, each with 750 event-stratified outer × 100 inner resamples. Frozen verifier PASS 60/60. 95 % percentile intervals:

  | Run | Interval |
  |---|---|
  | 0 | [−0.0086, +0.0519] |
  | 1 | [−0.0090, +0.0469] |
  | 2 | [−0.0093, +0.0484] |
  | 3 | [−0.0091, +0.0518] |
  | 4 | [−0.0091, +0.0494] |

- **Monte Carlo precision:** between-run endpoint SD 0.00025 (lower) and 0.00217 (upper); largest endpoint MC-SE 4.2 % of interval width (EP-1 met).
- **Meaning:**
  - The data are compatible with a loss of about 0.009 and a gain of about 0.05 in C.
  - The interval's width reflects **sampling uncertainty with 136 events**, not Monte Carlo noise.
  - A zero-including interval does **not** establish equivalence or the absence of a pathway contribution.
- **Diagnostics:** the bootstrap distribution sits above θ̂ (warnings A7/A8 in every run; median ΔC about +0.012, 20.2 % of pooled draws below zero). This is unexplained. Percentile, basic and bias-corrected intervals all include zero in all 15 run-by-method cases (FR §8).
- **Incremental predictive value: not demonstrated.**

### 3D. Robustness and sensitivity

| Analysis | Evidence | Interpretation |
|---|---|---|
| Median vs representative-IQR (S-7) | Corrected ΔC +0.0043 under both; coefficient directions agree; M2 corrected slope 0.862 vs 0.832 (FR §11) | The discrimination result is robust to probe summarisation. **Not replication**: same arrays, same patients |
| Four-set vs five-set | Corrected C 0.6059 vs 0.6066; only the five-set E2F interval excludes 1, under collinearity | The five-set model adds a parameter, not information (`PHASE6C` §0.5 predicted this structurally) |
| Batch-stratified fit | Molecular coefficients moved ≤ 0.027; stage 0.662 → 0.520; chemotherapy 0.137 → 0.029 | Clinical terms are more batch-sensitive than the molecular block in this fit |
| Leave-one-batch-out (9 refits) | Apparent C 0.614–0.659 (0.638 on all data); Wnt and p53 change sign | Coefficients depend on cohort composition. Batch size varies about 7-fold (14–103 patients), so batch identity and sample size are confounded |
| Batch–partition confounding | Cramér's V = 0.49; residual ComBat-batch association with PC3 (ε² = 0.061, BH q = 2.8×10⁻⁵) (FR §12) | Part of the score variation may be technical. Batch effects are characterised, **not resolved** |
| Calibration | Corrected slopes 0.833 (M1), 0.862 (M2), 0.836 (five-set) | Residual overfitting: predictions would need about 14–17 % shrinkage. Log-hazard slope only; single layer, no intervals |
| Informative censoring | 46 deaths without recurrence are censored. Those patients were older (median 75.0 vs 68.0), less often stage III (32.6 % vs 38.7 %) and less often treated (23.9 % vs 42.1 %) | The cause-specific hazard is a valid aetiological quantity but cannot become an absolute recurrence probability |
| Time unit | `rfs.delay` unit undocumented (P2/U1) | Only rank-based quantities (C, HR) are reportable; no horizon-based metrics |

---

## 4. Biological interpretation

The rules applied here are those of `docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md`:
- A Hallmark score is the mean training-standardised **mRNA** level of a fixed gene set in **bulk tissue**.
- It is not pathway activation, protein activity, signalling flux or mutation status, and no mechanism can be inferred from it.
- References R1–R8 are the verified sources listed in that audit.

| Pathway | This dataset (observation) | Literature context (hypothesis only) | What cannot be concluded |
|---|---|---|---|
| EMT | HR 1.070 (0.951–1.202), consistently > 1, never distinguishable from 1 | Mesenchymal/stromal expression characterises CMS4, which had worse relapse-free survival (R1). GSE39582 contributed to the CMS derivation, so this context is **not independent** of this cohort. Mesenchymal-subtype transcripts in bulk tumours are mostly stromal (R2, R3) | Whether the Hallmark EMT score reflects stroma or tumour-cell EMT; this was not tested |
| E2F targets | HR 0.906 (0.807–1.016), consistently < 1 | **Mixed evidence:** lower proliferation-signature expression in relapsed primaries (R4), but high Ki-67 associated with poorer survival in a meta-analysis (R5) | An association (the interval includes 1); any mechanism |
| Wnt/β-catenin | HR 1.016; changes sign under leave-one-batch-out | Wnt-pathway genes are genomically altered in about 93 % of tumours (R6), but genomic alteration is not activation, and Wnt/MYC expression varies between subtypes (CMS2; R1) | Any effect. The original "near-universal, so poor discriminator" rationale (`PHASE6C` §4) was withdrawn as unsupported |
| p53 pathway | HR 0.995; changes sign twice | TP53 is mutated in about 60 % of non-hypermutated tumours (R6); target-gene expression is not mutation status | TP53 status (not analysed) |
| PI3K/AKT/mTOR (sensitivity) | HR 1.043 (0.891–1.221); collinear with p53 | PIK3CA mutated in about 15 % of non-hypermutated tumours (R6) | Any effect; it destabilised the other estimates |

**Overall biological reading:**
- The E2F (< 1) and EMT (> 1) directions are **hypotheses consistent with part of the prior literature** (EMT/stroma: R1–R3; proliferation: mixed, R4 vs R5).
- The dataset demonstrates **neither an association nor a mechanism**, and no explanation is offered for the null and unstable Wnt and p53 results.
- One exploratory finding bears on interpretation: the dominant transcriptional axis in this cohort is **MMR/CIMP/CIN/BRAF/tumour location**, not stage or partition (`PHASE8D_CONSOLIDATED_EDA_REVIEW.md`, executive summary).
- So biologically meaningful variation exists, but the frozen pathway panel was not designed around that axis, and MMR status was not a model covariate.

---

## 5. Unexpected findings and lessons

| # | Finding | Evidence | Lesson |
|---|---|---|---|
| U-1 | **The original aim was infeasible, and this was found before any model was built.** Available external GPL570 cohorts had 7–40× too few events, and two candidate cohorts shared about 129 patients | `PHASE4_PREACQUISITION.md` D1, D3; `PHASE5` §8 | Patient-level overlap checks and precision-based sample-size assessment should precede "external validation". Total GEO sample counts badly overstate usable events |
| U-2 | **Transcriptome-wide selection at this sample size is self-deception.** A null simulation gave apparent C 0.92 at true C 0.50 | `PHASE5` §6.2 | This is why the molecular question was restricted to four pre-specified scores |
| U-3 | **A marginal apparent gain was mostly optimism.** Apparent ΔC +0.0140 vs corrected +0.0043 | FR §7.3 | A single train-test or apparent C-index would have overstated the molecular contribution about threefold |
| U-4 | **A knife-edge "significant" interval belonged to the wrong quantity.** The historical paired out-of-original ΔC interval was [+0.0002, +0.0240] (median) but [−0.0001, +0.0224] under S-7 | FR §8.1 | It mixes in-sample and out-of-sample information and is not the corrected ΔC's uncertainty. Choosing the estimand before looking prevented a spurious positive claim |
| U-5 | **Computational stability is not statistical stability.** Monte Carlo checks passed with large margins and 0 of 7,500 fits failed, yet leave-one-batch-out flipped the signs of Wnt and p53 | `PHASE9_EXECUTION_REPORT.md` §15 | Reproducible numbers can still be fragile scientifically |
| U-6 | **Batch structure survives ComBat, and its effect depends on the analysis.** Phase 8D found batch location effects explained by MMR composition (0 of 10 score tests survived conditioning), but residual AFFX control-probe structure, PC2 dispersion and a PC3–batch association remain. In Phase 9, coefficients moved under batch omission and stage attenuated under stratification | `PHASE8D` summary; FR §12 | "Batch-corrected" deposits are not batch-free. Partition-confounded batches undermine the CIT partition as a test set |
| U-7 | **Pathway correlation shapes interpretation.** The ten-candidate panel had condition number 77 and six sets on one axis. Adding PI3K produced the only "significant" pathway result through collinearity | `PHASE6C` §3, §5; FR §9.1 | Per-pathway coefficients are interpretable only in well-conditioned panels. A significant coefficient in a collinear model is a warning, not a finding |
| U-8 | **The bootstrap distribution sat above θ̂** (A7/A8 in every run), opposite to the dilution mechanism the plan anticipated | FR §8 | Diagnostic intervals must be reported, not hidden. The qualitative result held across all interval types |
| U-9 | **A seed-reuse defect was found in the as-run nested driver.** Cross-validation reused the resample's random state. The whole campaign was re-run under corrected streams, and the as-run intervals were superseded (pooled range [−0.0091, +0.0510], close to the corrected ones) | FR §6, §14 (owner Decisions 1–2) | Random-stream independence must be designed and tested. Corrected and as-run results agreed, but only the corrected ones are reportable |
| U-10 | **A pre-execution gate could be bypassed through a shared cache.** An adversarial test found it; it was remediated before execution | `PHASE9_DECISION13_FREEZE_PROCEDURE.md` (gate-cache remediation) | Safety checks need adversarial tests too |
| U-11 | **Endpoint ambiguity limited every claim.** The time unit is undocumented, and 46 competing deaths are censored informatively | `PHASE1` P1/P2; `PHASE6A` | Rank-based discrimination only; no absolute risk, cumulative incidence, time-dependent AUC or clinical-utility metric |
| U-12 | **Eight arrays show expression-inferred sex discordant with recorded sex**, across five batches and consistent across 10–11 of 11 Y/XIST markers; 3 of 8 carry a QC soft flag (Fisher OR 9.67, p = 0.011) | `PHASE8B` §7 | A positive-control check exposed possible sample-metadata issues. They were not excluded or recoded, and this is a per-sample inference that needs data-governance care |
| U-13 | **The independent audit's scope is narrower than its label.** It recomputed 138 reported quantities from recorded per-fit values (agreement within 2.2×10⁻¹⁶) and confirmed the records were complete and hash-linked. It did **not** refit models; fit correctness rests on the verifier's bit-identical recomputation of 20 replicates (V5) | `INDEPENDENT_AUDIT_REPORT.md` line 77 | It establishes that the reported statistics follow exactly from the records. It does not establish external replication or model correctness beyond V5 |
| U-14 | **Biological rationale written at panel selection needed later correction.** The 8 November citation audit narrowed or withdrew statements on proliferation, Wnt and EMT/stroma | `BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md` §4 | Panel rationale should be checked against verified sources at the time of freezing |

---

## 6. What we built

| Component | What it is | Why it is technically meaningful | Evidence |
|---|---|---|---|
| Faithful data layer | Single SOFT parser; Parquet cache verified bit-identical to the deposit; SHA-256 raw-data tracking | Separates "we read the data correctly" from "the data are valid", which most pipelines conflate | `PHASE0_CHANGELOG.md`; 21 tests |
| Cohort and endpoint construction | Verbatim source columns plus derived columns; superset storage (S3 ⊂ S1 ⊂ S0); `rfs_*_60units` naming refuses to assert the unit | Uncertainty is visible in the artifact itself; no exclusion is silent | `PHASE3_DATASET.md`; 57/57 tests at Phase 3 |
| Probe-to-gene representations | Median (43,163 probes → 21,815 genes) and representative-IQR | Two frozen, auditable representations allow a pre-specified sensitivity analysis (S-7) | `PHASE2`, `PHASE2B` |
| QC protocol | Soft flags on median, IQR and correlation fences; 29 of 453 flagged, none excluded without technical evidence | Avoids outcome-adjacent array removal | `PHASE7B`, `PHASE7C` |
| Outcome-blind pathway panel | Correlation-cluster representatives; frozen JSON specs with hashes; 22 freeze tests | Pre-specification is provable, not asserted | `PHASE6C`; `tests/test_phase6c_m2.py` |
| Leakage-controlled scoring | Per-fold gene standardisation (ddof pinned); scoring and λ tuning inside each training fold | Production scorer matches an independent re-implementation to 3×10⁻¹⁵; fold isolation holds under deletion and poisoning | `PHASE8D` executive summary |
| Penalised Cox design | Unpenalised clinical block (scale-equivariant); ridge only on pathway scores; M1 nested in M2 | Makes ΔC a well-defined, conservative comparison | `PHASE6B` §3.1.1 (P6B-A1) |
| Patient-level nested bootstrap | 5 runs × 750 outer × 100 inner, event-stratified; paired M1/M2 on the same resamples; about 765,000 model fits plus 60 CV fits per M2 fit | Puts a valid interval on an optimism-corrected quantity, which a single bootstrap cannot | FR §6, §15 |
| Deterministic random-stream architecture | Separate NumPy `SeedSequence` streams for outer, inner, inner-CV and apparent-CV | Removes the stream reuse that invalidated the as-run campaign; makes every replicate reproducible | FR §6 |
| Checkpointing and recovery | Append-only, hash-linked, crash-consistent checkpoints; Path C top-up from 500 to 750 when EP-1 required it | Allowed multi-hour runs to be resumed and extended without re-running or altering completed records | FR §6, §15; `scripts/45` |
| Frozen protocol and verification | 27-file code freeze (`618f20eb…a790`); configuration hash bound into every record; pre-execution gate; frozen post-run verifier | Any change to code or configuration after freezing is detectable; reporting was gated on verifier PASS | Freeze record; `L00750.report.json` PASS 60/60 |
| Independent numerical audit | Standard-library-only re-implementation (`scripts/48`) | Rebuilt 3,750 outer replicates and 375,000 inner pairs from raw checkpoint bytes; 138 quantities within 2.2×10⁻¹⁶. Scope as in U-13 | `reports/phase9/independent_audit/` |
| θ̂ alignment audit | Showed the point estimate is the same functional the interval bootstraps | Closed a potential estimand mismatch; bounded one unmeasured RNG link analytically (worst case 1.6×10⁻³) | `PHASE9_THETA_ALIGNMENT_AUDIT.md` |
| Automated tests | 1,097 tests, all passing | Cover freezes, leakage, endpoint construction, verifier logic and figure traceability | full suite, 8 Nov 2026 |
| Results figures | 9 main + 2 supplementary figures, generated deterministically from frozen outputs (`scripts/49`); write-once output; `models_fitted: 0` | Every plotted value is traceable to a hashed source; 12 dedicated tests | `reports/figures/final/`; `docs/SCIENTIFIC_RESULTS_AND_FIGURES.md` |
| Version control and archive | Private GitHub repository; release `phase9-corrected-campaign-2026-11-08` with a 579,716,096-byte checkpoint archive (SHA-256 `0812b3ba…2aa2f`); 52-file manifest | Checkpoints too large for Git are preserved verifiably outside history | `docs/LARGE_ARTIFACT_ARCHIVE.md` |

**Verification evidence (latest, 8 November 2026):**
- full suite 1,097 passed;
- frozen code 27 of 27 hashes match;
- archive manifest 52 of 52 OK;
- verifier PASS 60/60;
- independent audit: 138 quantities within 2.2×10⁻¹⁶;
- Phase 8D cross-phase corroboration 15 of 15.

**Portability limitation:** the frozen launcher and verifier compare the interpreter path and virtual-environment prefix, so they are not portable to another machine. The standard-library audit (`scripts/48`) is portable (`PUBLIC_RELEASE_REPRODUCIBILITY_AUDIT.md`).

---

## 7. What we did not achieve

### 7.1 Scientific limitations: not achievable with this dataset and design

| Goal | Why not | Source |
|---|---|---|
| External validation | No external cohort with enough events; GSE17536 overlaps GSE14333 | `PHASE5` §8 |
| Independent-cohort replication | Same reason; S-7 and the CIT partition are not replication | FR §11, §13 |
| Transportability (the original aim) | Requires external data | `PHASE5` §6.2 |
| Absolute recurrence risk | Time unit undocumented (U1); competing deaths (U2) | FR §13.6 |
| Competing-risk cumulative incidence | S-1 (Fine–Gray) specified but not run; would also need U1 for horizons | `PHASE6A`; FR §13.8 |
| Clinical utility | Needs absolute risk at a horizon; no decision-curve or reclassification analysis | FR §13.9 |
| Causal biological mechanism | Bulk mRNA scores and an observational design cannot support it | FR §9.2 |
| Reliable pathway-specific prognostic associations | No primary interval excludes 1; signs unstable across batches | FR §9, §12 |
| Resolution of batch sensitivity | Characterised only; the batch-perturbation analysis was not performed | FR §12 |
| Resolution of the time unit | No source documents it | `PHASE1` P2 |
| Generalisation to other platforms or populations | GPL570 data processed like GSE39582 only | `PHASE9_SCIENTIFIC_INTERPRETATION_REVIEW.md` |

### 7.2 Planned analyses that were not carried out

- M3, the exploratory transcriptome-wide model.
- Endpoint and cohort sensitivity analyses S-1 to S-6.
- The fRMA pilot.
- The exploratory external direction-of-effect check (`PHASE5` §8.3).
- The batch-perturbation sensitivity analysis.
- A likelihood-ratio test of M2 vs M1, listed in `PHASE6B` §3.2.3; no record of it appears in the Phase 9 results.
- Calibration, coefficients, S-7 and batch analyses were not re-estimated under the corrected random streams. Only ΔC has a nested interval.

### 7.3 Unfinished engineering and publication tasks

- Code licence not selected.
- Repository private. Public release pending the data-governance decisions on per-sample tables (D1, D7 in `PUBLIC_RELEASE_REPRODUCIBILITY_AUDIT.md`) and remediation items F1–F15.
- The public-release audit is uncommitted.
- Frozen launcher and verifier not portable.
- No external DOI deposit.

---

## 8. Overall assessment

1. **Did we answer our original research question?**
   - The **original** question, transportability, was **not answered**: Phases 4–5 showed it cannot be answered with available public data.
   - The **revised** question, whether four pre-specified pathway scores add discrimination beyond clinical covariates, **was answered**: no demonstrated improvement (ΔC_corrected +0.0043; all five intervals include zero).

2. **Strongest positive achievement.** A null result that can be trusted:
   - the panel was fixed outcome-blind;
   - leakage was controlled and independently checked;
   - the headline quantity has a verified nested-bootstrap interval;
   - every reported number is reproducible from hashed records.

3. **Most important scientific finding.** In this cohort, the apparent molecular gain (+0.014 in C) is mostly optimism. After correction, the gain is +0.004, with uncertainty spanning a small loss to a modest gain.

4. **Most important null finding.** No pathway in the primary model shows an adjusted association with recurrence. The only "significant" pathway result (five-set E2F) arises in a collinear sensitivity model and was correctly not promoted.

5. **Most technically sophisticated achievement.** The corrected nested-bootstrap campaign:
   - per-replicate deterministic random streams;
   - paired inner optimism estimation;
   - hash-linked checkpoints with a rule-driven top-up;
   - frozen code and configuration binding;
   - a gated verifier;
   - a separate standard-library audit.

6. **What a conventional train-test split would have missed:**
   - the threefold optimism in the apparent gain;
   - that a narrowly positive interval belonged to a different quantity;
   - that the interval's width is driven by 136 events;
   - that coefficients flip sign when one batch is removed;
   - that the clinical model itself overfits (slope 0.833);
   - that the CIT "validation" partition is confounded with batch and is not a test set.

7. **What can be claimed in a portfolio:**
   - designed and executed a pre-specified, leakage-controlled penalised-Cox evaluation of pathway scores on 452 patients (136 events);
   - built a reproducible nested-bootstrap pipeline (about 765,000 model fits) with frozen code, a verifier and an independent numerical audit;
   - found that the pathway scores did not demonstrably improve optimism-corrected discrimination;
   - identified and documented batch instability, informative censoring and endpoint limitations;
   - showed that external validation was infeasible before building models (precision assessment, patient-overlap detection).

8. **What must never be claimed:**
   - that the pathways improve prediction, or are prognostic biomarkers;
   - statistical significance of any gain;
   - equivalence or "no effect";
   - clinical utility, absolute risk or patient stratification;
   - external validation, transportability or generalisability;
   - a causal or mechanistic biological finding;
   - that the five runs are independent experiments or cohorts;
   - that the independent audit is external replication or model refitting;
   - that the CIT partition is external validation;
   - that the five-set E2F result is a finding.

9. **Logical next research project.** Two directions follow from the gaps:
   - **External test, if data allow:** obtain a cohort with about 340 or more stage II–III recurrence events, documented time units and death handling (`PHASE5` §8.5); preprocess uniformly (fRMA or raw CEL); test the frozen M1/M2 unchanged.
   - **Within-cohort follow-up (new, pre-registered analyses):**
     - stromal-content deconvolution to test whether the Hallmark EMT score is stroma-driven;
     - models that include MMR/CIMP/BRAF, the dominant expression axis;
     - competing-risk analysis, if the time unit can be documented;
     - the batch-perturbation analysis.

**Overall verdict.**
- **As a test of molecular added value:** a rigorous, honest study that answered its revised question with a null result and documented exactly why stronger claims cannot be made.
- **As the transportability study it set out to be:** it did not succeed, and the records show this was established early and deliberately, not hidden.
- Calling the project "successful" in the sense of a positive predictive finding would be wrong. Its value lies in the reliability of its negative answer and the infrastructure that makes that answer verifiable.

---

## 9. Portfolio translation

### A. One-sentence achievement
I built a pre-specified, leakage-controlled penalised-Cox pipeline with a verified nested bootstrap (about 765,000 model fits) to test whether four Hallmark pathway scores improve recurrence-risk discrimination beyond clinical covariates in 452 stage II–III colon cancers (GSE39582), and found no demonstrated improvement: optimism-corrected ΔC = +0.004, with 95 % intervals of about −0.009 to +0.05.

### B. Resume entry
**Colon-cancer recurrence modelling: optimism-corrected value of transcriptomic pathway scores** (Python, lifelines, NumPy)
- Built a reproducible survival-analysis pipeline on GEO GSE39582 (585 Affymetrix arrays; 452 stage II–III patients, 136 recurrences): probe-to-gene mapping, QC, cohort and endpoint audit, and a four-pathway MSigDB Hallmark panel frozen before any outcome analysis.
- Compared a clinical Cox model with a ridge-penalised pathway extension using fold-internal scoring and tuning, and quantified the optimism-corrected gain in Harrell's C with a nested bootstrap (5 seed-replicate runs × 750 × 100 resamples, about 765,000 fits), hash-frozen code, checkpointing and a post-run verifier (60/60 checks).
- Result: no demonstrated improvement (ΔC +0.004; all intervals include zero). Apparent gain was about threefold optimism; coefficients were unstable across processing batches. Showed beforehand that external validation was infeasible (7–40× event shortfall; patient overlap between public cohorts).

### C. Interview explanation (about 4 minutes)

**The problem.** Stage II–III colon cancer recurs in a substantial minority of patients after surgery, and stage alone stratifies risk imperfectly. Many papers report that expression signatures improve prognosis. I wanted to know whether that holds up once you correct for optimism, the tendency of a bigger model to fit its own data better than new patients.

**Why this dataset and approach.** GSE39582 is a large public Affymetrix colon-cancer cohort with recurrence data. I originally wanted to test transportability to other cohorts. Before building any model, I checked feasibility. The available external cohorts had 48 and 18 recurrence events, against roughly 340 needed for precise validation, and two cohorts turned out to share about 130 patients. So I narrowed the project to a rigorously internally validated question. A simulation also showed that selecting genes transcriptome-wide at 136 events produces an apparent C of 0.92 when the truth is 0.5. So I pre-specified four Hallmark pathway scores (EMT, E2F targets, Wnt, p53), chosen outcome-blind from correlation structure and literature.

**What I built.**
- A clinical Cox model compared with the same model plus ridge-penalised pathway scores.
- Every scoring and tuning step happens inside each training fold, so there is no leakage.
- For uncertainty, a nested bootstrap: for each of 750 resamples I re-estimated the optimism correction with 100 inner resamples, repeated in five runs differing only in random seed. That is about 765,000 model fits.
- Engineering controls: deterministic per-replicate random streams, hash-linked checkpoints, frozen code, a verifier that had to pass before any result could be reported, and a separate re-implementation that recomputed every reported number from the raw checkpoints.

**What I found.**
- The pathway model's apparent C advantage of 0.014 was mostly optimism; the corrected gain was 0.004.
- Every interval included zero, ranging from a small loss to a gain of about 0.05.
- No pathway hazard ratio was distinguishable from 1 in the primary model.

**What surprised me.**
- An earlier, different bootstrap quantity had an interval that just excluded zero. It would have been easy to report that as a positive result, but it wasn't the right estimand.
- Leaving out one processing batch at a time flipped the sign of two pathway coefficients, even though the batches were "ComBat-corrected".
- I also found and fixed a random-seed reuse bug in my own first nested run, and re-ran the whole campaign.

**Limitations.**
- Internal validation only.
- The recurrence time unit is undocumented, so no absolute risk or clinical-utility metrics.
- Deaths are censored, and informatively.
- Batch effects are characterised, not resolved.

**What I learned.**
- Validation design matters more than model complexity.
- Separating the estimand from convenient-looking statistics prevents false positives.
- Computational reproducibility and statistical stability are different things.
- A well-verified null result is a legitimate scientific outcome.

---

## 10. Evidence sources and unrecoverable information

### 10.1 Sources used

**Git history:** 16 commits on `phase9-corrected-pipeline`, from the root `1cd2bba` (Phase 0, including the baseline scripts) to `51b283e`; and the local-only `main` commits `6fda652` (baseline) and `1b8c1c2`, inspected read-only.

**Early phases:**
- `docs/PHASE0_CHANGELOG.md`
- `docs/PHASE1_COHORT_ENDPOINT_INVESTIGATION.md` (status summary A1–A10, P1–P3, U1–U7)
- `docs/PHASE2_GENE_REPRESENTATION.md` §4
- `docs/PHASE2B_SELECTION_SENSITIVITY.md`
- `docs/PHASE3_PLAN.md`, `docs/PHASE3_DATASET.md` §1
- `docs/PHASE4_PLAN.md`, `docs/PHASE4_PREACQUISITION.md` §0 (D1–D5), `docs/PHASE4_FRMA_PILOT.md`
- `docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md` §6–§8
- `docs/PHASE6A_ENDPOINT_RESOLUTION.md` (S-1 to S-7)
- `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md` §1–§3
- `docs/PHASE6C_PATHWAY_PANEL_REVIEW.md` §0.5, §3–§6
- `docs/PHASE7B_QC_PROTOCOL.md`, `docs/PHASE7C_QC_REPORT.md`
- `docs/PHASE8B_PCA_BATCH_INVESTIGATION.md` §7
- `docs/PHASE8D_CONSOLIDATED_EDA_REVIEW.md` executive summary

**Phase 9:**
- `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md`
- `docs/PHASE9_EXECUTION_REPORT.md` §1, §15
- `docs/PHASE9_S7_SENSITIVITY_REPORT.md`
- `docs/PHASE9_DELTAC_INTERVAL_PLAN.md`
- `docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md`
- `docs/PHASE9_THETA_ALIGNMENT_AUDIT.md`
- `docs/PHASE9_CORRECTED_FINAL_REPORT.md` §1, §6–§15
- `docs/OWNER_DECISION_REGISTER.md`

**Frozen results:**
- `reports/phase9/execution/phase9_execution_results.json` (apparent, monte_carlo, delta_C, batch_sensitivity)
- `reports/phase9/execution_s7/phase9_s7_results.json`
- `reports/phase9/corrected_rerun_report/corrected_interval_results.json`
- `reports/phase9/corrected_rerun/verification/L00750.report.json`
- `reports/phase9/independent_audit/INDEPENDENT_AUDIT_REPORT.md`

**Figures and later documentation:**
- `reports/figures/final/FIGURE_CAPTIONS.md`
- `docs/SCIENTIFIC_RESULTS_AND_FIGURES.md`
- `docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md`
- `docs/PUBLIC_RELEASE_REPRODUCIBILITY_AUDIT.md`
- `docs/FINAL_PROJECT_HANDOVER.md`

### 10.2 Information that could not be recovered

1. **Why GSE39582 was chosen.** The choice predates every document; the baseline commit has no README or plan.
2. **A written statement of the original objectives.** The transportability aim is inferred from the repository name, the Phase 0 configuration and changelog, Phase 2 §4 and the Phase 4 framing. No project charter exists.
3. **Original success criteria or expected effect sizes.** None were recorded at the start; criteria appear only phase by phase.
4. **The rationale for not fitting M3 or running S-1 to S-6.** The records state they were not authorised or performed, but give no explicit decision explaining why.
5. **Whether a likelihood-ratio test was ever computed.** It is listed in `PHASE6B` §3.2.3, but there is no record in the Phase 9 results.
6. **The primary sources behind several Phase 6C literature grades.** These were written as summaries (for example "recurrence-signature literature (§8)"). Only the subset re-checked on 8 November 2026 is verified (R1–R8).

### 10.3 Integrity statement

- This retrospective adds one new file, `docs/PROJECT_RETROSPECTIVE_AND_ACHIEVEMENTS.md`.
- No model was fitted, no analysis was run, and no number was recomputed; every value is quoted from the records above.
- No model code, frozen output, checkpoint, configuration, figure, archive or existing report was modified.
- Nothing was committed, pushed or published.
