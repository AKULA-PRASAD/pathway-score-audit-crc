# Additional analyses: methods and results (registered exploratory work, X1–X9)

**Date:** 8 November 2026.
**Plan:** `docs/EXPLORATORY_ANALYSIS_PLAN.md`, registered at 2026-11-08T21:29:20Z before any of these analyses ran. Its SHA-256 (`2943119e…`) is recorded in `reports/exploratory/plan_registration.json` and re-checked by `x_validate.py`.
**Status of every result here:** post-hoc exploratory or sensitivity work, using **the same data as the primary analysis**, carried out **after** the primary results were known. No result here is confirmatory. An independent re-implementation (§Independent review) reproduced the key quantities.

**Nothing here changes the frozen primary analysis.** The frozen primary results are:
- M1 corrected C 0.6016;
- M2 four-set corrected C 0.6059;
- corrected ΔC +0.004316, with all five 95 % intervals including zero;
- calibration slopes 0.833 / 0.862 / 0.836;
- every primary pathway HR interval includes 1.

**Reproduce:** `bash scripts/exploratory/run_all.sh`. This needs the project `.venv`, plus Rscript with base `survival`. Outputs go to `reports/exploratory/`, and the validation checks to `reports/exploratory/validation_checks.json` (22/22 PASS).

## Common inputs

| Item | Value |
|---|---|
| Population | Frozen A-19 analysis cohort: **452 patients, 136 recurrences**, loaded with `src/phase9_corrected/data.py::load_context` |
| Pathway scores | Frozen four-set definitions and frozen scorer (`src/pathway_scoring.score_panel`), standardised on all 452 patients. Check: on the 453-patient S3 cohort they reproduce `reports/phase8c/scores_4set_primary_reproduced_full_cohort.csv` to 2.2 × 10⁻¹⁶ |
| Clinical design | Frozen M1 design (`src/phase9_execution.build_clinical`) |
| Endpoint | `rfs_time_full` / `rfs_event_full`; time unit **provisional: months** |
| Downloads | `data/external/exploratory/` with `PROVENANCE.json` (URL, release, date, SHA-256, licence, purpose) |

## X1. Molecular anchoring of the scores (outcome-free)

- **Question:** do the expression scores track the genomic and molecular features expected from their biology?
- **Method:**
  - Mann–Whitney U, reporting the rank-biserial correlation and the difference in medians in score-SD units.
  - Kruskal–Wallis with ε² for the six CIT subtypes.
  - BH-FDR across 24 tests.
  - Pre-specified secondary test: TP53 within pMMR tumours, using a batch-stratified permutation (10,000).
- **Data:** `scripts/exploratory/x1_molecular_anchoring.py` → `reports/exploratory/x1_molecular_anchoring/`.
- **Missingness** (patients with a valid value, of 452): TP53 286, KRAS 431, BRAF 398, CIMP 385, CIN 363, CIT subtype 452.

| Feature (n present / absent) | EMT | E2F | Wnt | p53 |
|---|---|---|---|---|
| TP53 mutant (156 / 130) | −0.06 | **+0.19** (q 0.011) | +0.14 (q 0.074) | **−0.30** (Δmedian −0.45 SD, q 4.1 × 10⁻⁵) |
| KRAS mutant (161 / 270) | 0.00 | **−0.13** (q 0.044) | +0.03 | +0.11 |
| BRAF mutant (41 / 357) | **+0.24** (q 0.022) | +0.17 | **−0.46** (Δmedian −0.74 SD, q 4.8 × 10⁻⁶) | **+0.48** (Δmedian +0.77 SD, q 1.5 × 10⁻⁶) |
| CIMP+ (75 / 310) | +0.14 | **+0.20** (q 0.011) | **−0.37** (q 2.2 × 10⁻⁶) | **+0.40** (q 2.5 × 10⁻⁷) |
| CIN+ (282 / 81) | −0.04 | −0.06 | **+0.23** (q 0.0033) | **−0.44** (q 1.4 × 10⁻⁸) |
| CIT subtype (6 groups, ε²) | **0.42** | **0.21** | **0.08** | **0.16** (all q < 10⁻⁶) |

Rank-biserial effects; bold = BH q < 0.05 across the registered 24 tests. Because the CIT subtypes are circular (see below), BH is also reported across the 20 non-CIT tests (deviation 7). Every association above keeps q < 0.05 except KRAS–E2F, which becomes q = 0.050.

**Secondary test.** Within pMMR tumours (147 TP53-mutant, 106 wild-type), the p53 score is lower in mutants by 0.39 SD (batch-stratified permutation p = 0.0003). The TP53 association is therefore not an artefact of MMR status: dMMR tumours are mostly TP53 wild-type (22 of 29 dMMR tumours with TP53 status).

**Interpretation:**
- The p53 *target-gene expression* score is lower in TP53-mutant tumours, as expected if missense mutations impair p53 transactivation.
- The Wnt score is lower, and the p53 score higher, in BRAF-mutant and CIMP+ tumours. This is consistent with the serrated, TP53-wild-type, lower-WNT-target route described in the CRC literature. It is an association, not a mechanism.
- EMT differs most across CIT subtypes, highest in C4 (median +1.68 SD).
- **Caveats:**
  - CIT subtypes were *derived from the same expression data*, so the subtype associations are partly circular and must not be read as independent validation.
  - BRAF mutation, CIMP+, CIN− and dMMR form one correlated (serrated, right-sided, often MSI) cluster. The repeated Wnt-down / p53-up pattern across BRAF, CIMP and CIN is **one finding, not three**.
  - The higher p53 score in that cluster may partly reflect immune or interferon-responsive genes within the Hallmark p53 set: p53 correlates 0.32 with ImmuneScore and 0.36 with cytotoxic lymphocytes (X6).
  - The TP53 result is the cleanest anchor (it survives within pMMR). Even so, it cannot separate p53 transactivation from CIN or tumour location.
  - None of these results is evidence of pathway function or activity.

## X2. Competing risks

**Event coding (gate passed):**
- 136 recurrences;
- 46 deaths without recurrence (`rfs_event_full = 0`, `os_event = 1`, equal times), reproducing Phase 6A;
- 270 censored.

**Aalen–Johansen cumulative incidence** (95 % CI; time unit provisional):

| Group | Recurrence at 36 | Recurrence at 60 | Death without recurrence at 60 |
|---|---|---|---|
| All (452) | 26.6 % (22.7–31.1) | 30.7 % (26.6–35.5) | 6.6 % (4.6–9.6) |
| Stage II | 17.8 % (13.6–23.3) | 22.9 % (18.0–29.0) | 7.8 % (4.9–12.3) |
| Stage III | 37.5 % (31.2–45.1) | 40.5 % (34.0–48.2) | 5.2 % (2.7–9.8) |

Estimator: R `survival::survfit` multistate Aalen–Johansen, which is tie-aware and deterministic (deviation 5). For comparison, the naive 1 − Kaplan–Meier estimate with deaths censored gives 31.3 % at 60. Treating deaths as censoring overstates recurrence probability only slightly (0.6 percentage points at 60) because the competing event is uncommon. Death without recurrence is more frequent in stage II, where patients are older and less often treated.

**Fine–Gray subdistribution HRs vs cause-specific HRs** (recurrence; clinical + four scores; unpenalised):

| Term | Cause-specific HR (95 % CI) | Fine–Gray sHR (95 % CI) | Frozen ridge M2 HR |
|---|---|---|---|
| Stage III | 1.92 (1.26–2.92) | 1.91 (1.23–2.96) | — |
| EMT | 1.12 (0.93–1.35) | 1.12 (0.94–1.34) | 1.070 (0.951–1.202) |
| E2F | 0.82 (0.69–0.98) | 0.83 (0.70–0.98) | 0.906 (0.807–1.016) |
| Wnt | 1.04 (0.88–1.23) | 1.05 (0.89–1.22) | 1.016 (0.906–1.138) |
| p53 | 0.99 (0.81–1.22) | 0.98 (0.79–1.21) | 0.995 (0.882–1.123) |

**Cross-check:** the unpenalised clinical model reproduces the frozen M1 stage HR, 1.966 (1.296–2.982).

**Interpretation:**
- Subdistribution and cause-specific estimates agree closely, so the competing-risk structure does not change any pathway interpretation.
- The scores show no association with death without recurrence (all intervals include 1). Those models have 46 events for 10 terms (about 4.6 per term), so they are unstable and only descriptive.
- In the **unpenalised** models, the E2F interval excludes 1. This does **not** overturn the primary result, for four reasons:
  - The registered primary model is ridge-penalised; its E2F HR is 0.906 (0.807–1.016). The unpenalised estimate is the same data without shrinkage.
  - It is one of four unadjusted coefficient tests; the Bonferroni-adjusted p is about 0.11.
  - Adding dMMR to the model (409 patients with MMR status) moves it to 0.871 (0.724–1.049) (X7, deviation 8). dMMR tumours have higher E2F scores and a favourable prognosis.
  - The primary estimand is discrimination, and the corrected ΔC includes zero.

## X3. Batch estimability and frailty sensitivity

**Batch structure** (9 batches):

| Batch | n | Recurrences | Deaths w/o recurrence | Stage III % | Chemo % | CIT partition n |
|---|---|---|---|---|---|---|
| batch01 | 33 | 8 | 3 | 42 | 9 | 0 |
| batch02 | 14 | 6 | 0 | 100 | 100 | 0 |
| batch03 | 102 | 24 | 13 | 26 | 41 | 8 |
| batch04 | 103 | 52 | 2 | 57 | 59 | 27 |
| batch05 | 48 | 6 | 6 | 35 | 31 | 4 |
| batch07 | 81 | 13 | 14 | 37 | 37 | 37 |
| batch08 | 41 | 16 | 3 | 66 | 71 | 4 |
| batch09 | 15 | 4 | 2 | 40 | 20 | 14 |
| batch10 | 15 | 7 | 3 | 40 | 20 | 5 |

Batches differ strongly in case mix: stage III ranges from 26 % to 100 %, and batch 4 holds 52 of 136 recurrences.

**Leave-one-batch-out diagnostics** (frozen ridge fits, read only):
- 34 of 36 coefficient changes lie within ±1.96 full-cohort SE.
- The two exceptions are EMT (+2.3 SE) and E2F (−2.5 SE). Neither changes sign.
- The sign changes (Wnt 1, p53 2) occur where the full-cohort coefficient is within 0.3 SE of zero (Wnt β 0.015, SE 0.058; p53 β −0.005, SE 0.062). They reflect near-zero estimates, not instability of a real effect.
- Only for EMT does the size of the change track the omitted batch's events (Spearman 0.73, p = 0.026).

**Frailty sensitivity models** (unpenalised clinical + four scores):

| Model | Batch variance | EMT HR | E2F HR | Wnt HR | p53 HR |
|---|---|---|---|---|---|
| No frailty | — | 1.12 | 0.82 | 1.04 | 0.99 |
| Gamma frailty | θ̂ = 0.13; profile 95 % ≈ 0.03–0.50 (refined grid, deviation 6); LRT 11.1, boundary p = 0.0004 | 1.14 | 0.82 | 1.02 | 1.03 |
| Gaussian frailty | 0.14 | 1.14 | 0.82 | 1.02 | 1.03 |

**Interpretation:**
- There is clear between-batch heterogeneity in the **baseline** recurrence hazard after clinical adjustment. It is estimable with nine batches, but θ is imprecise (interval spans roughly an order of magnitude).
- Allowing for it leaves every pathway HR essentially unchanged.
- Batch-score location associations stratified by `combat_mod` (Phase 8C) were not significant (p 0.38–0.98).
- **Classification by the pre-registered rule:** not "resolved", because θ's upper bound is above 0.05 and two LOBO changes exceed 1.96 SE. It meets "partially mitigated": the frailty model is estimable and bounded, and coefficients are stable under it.
- The source of the baseline heterogeneity cannot be identified: technical batch, centre, recruitment era, or case mix beyond the measured covariates. Batch is confounded with all of these.

## X4. MMR status and `combat_mod` provenance

**Provenance:**

| Source | MMR field | ComBat model field |
|---|---|---|
| GEO SOFT (`!Sample_characteristics_ch1`) | `mmr.status` | **absent**. Processing text: "RMA normalisation … + ComBat (R sva package)", with no model-matrix description |
| ArrayExpress SDRF (E-GEOD-39582) | `Characteristics [mmr.status]` | `Characteristics [normalizationcombatmod]` = modcan1 / modcan2 |
| Marisa 2013 Table S1 | `MMR.Status` | absent |

- MMR agrees 100 % across GEO, SDRF and Table S1 for the 452 patients.
- In the whole deposit:
  - modcan1 ↔ dMMR in 75/75 known samples, and modcan2 ↔ pMMR in 444/444;
  - **47 samples with MMR missing nevertheless have a `combat_mod` value** (42 modcan2, 5 modcan1).
- **What is now established:**
  - `combat_mod` is the ComBat `mod` covariate, i.e. the biological covariate protected during batch correction. This follows from the SDRF field name `normalizationcombatmod`.
  - It is perfectly concordant with MMR where MMR is recorded.
  - It is not simply a copy of the published MMR column, because it is defined where MMR is missing.
  - Its exact source (for example an internal MSI assay or an earlier classification) is **not documented**. It remains a separate variable and is not relabelled as MMR.
  - In the 452-patient cohort, the 43 MMR-missing patients were assigned to modcan2 (38) or modcan1 (5).
  - **Consequence:** because ComBat protected this covariate while removing batch effects, and the groups are unbalanced across batches, score differences between MMR groups in the deposited data may be *inflated* (Nygaard, Rødland & Hovig, *Biostatistics* 2016;17:29–39). Every score–MMR association, and any MMR-adjusted model, should be read with that in mind.

**Robustness of score–MMR associations** (n = 409 with MMR; 59 dMMR, 350 pMMR):

| Score | ε² | Batch-stratified permutation p (BH q) | ε² range leaving one batch out | Within CIT C2 (50 dMMR) ε², p |
|---|---|---|---|---|
| EMT | 0.036 | 0.0001 (0.0002) | 0.023–0.050 | −0.003, 0.37 |
| E2F | 0.027 | 0.031 (0.031) | 0.008–0.034 | 0.026, 0.09 |
| Wnt | 0.021 | 0.010 (0.014) | 0.007–0.027 | −0.005, 0.43 |
| p53 | 0.167 | 0.0001 (0.0002) | 0.084–0.210 | **0.075, 0.011** |

- dMMR is concentrated in batch 5 (30 of 59) and in CIT subtype C2 (50 of 59).
- **Only the p53 association persists within C2.** The EMT, E2F and Wnt associations with MMR largely mirror subtype composition.
- Power for MMR-stratified outcome models is very low: there are only 10 recurrences among dMMR tumours.
- The frozen Phase 8C values (ε² 0.0394 / 0.0251 / 0.0187 / 0.1716, n = 410) came from the 453-patient S3 cohort; the 452-patient values here differ slightly for that reason.

## X5. Composition estimates (computational, not measured)

**ESTIMATE** (Yoshihara 2013; v1.0.13):
- Re-implemented in Python from the package source; the package was **read, not installed or executed**.
- **Validation gate:** the bundled example `sample_input.gct` → `sample_estimate.gct` is reproduced to 1.4 × 10⁻¹⁴ relative error.
- Coverage: 10,035 of 10,412 common genes; stromal set 136/141; immune set 137/141.
- Observation: the bundled `sample_purity.gct` disagrees with the package's own purity formula by up to 4.5 %; it appears to predate that formula. Purity here follows the package code. As a monotone transform of the ESTIMATE score, it does not affect rank-based results.

**MCP-counter** (Becht 2016; v1.2.0):
- Native HG-U133 Plus 2.0 probe-set mode on the deposited probe matrix; all marker probe sets are present.
- Gene-symbol mode agrees, Pearson r ≥ 0.95 for every population except CD8 T cells (two probes, r = 0.12).

**Descriptives (n = 452):** ESTIMATE purity has median 0.74 (range 0.28–0.95), and none was out of bounds. The purity formula was calibrated on unbatch-corrected Affymetrix data, and here 10,035 of 10,412 common genes were used. **Purity is therefore treated as an uncalibrated, relative index.**

## X6. Scores versus composition; EMT marker overlap

**Overlap with stromal and immune marker sets:**

| Pathway | ESTIMATE stromal | ESTIMATE immune | MCP fibroblasts | MCP endothelial | Stromal union |
|---|---|---|---|---|---|
| EMT (198) | 26 (13.1 %) | 5 (2.5 %) | 6 (3.0 %) | 0 | **29 (14.6 %)** |
| E2F (195) | 0 | 0 | 0 | 0 | 0 |
| Wnt (41) | 0 | 0 | 0 | 0 | 0 |
| p53 (194) | 0 | 4 (2.1 %) | 0 | 0 | 0 |

**Spearman ρ (bootstrap 95 % CI, 2,000 resamples, n = 452):**
- **EMT** vs MCP fibroblasts **0.95** (0.94–0.96); vs ESTIMATE stromal **0.91**; vs endothelial 0.71; vs purity −0.80.
- **EMT-ex** (169 genes; the 29 stromal-marker genes removed): correlation with the original EMT score 0.996. Against stromal estimates: fibroblasts **0.93** (0.91–0.94), ESTIMATE stromal **0.89** (0.86–0.91).
- **E2F**: modest negative correlation with stromal content (−0.32) and positive with purity (+0.30).
- **p53**: modest positive correlation with immune and cytotoxic estimates (0.32 / 0.36).
- **Wnt**: |ρ| ≤ 0.18 with every estimate.
- The two independent estimators agree with each other (ESTIMATE stromal vs MCP fibroblasts ρ = 0.93).

**Registered rule met** (|ρ| ≥ 0.7 for EMT-ex). The rule needs two qualifications:
- All estimates and scores come from the same arrays, so shared technical variation (RNA quality, sampling cellularity, batch) inflates these correlations.
- The 169 genes remaining in EMT-ex are largely fibroblast-expressed and co-regulated, so "no shared genes" does not mean an independent measurement (EMT-ex vs EMT ρ 0.996).

With those caveats, the bulk EMT score co-varies almost entirely with estimated stromal and fibroblast content, even without shared genes. This is *consistent with* a predominantly stromal source of the bulk EMT signal, as reported by Calon 2015 and Isella 2015. It does **not** establish the cellular source, and it does not exclude a tumour-cell EMT component.

## X7. Composition-adjusted Cox models (exploratory)

Unpenalised Cox models; n = 452, 136 recurrences. HR per SD, Wald 95 % CI:

| Score | Frozen ridge M2 | R0 unpenalised | R1 + ESTIMATE stromal, immune | R2 + MCP fibroblasts, endothelial, cytotoxic |
|---|---|---|---|---|
| EMT | 1.070 (0.951–1.202) | 1.12 (0.93–1.35) | 1.81 (0.98–3.36) | 1.54 (0.86–2.76) |
| E2F | 0.906 (0.807–1.016) | 0.82 (0.69–0.98) | 0.78 (0.65–0.93) | 0.80 (0.66–0.96) |
| Wnt | 1.016 (0.906–1.138) | 1.04 (0.88–1.23) | 0.97 (0.81–1.16) | 0.99 (0.83–1.19) |
| p53 | 0.995 (0.882–1.123) | 0.99 (0.81–1.22) | 0.97 (0.76–1.24) | 1.01 (0.81–1.26) |
| StromalScore | — | — | 0.61 (0.29–1.29) | — |
| MCP fibroblasts | — | — | — | 0.80 (0.45–1.43) |

- Composition alone (C1 = clinical + stromal + immune): stromal HR 1.29 (1.01–1.65), immune HR 0.82 (0.64–1.07).
- Variance inflation factors in R1: EMT 12.8, StromalScore 19.0. In R2: EMT 11.3, fibroblasts 10.5.
- Likelihood-ratio tests for added blocks (unpenalised, apparent; deviation 8):

| Comparison | χ² | df | p |
|---|---|---|---|
| Four scores added to clinical model (R0 vs C0) | 9.19 | 4 | 0.056 |
| ESTIMATE stromal + immune added (R1 vs R0) | 5.02 | 2 | 0.081 |
| MCP block added (R2 vs R0) | 3.56 | 3 | 0.31 |
| Composition added to clinical model (C1 vs C0) | 3.99 | 2 | 0.14 |

- MMR-adjusted sensitivity (deviation 8; 409 patients with MMR status, 129 recurrences), HR (95 % CI):

| Score | R0 in this subset | R0 + dMMR |
|---|---|---|
| EMT | 1.14 (0.94–1.39) | 1.17 (0.96–1.42) |
| E2F | 0.85 (0.71–1.02) | 0.87 (0.72–1.05) |
| Wnt | 1.03 (0.87–1.23) | 0.98 (0.82–1.18) |
| p53 | 0.98 (0.79–1.20) | 1.03 (0.83–1.28) |
| dMMR | — | 0.49 (0.23–1.01) |

The ComBat caveat in X4 applies to this model.

**Exploratory discrimination** (single-level bootstrap, B = 200, event-stratified, every step refitted in each resample; no failures):

| Model | Apparent C | Optimism | Corrected C |
|---|---|---|---|
| C0 clinical | 0.6242 | 0.020 | 0.604 |
| C1 clinical + composition | 0.6271 | 0.026 | 0.601 |
| R0 clinical + scores | 0.6418 | 0.032 | 0.610 |
| R1 R0 + ESTIMATE | 0.6516 | 0.034 | 0.617 |
| R2 R0 + MCP | 0.6468 | 0.041 | 0.606 |

**Interpretation:**
- **EMT:** with composition in the model, EMT and the stromal estimate split one shared signal with opposite-signed, inflated coefficients (VIF > 10). On its own the stromal estimate has HR 1.29 (C1); next to EMT it flips to 0.61. That is one latent stromal factor divided between two near-collinear terms. **The EMT and stromal coefficients in R1/R2 must not be interpreted individually.** This is collinearity and possible over-adjustment, not evidence that adjustment revealed a hidden EMT effect.
- **E2F:** its estimate is essentially unchanged by composition adjustment. It shares no genes with the stromal sets and correlates only modestly with them.
- **Composition:** it adds no discrimination beyond clinical variables.
- **Discrimination differences:** all corrected-C differences are ≤ 0.013, carry no interval, and are of the same order as the Monte Carlo and sampling uncertainty. None is evidence of improvement, and none is comparable with the frozen nested-bootstrap ΔC estimand.
- **Over-adjustment and measurement error:** composition may lie on the path between tumour biology and recurrence (stromal reaction is itself biology). The estimates are noisy proxies computed from the same arrays as the scores.
- **Cross-check:** C0's apparent C reproduces the frozen M1 apparent C, 0.6242.

## X8. External-cohort eligibility audit

See `docs/EXTERNAL_COHORT_ELIGIBILITY_AUDIT.md`.
- GEO search: 554 hits, 85 screened, 14 assessed. An ArrayExpress search found no additional eligible cohort.
- **No cohort meets the four compatibility criteria**, so no frozen-model validation was run.
- The largest relevant cohort has 124 stage II–III events (GSE24551, exon arrays). Its approximate calibration-slope 95 % CI width is 0.66, above the Riley relaxed target of 0.4 (about 337 events).

## X9. Recurrence-time unit

See `reports/exploratory/x9_time_unit/x9_summary.json`. Six further sources were searched:
- the GEO series and SuperSeries records;
- the ArrayExpress study record;
- the CIT programme website, which no longer resolves;
- two literature searches.

The publication documents the time origin (surgery) and the event (first recurrence). **No primary source states the unit of the deposited column.** Status is unchanged: `[PROVISIONAL: months]`.

## Deviations from the registered plan

1. **X4 sample size.** The registered population is the 452-patient A-19 cohort, which gives n = 409 with MMR. The frozen Phase 8C values (n = 410) used the 453-patient S3 cohort. Both are reported.
2. **X5 purity.** The plan's gate (scores vs `sample_estimate.gct`) passed. An additional purity comparison against `sample_purity.gct` showed that the bundled file predates the package's own formula. This is recorded as an observation, not a gate.
3. **X6 EMT-ex.** It is computed with the frozen scoring formula but without the frozen minimum-coverage rule, which would otherwise drop a reduced set. The coverage remained 169 of 198 genes (85 %).
4. **X8.** ArrayExpress was searched through the BioStudies API. A query for "PETACC" returned an error page and is not counted; PETACC-3 expression data were therefore not assessed from the repository.
5. **X2 estimator (after independent review).** lifelines' `AalenJohansenFitter` breaks ties between causes by random jitter. With 109 tied event times among 134 distinct times, the first-run cumulative incidences depended on that jitter: stage III death without recurrence at 60 was 4.4 % vs 5.2 % deterministic. It was replaced by R `survival::survfit` multistate (tie-aware, deterministic). The tables and figure XF2 were regenerated.
6. **X3 profile grid.** Refined to 0.0025 steps around both interval bounds. The bounds are 0.0275 and 0.50 (grid points inside the 95 % cut-off).
7. **X1 multiplicity.** BH is also reported across the 20 non-CIT tests, because the CIT tests are circular.
8. **X7 additions (after independent review).** Likelihood-ratio tests for added blocks, and an MMR-adjusted sensitivity model. Neither was in the registered plan; both are labelled post-review sensitivity.

## Independent review

A separate reviewer wrote independent implementations, sharing only the frozen loaders. The scripts and outputs are kept outside the repository in the session scratch area. The reviewer recomputed:
- X1 (TP53, BRAF tests);
- X2 (event counts, Aalen–Johansen at 36 and 60);
- X5 (ESTIMATE for all 452 to ≤ 5 × 10⁻¹⁰; MCP fibroblasts to 1.8 × 10⁻¹⁵);
- X6 (EMT–fibroblast ρ 0.94981; 26 shared genes);
- X7 (R0 and R1 HRs, VIFs);
- X8 (GSE14333 185/48; GSE103479 156/42).

**Everything matched except the original X2 cumulative incidences**, which were corrected (deviation 5) and now agree with the reviewer's estimate (0.3069 at 60). The reviewer found no leakage in the X7 bootstrap. Its interpretive critiques are addressed in X1, X2, X4, X5, X6 and X7 above, and in `docs/SCIENTIFIC_INTERPRETATION_AND_CONCLUSION.md`.
