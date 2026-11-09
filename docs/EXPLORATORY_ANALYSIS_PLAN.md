# Registered exploratory analysis plan (post-primary)

**Status:** registered on 8 November 2026, before any of the analyses below was run. Its SHA-256 is recorded in `reports/exploratory/plan_registration.json` at registration time.
**Owner approval:** Tiers 1–3 (analyses X1–X9) approved by the project owner on 8 November 2026. Single-cell reference integration was explicitly excluded.
**Relationship to the primary analysis:** every analysis here is **exploratory or sensitivity work, registered after the primary results were known**. None replaces, refits or reinterprets the frozen primary analysis (M1, M2 four-set primary, five-set sensitivity, S-7, nested-bootstrap ΔC estimand and intervals, calibration). The primary conclusion is fixed: corrected ΔC +0.004316, all five 95 % intervals include zero.

## 0. Common definitions

| Item | Definition |
|---|---|
| Population | The frozen A-19 analysis cohort: 452 patients, 136 recurrences, loaded read-only with `src/phase9_corrected/data.py::load_context` |
| Expression | `data/processed/analysis/expression_median_S0.parquet` (deposited RMA + ComBat, median of probes per gene), unchanged. MCP-counter (X5) also uses the deposited probe matrix `data/processed/probe_matrix.parquet` in its native HG-U133 Plus 2.0 probe-set mode |
| Pathway scores | The frozen four-set definitions (`config/pathways/m2_panel_4set.json`), computed with the frozen `src/pathway_scoring.score_panel`. Because these analyses are cohort-level descriptions, scores are standardised on all 452 patients (not inside folds). Checked against `reports/phase8c/scores_4set_primary_reproduced_full_cohort.csv` |
| Endpoint | `rfs_time_full` / `rfs_event_full` (time unit **provisional: months**, undocumented) |
| Clinical design | Frozen M1 design (`src/phase9_execution.build_clinical`): age restricted cubic spline (3 knots, 10/50/90th percentiles), sex, stage III vs II, proximal vs distal, adjuvant chemotherapy |
| Code and outputs | `scripts/exploratory/`, `reports/exploratory/<analysis>/`, downloads under `data/external/exploratory/` with a provenance record. No frozen file is written |
| Software | Project `.venv` (Python, lifelines 0.30.3, scipy); R 4.5.2 with the base `survival` package only |
| Seeds | Every resampling uses a fixed, recorded seed |

## X1. Molecular anchoring (outcome-free)

- **Question:** do the expression scores track the genomic or molecular features expected from their biology?
- **Annotations (GEO deposit):** `tp53.mutation`, `kras.mutation`, `braf.mutation` (M vs WT), `cimp.status`, `cin.status` (+ vs −), `cit.molecularsubtype` (6 levels). Missing or "N/A" values are excluded and counted.
- **Tests:** Mann–Whitney U for two-group annotations. Effect sizes are the rank-biserial correlation and the difference in medians (in score SD units). Kruskal–Wallis with ε² for CIT subtype. **4 scores × 6 annotations = 24 primary tests, Benjamini–Hochberg FDR across all 24.**
- **Directional expectation (stated in advance):**
  - p53 score lower in TP53-mutant tumours (loss of p53 target transactivation).
  - Other comparisons are two-sided and without a stated direction.
- **Pre-specified secondary analysis:** because MMR status is associated with both TP53 status and the p53 score, the TP53 comparison is repeated within pMMR tumours, with a batch-stratified permutation test (10,000 permutations).
- **Interpretation rule:** an association links expression to mutation or subtype. It is **not** evidence of pathway function or activation.

## X2. Competing risks

- **Event coding:**
  - recurrence = `rfs_event_full == 1`;
  - death without recurrence = `rfs_event_full == 0` and `os_event_source == 1` with `os_time_source == rfs_time_full`;
  - everyone else is censored.
  - The coding is verified to reproduce the 46 deaths without recurrence in Phase 6A before any estimate is made. Any mismatch stops the analysis.
- **Estimands:**
  - cumulative incidence of recurrence and of death without recurrence (Aalen–Johansen), overall and by stage, at 36 and 60 time units, with 95 % confidence intervals;
  - subdistribution hazard ratios (Fine–Gray; R `survival::finegray` + weighted `coxph`) for the M1 design and for an unpenalised M1 + four-score design;
  - reference comparison: the cause-specific Cox estimates (unpenalised) for the same designs, with both event types.
- **Labels:** all absolute quantities carry "time unit provisional: months". Subdistribution HRs are compared with cause-specific HRs only qualitatively. The primary conclusion is not changed by this analysis.

## X3. Batch estimability and frailty sensitivity

- **Descriptive:** per-batch n, recurrences, deaths without recurrence, stage III %, chemotherapy %, partition composition and score means.
- **Instability diagnostics** from the frozen leave-one-batch-out (LOBO) results (read, not refitted):
  - coefficient change per pathway as a multiple of its full-cohort standard error;
  - Spearman correlation of the change with the omitted batch's events.
- **Frailty sensitivity models (R `survival`):** unpenalised clinical + four-score design plus `frailty(batch, distribution = "gamma")`, and a Gaussian frailty as a check. Reported:
  - the frailty variance θ, with a profile-likelihood 95 % interval from a θ grid;
  - the likelihood-ratio comparison with the no-frailty model;
  - coefficient changes.
- **Classification rules (fixed now):**
  - *Resolved to a defensible degree* requires all three of: θ's upper 95 % bound below 0.05; no pathway coefficient whose LOBO range crosses zero beyond sampling variation (all LOBO changes within 1.96 full-cohort SE); and no residual batch–score association after stratification.
  - *Partially mitigated* requires a frailty model that is estimable, with a bounded θ, and coefficients stable under it, while LOBO sign changes remain.
  - Otherwise *characterised, not resolved*.
  - *Not estimable* applies if θ cannot be bounded with 9 batches.

## X4. MMR status and `combat_mod` provenance

- **Provenance:** trace `mmr.status` in GEO SOFT, the ArrayExpress SDRF and Marisa et al. 2013 Table S1, and `combat_mod` in each source. Record the field names and any methods text describing the ComBat model matrix.
- **Robustness of score–MMR associations:**
  - ε² overall;
  - a batch-stratified permutation test (10,000 permutations);
  - ε² recomputed leaving out each batch in turn;
  - ε² within CIT subtypes where both MMR groups have ≥ 10 patients.
- **Power:** dMMR / pMMR counts overall and per batch.
- **Rule:** `combat_mod` is never relabelled as MMR status.

## X5. Composition estimation (computational estimates, not measurements)

- **ESTIMATE** (Yoshihara et al., *Nat Commun* 2013): re-implemented in Python from the package source (`estimateScore.R`, v1.0.13). The input is restricted to the package's common-gene list, as `filterCommonGenes` does. Outputs are StromalScore, ImmuneScore, ESTIMATEScore and Affymetrix purity, where purity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore).
  - **Validation gate:** the re-implementation must reproduce the package's bundled example output (`sample_input.gct` → `sample_estimate.gct`) to within 1e-6 relative error, or X5–X7 stop.
- **MCP-counter** (Becht et al., *Genome Biol* 2016): the arithmetic mean of log2 expression of each population's HG-U133 Plus 2.0 marker probe sets (the native mode), on the deposited probe matrix. Gene-symbol mode on the median gene matrix is computed as a check.
- **Reported:** coverage of each marker set and per-sample estimates. These are arbitrary-unit abundance scores (MCP-counter) and enrichment scores (ESTIMATE), **not fractions and not measurements**.

## X6. Scores versus estimated composition; EMT marker overlap

- **Spearman ρ with bootstrap 95 % intervals** (2,000 resamples, seed 20261108) between each of the four scores and StromalScore, ImmuneScore, purity, MCP fibroblasts, endothelial cells, cytotoxic lymphocytes, T cells and monocytic lineage. These are descriptive; no hypothesis tests.
- **Overlap:** the number and proportion of Hallmark EMT genes (198) in the ESTIMATE stromal set, the ESTIMATE immune set, and the MCP fibroblast and endothelial markers. The same is done for the other three pathways.
- **EMT-ex (exploratory sensitivity score):** the EMT score recomputed after removing every gene in the ESTIMATE stromal set ∪ MCP fibroblast markers ∪ MCP endothelial markers. The original EMT score is unchanged.
  - Reported: correlation of EMT-ex with the original EMT score and with the stromal estimates.
- **Interpretation rule:** if EMT-ex remains strongly correlated with the stromal estimates (|ρ| ≥ 0.7) despite sharing no genes with them, the bulk EMT signal co-varies with stromal content. This is *consistent with* a substantial stromal contribution. It does **not** establish the cellular source in either direction.

## X7. Composition-adjusted Cox models (exploratory outcome analysis)

- **Models**, all unpenalised Cox on the A-19 cohort with the frozen clinical design:
  - C0 = clinical;
  - C1 = clinical + ESTIMATE StromalScore + ImmuneScore;
  - R0 = clinical + four scores (an unpenalised analogue of M2, used as the reference for coefficient changes);
  - R1 = R0 + StromalScore + ImmuneScore;
  - R2 = R0 + MCP fibroblasts + endothelial cells + cytotoxic lymphocytes.
- **Estimates:**
  - HR (Wald 95 % CI) for each pathway score in R0, R1 and R2, and the change in log-HR;
  - HRs for the composition terms;
  - the frozen ridge M2 HRs are shown alongside for reference.
- **Exploratory discrimination:** single-level Harrell bootstrap optimism correction, B = 200, seed 20261108. Knots, standardisation and all fits are redone within each resample. Reported: apparent C, optimism and corrected C for C0, C1, R0, R1 and R2.
  - **Not comparable** with the frozen nested-bootstrap ΔC estimand. No interval is claimed and no model is promoted.
- **Interpretation rules:**
  - Composition adjustment can over-adjust (composition may lie on the biological path between pathway programmes and recurrence) and carries measurement error (estimates are proxies).
  - A changed coefficient is therefore read as "the score's association is not separable from estimated composition", not as confounding removed.

## X8. External-cohort eligibility audit

- **Search:** NCBI GEO DataSets (E-utilities) with recorded queries and date, plus ArrayExpress/BioStudies, plus cohorts named in the literature already reviewed (Phase 4).
- **Screening:** human primary colorectal tumours; per-sample gene expression; per-sample recurrence or disease-free outcome with time; stage available.
- **Eligibility assessment for each screened cohort:**
  - independence from GSE39582/CIT (shared institutions, overlapping deposits);
  - endpoint definition and death handling;
  - stage II–III subset and events;
  - follow-up;
  - availability of all five M1 predictors;
  - platform (exact reproduction of the frozen scores requires HG-U133 Plus 2.0 or a platform covering the 607 genes; a different platform changes the measurement);
  - expected precision. The approximate standard error of C and of the calibration slope are derived from the event count, following the Phase 5 precision assessment (Riley et al. 2021 / Collins et al. 2016). This replaces a fixed event threshold.
- **Frozen-model validation** is run only if a cohort is independent, has a compatible recurrence endpoint, provides all M1 predictors and supports exact score reproduction. Otherwise no validation is run, and the audit documents why. An exploratory direction check is reported only if labelled as such and never pooled.

## X9. Recurrence-time unit

A further documented search beyond the 8 sources in `PHASE6A_ENDPOINT_RESOLUTION.md`:
- the GEO series record text;
- the ArrayExpress experiment description;
- later CIT-programme publications that reuse GSE39582 (searched by DOI and title);
- the CIT programme website.

The unit is upgraded only by an explicit statement from a primary source.

## Multiplicity and reporting

- X1 uses BH-FDR across its 24 primary tests.
- X4 tests are reported with their unadjusted and BH-adjusted values within X4.
- X6 and X7 are estimation-only.
- Every result is reported with its sample size, event count where relevant, and label (exploratory / sensitivity).

## Deviations

Any deviation from this plan is recorded in `docs/ADDITIONAL_ANALYSES_METHODS_AND_RESULTS.md` §Deviations, with the reason.
