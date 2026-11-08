# Phase 6D — Definitive Missingness Audit

**Date:** 2026-11-02 · **Status:** Audit complete. **No imputation, no modelling, no feature selection, no normalization.**
**Reproduce:** `.venv/bin/python scripts/28_definitive_missingness_audit.py` → `reports/phase6d/`

**Outcome handling:** recurrence and survival **values were never read**. Endpoint completeness is **cited** from the verified Phase 6A result rather than recomputed, so no outcome column was opened by this audit. No outcome informed any variable selection or missingness rule.

---

## 0. Headline

> **Multiple imputation is NOT indicated for any planned model.**
>
> M1, M2 (4-set), M2 sensitivity (5-set) and M3 all have **452 / 453 complete cases (99.8%)** — complete-case analysis loses **one patient**. Pathway inputs have **zero** missing values.
>
> The substantial missingness (tp53 36.9%, cin 19.6%, cimp 14.8%) is confined to **optional molecular covariates that no planned model uses**, and its mechanism is **assay availability determined by processing batch** — not a missing-data problem to be imputed around.

---

## 1. Missing-value encoding — the codes are NOT interchangeable, and the deposit conflates them

The brief's caution was warranted. The GSE39582 deposit **documents** a distinction in its own `!Sample_description`:

> `ND: not determined`
> `NA: Not applicable`
> *"gene.mutation=WT (wild type), the exon, protein and DNA columns are NA (not applicable)."*

**But `ND` is used zero times across all 37 metadata columns.**

| Finding | |
|---|---|
| Occurrences of `ND` anywhere | **0** |
| Code actually used | **`N/A`**, for *both* meanings |
| Empty strings | 0 |
| Literal `NA` | 0 |
| Consequence | **"not determined" and "not applicable" are encoded identically and cannot be separated from the encoding alone** |

**My preliminary (Phase 6C-era) audit conflated them**, because `pandas` parses `N/A` to `NaN` by default. This audit reads with `keep_default_na=False` and separates the two meanings **structurally**, by verified logical rule.

### 1.1 Structural separation — verified

| Rule | Verification |
|---|---|
| `mutation = WT` → detail columns `N/A` (**not applicable**) | tp53 **161/161**, kras **328/328**, braf **461/461** — holds exactly |
| `chemotherapy.adjuvant = N` → `type = N/A` (**not applicable**) | **316 / 316** |
| `chemotherapy.adjuvant = Y` but `type = N/A` (**truly missing**) | **66** |

So `chemotherapy.adjuvant.type`'s 382 `N/A` values split into **316 structural** and **66 genuinely missing** — a distinction invisible to naive `isna()`. *(That variable is excluded from all models by Phase 1 regardless.)*

---

## 2. Variable-level missingness (S3, n = 453)

| Variable | group | observed | missing | % | disc % | val % | **diff pp** |
|---|---|---:|---:|---:|---:|---:|---:|
| **age at diagnosis** | **M1 primary** | 453 | **0** | **0.00** | 0.00 | 0.00 | 0.00 |
| **sex** | **M1 primary** | 453 | **0** | **0.00** | 0.00 | 0.00 | 0.00 |
| **tnm.stage** | **M1 primary** | 453 | **0** | **0.00** | 0.00 | 0.00 | 0.00 |
| **tumor.location** | **M1 primary** | 453 | **0** | **0.00** | 0.00 | 0.00 | 0.00 |
| **chemotherapy.adjuvant** | **M1 primary** | 452 | **1** | **0.22** | 0.28 | 0.00 | −0.28 |
| tnm.t / tnm.n / tnm.m | M1 optional | 434 | 19 | 4.19 | 3.11 | 8.08 | +4.97 |
| mmr.status | M1 optional | 410 | 43 | 9.49 | 9.04 | 11.11 | +2.07 |
| kras.mutation | molecular optional | 432 | 21 | 4.64 | 5.37 | 2.02 | −3.35 |
| braf.mutation | molecular optional | 399 | 54 | 11.92 | 14.41 | 3.03 | −11.38 |
| cimp.status | molecular optional | 386 | 67 | 14.79 | 17.23 | 6.06 | −11.17 |
| cin.status | molecular optional | 364 | 89 | 19.65 | 21.47 | 13.13 | −8.34 |
| **tp53.mutation** | molecular optional | 286 | **167** | **36.87** | **43.50** | **13.13** | **−30.37** |

**Every M1 primary predictor is complete except one chemotherapy value.**

**Endpoint and provenance variables:** endpoint missingness in S3 is **zero by construction** (Phase 1 eligibility requires documented RFS; verified in Phase 6A §4.2 as an exhaustive 136 + 46 + 271 = 453 partition with zero missingness in all four endpoint fields). Provenance fields — sample ID, cohort membership, dataset, batch — are complete by construction. **Neither category was re-opened here.**

**Pathway inputs (M2/M3):** the S3 expression matrix has **0 NaN values and 0 zero-variance genes**. Pathway inputs are complete.

---

## 3. Missingness patterns

20 distinct patterns over 14 audited variables; **203 samples (44.8%) complete on all 14**.

**Co-missingness — assays travel together:**

| pair | both missing | Jaccard |
|---|---:|---:|
| **tnm.t ↔ tnm.n ↔ tnm.m** | 19 | **1.000** |
| cimp.status ↔ braf.mutation | 52 | 0.754 |
| cimp.status ↔ cin.status | 60 | 0.625 |
| cin.status ↔ braf.mutation | 50 | 0.538 |

TNM sub-staging is missing as a **block** (Jaccard exactly 1.000 — the same 19 records), consistent with a single absent pathology record rather than three independent omissions. The molecular assays cluster similarly.

---

## 4. Mechanism — evidence, hypothesis, and what remains unresolved

### 4.1 Observed evidence

`tp53.mutation` missingness **stratified by ComBat processing batch**:

| batch | n | **tp53 missing %** | cin % | cimp % |
|---|---:|---:|---:|---:|
| batch01 | 33 | **0.0** | 15.2 | 0.0 |
| **batch02** | 14 | **100.0** | 0.0 | 0.0 |
| **batch03** | 102 | **100.0** | 19.6 | 16.7 |
| batch04 | 103 | **0.0** | 55.3 | 43.7 |
| **batch05** | 49 | **100.0** | 2.0 | 8.2 |
| batch07 | 81 | 1.2 | 2.5 | 0.0 |
| batch08 | 41 | 2.4 | 0.0 | 0.0 |
| batch09 | 15 | **0.0** | 26.7 | 0.0 |
| batch10 | 15 | **0.0** | 0.0 | 6.7 |

**`tp53.mutation` missingness is near-deterministic in batch: three batches are 100% missing, four are 0%.** `cin.status` ranges 0–55.3% and `cimp.status` 0–43.7% across batches. By contrast, stratifying by **stage** moves tp53 missingness only 29.0%→43.1%, and by **sex** or **location** even less.

### 4.2 Assessment against the three mechanisms

| Mechanism | Verdict | Basis |
|---|---|---|
| **MCAR** | **Inconsistent with the evidence** | Missingness varies from 0% to 100% across batches and by 30.4 pp between subsets. Not plausibly random |
| **MAR given observed variables** | **Consistent, conditional on batch being included** | Batch is observed and recorded in the Phase 3 bundle. Conditional on batch, missingness is close to deterministic — which is the strongest form of "explained by observed data" |
| **MNAR / structural** | **Partly — and this is the better description** | A batch with 100% missing tp53 was not *sampled* from a distribution; **the assay was not run for that batch**. This is structural unavailability, not a stochastic missing-data process |

> **Stated as the brief requires: neither MAR nor MNAR is declared proven.** These are assessments of consistency with observed patterns, not statistical proofs. A near-deterministic relationship with batch is strong evidence that the mechanism is *assay availability*, but the audit cannot exclude that the decision to assay a batch was itself influenced by unobserved factors.

### 4.3 Hypothesis — stated as hypothesis

**Hypothesis:** the molecular assays (TP53 sequencing, CIMP, CIN/CGH) were performed on a **subset of the CIT collection selected by processing batch**, most plausibly because batches were accrued at different centres or times and not all underwent the full molecular workup. The GSE39582 overall design notes that of 598 collected tumours, **463 were also analysed by CGH array** — a documented partial-assay design consistent with this pattern.

**Unresolved:** the deposit does not state which samples received which assay, or why. Resolving this would require the CIT programme's own documentation. **It is not needed for any planned model** (§5).

### 4.4 A consequence worth recording

Batch is **confounded with the discovery/validation partition** (Phase 4: Cramér's V = 0.488, two discovery-only batches). Since molecular missingness is near-determined by batch, **molecular missingness is correlated with partition as a downstream consequence** — the −30.4 pp tp53 differential is a symptom of the batch structure, not an independent phenomenon.

---

## 5. Model-specific assessment

| Model | predictors | **complete cases** | % | lost | disc | val |
|---|---:|---:|---:|---:|---:|---:|
| **M1 primary (clinical only)** | 5 | **452** | **99.8** | **1** | 353 | 99 |
| **M2 primary (4-set) + M1** | 9 | **452** | **99.8** | **1** | 353 | 99 |
| **M2 sensitivity (5-set) + M1** | 10 | **452** | **99.8** | **1** | 353 | 99 |
| **M3 exploratory (21,815 genes) + M1** | 21,820 | **452** | **99.8** | **1** | 353 | 99 |
| M1 + optional TNM/MMR | 9 | 390 | 86.1 | 63 | 310 | 80 |
| M1 + all molecular *(NOT PROPOSED)* | 14 | 203 | 44.8 | 250 | 139 | 64 |

**All four planned models lose exactly one patient to complete-case analysis**, because pathway and transcriptome inputs are complete and the only missing M1 value is a single chemotherapy record.

**Optional molecular covariates are NOT added to M1** (per the brief and Phase 6B §3.1). The 44.8% figure is reported only to show what including them *would* cost.

**M2 frozen specification is untouched:** memberships, fold-specific scoring, and the ≥80%/≥20-gene coverage rule are unchanged. **No expression value or pathway member is imputed.** The coverage rule is evaluated per fold at scoring time, not here.

---

## 6. Decision table

| Variable / input | group | missing | **Imputation indicated?** | Basis |
|---|---|---:|---|---|
| age at diagnosis | M1 primary | 0 | **UNNECESSARY** | complete |
| sex | M1 primary | 0 | **UNNECESSARY** | complete |
| tnm.stage | M1 primary | 0 | **UNNECESSARY** | complete |
| tumor.location | M1 primary | 0 | **UNNECESSARY** | complete |
| **chemotherapy.adjuvant** | M1 primary | **1** | **UNNECESSARY** | 1 of 453. Complete-case loses 0.22%; imputing one value from four observed covariates adds no information |
| tnm.t / n / m | M1 optional | 19 | **UNRESOLVED** | only if these enter M1 — a decision not yet made. Missing as a block (Jaccard 1.000) |
| mmr.status | M1 optional | 43 | **UNRESOLVED** | same condition |
| kras / braf / cimp / cin / **tp53** | molecular optional | 21–167 | **NOT APPLICABLE** | **no planned model uses them.** If ever added, MI would be indicated *and* would have to condition on batch |
| **pathway score inputs (M2)** | expression | **0** | **UNNECESSARY** | 0 NaN, 0 zero-variance genes. Frozen exclusion + coverage rule governs, **never substitution** |
| transcriptome (M3) | expression | 0 | **UNNECESSARY** | as above |
| endpoint fields | endpoint | **0** | **PROHIBITED** | zero missing by construction; imputation of endpoints is forbidden (Phase 6B §6.4.3) |
| provenance fields | provenance | 0 | **PROHIBITED** | forbidden regardless |

| Model | **Verdict** |
|---|---|
| **M1 primary** | **Imputation UNNECESSARY** — complete-case loses 1 patient |
| **M2 primary (4-set)** | **Imputation UNNECESSARY** |
| **M2 sensitivity (5-set)** | **Imputation UNNECESSARY** |
| **M3 exploratory** | **Imputation UNNECESSARY** |
| M1 + optional TNM/MMR | **UNRESOLVED** — depends on whether those covariates are adopted (86.1% complete case) |

---

## 7. Limitations

1. **Pre-expression-QC.** Array-level QC has not run. If arrays are later excluded, the denominator changes — though all four planned models are at 99.8%, so the conclusion is robust to small changes.
2. **Mechanism is assessed, not proven.** §4.2 states consistency with observed patterns; no statistical test can establish MAR or MNAR.
3. **The ND/NA conflation is irreducible.** The deposit never used `ND`, so "not determined" and "not applicable" are separable only where a documented logical rule exists (§1.1). For `tp53.mutation` itself no such rule exists — its `N/A` is presumed "not determined", but that is an inference.
4. **Why the assay subset was chosen is unknown**, and cannot be resolved from the deposit.
5. **Single-value imputation of the one chemotherapy record was not attempted**, since the recommendation is complete-case for that model.
6. **No outcome column was opened**, so no analysis of whether missingness relates to outcome was performed — deliberately, and it would require explicit approval.

---

## 8. A-19 status

**A-19 remains PROPOSED.** The audit supplies the evidence; the scientific decisions are the decision-maker's.

| Decision | Evidence from this audit | Status |
|---|---|---|
| Is MI needed for M1 / M2 / M2-sens / M3? | **No** — 452/453 complete for all four | **awaiting approval** |
| Primary strategy for planned models | **complete-case** (loses 1 patient) | **awaiting approval** |
| Do optional TNM/MMR covariates enter M1? | would reduce to 390/453 (86.1%) | **unresolved — decision-maker's** |
| If they do, is MI then indicated? | likely yes, conditioning on **batch** | **unresolved** |
| Per-variable imputation models | not required under the recommendation | **moot unless the above changes** |
| *m*, iterations, seeds | moot unless MI adopted | **moot** |
| **Pooling: Rubin vs averaged predictions** | moot unless MI adopted | **still unresolved in principle** |
| Exclusion rule for high-missingness variables | not needed — they are already excluded by model design | **moot** |

**Recommended, pending approval:** adopt **complete-case analysis** for M1, M2, M2-sensitivity and M3, report the single excluded patient explicitly, and retain the MI machinery in the protocol as a **contingency** that activates only if the optional molecular covariates are later adopted — in which case it must condition on batch.
