# Phase 7C — QC Execution Report (rev. 2, DECISIONS APPROVED)

**Date:** 2026-11-02 · **Revision 2** — three scientific decisions approved and recorded.
**Status:** Executed. ✅ **All 453 S3 samples retained. 0 samples excluded.**
**Reproduce:** `.venv/bin/python scripts/31_execute_qc.py` → `reports/phase7c/`
**Specification:** Phase 7B rev. 2, thresholds read from `config/analysis/QC_SPECIFICATION.json` and **not redefined in code**.

**Restrictions honoured:** no exclusions · no outcome column opened · no event counts computed · no normalization · no ComBat · no fRMA · no packages installed · no CEL files · no modification of frozen Phase 2 / Phase 3 / M2 artifacts.

---

## 0. APPROVED DECISIONS (2026-11-02)

Record: `config/analysis/QC_DECISIONS_RECORD.json` · Spec now at **revision 3**.

### D1 — H-C1 amended ✅

| | |
|---|---|
| **Original metric** | *"sample-ID set equality between the expression matrix and the frozen Phase 3 bundle"* |
| **Defect** | Required equality between artifacts of **different scope** — the bundle is S0-scoped (461) **by design**, the matrix holds 585. Equality was never achievable; the halt fired on the criterion's own wording, not on a data problem |
| **Original execution result** | **PRESERVED** — `as_written_set_equality: false`, 585 vs 461, `status_as_written: "FAIL (halt triggered)"` |
| **Corrected metric** | the frozen bundle's sample-ID set is a **SUBSET** of the matrix, with **zero bundle IDs missing** |
| **Corrected rule** | every bundle sample must be present under an exactly matching identifier; any missing ID **halts QC**. Matrix-only samples are expected and are **not** a failure |
| **Rationale** | The integrity property that matters is that no bundle sample is missing or renamed. The 124 matrix-only samples reconcile exactly as non-S0 by construction: 19 non-tumoural, 97 stage 0/I/IV or missing-stage, and **8** stage II–III lacking documented RFS (469 − 461 = 8, matching Phase 1) |
| **Broadens the hard list?** | **No** — still 9 criteria |
| **Changes any threshold?** | **No** |

### D2 — All 453 S3 samples retained ✅

| | |
|---|---:|
| Samples retained | **453** |
| Samples excluded | **0** |
| Soft-flagged | **29** |
| **Investigation records preserved** | **29** |
| Any soft-flagged sample authorised for exclusion | **NO** |

> **Explicit statement, as required: no independent evidence demonstrated technical invalidity for any soft-flagged sample.**
>
> Of the 29 flagged, **26 had no independent technical evidence** of invalidity. The remaining 3 carried only **further statistical flags** — an extreme PCA score, or a degradation contrast deviating from the unflagged median — not a demonstrated technical fault.
>
> The sole direct technical measure available on this matrix, the AFFX 3′/5′ degradation contrast, **did not discriminate** flagged from unflagged arrays (ACTB 0.792 vs 0.772; GAPDH 0.259 vs 0.188).

### D3 — Two findings recorded as OBSERVATIONS, not exclusion justifications ✅

| ID | Observation | Status |
|---|---|---|
| **OBS-1** | **batch05 flag-rate association** — 20.4% (10/49) vs 0.0% in batch02/09/10; χ² = 23.2, p = 0.003, **V = 0.226**. Batch-*associated*, not batch-*determined*; markedly weaker than the near-deterministic patterning Phase 6D found for assay missingness | **OBSERVATION ONLY** — `is_exclusion_justification: false`. No sample excluded, no batch correction applied |
| **OBS-2** | **GSM972311 clustering singleton** — 452/1 split; Σ(1−ρ) = 57.72, the cohort maximum and ~1.8× the fence. Yet median 4.8096 and IQR 2.8450 both sit **inside** their fences, PCA z-scores are unremarkable (PC1 −0.48, PC2 −2.19), and the degradation contrast shows **no deviation** | **OBSERVATION ONLY** — `is_exclusion_justification: false`. Retained |

**Preserved:** all QC thresholds unchanged · specification revision history (rev 1 → 2 → 3) · all 29 investigation records · the original H-C1 execution result.

---

## 0. Summary

| | |
|---|---:|
| Criteria run | **20 / 20** |
| **HARD failures** | **0** (9 passed) |
| **SOFT criteria flagging** | **3 / 3** |
| **Samples soft-flagged** | **29 / 453 (6.4%)** |
| Descriptive findings reported | **8** |
| **SAMPLES EXCLUDED** | **0** |
| Outcome columns opened / event counts computed | **No / No** |
| **Issues requiring approval** | **0** — all resolved (§0) |

---

## 1. Cohort integrity — 4 criteria, 0 failures

| ID | Result |
|---|---|
| **H-C1** | ⚠️ **PASS with a SPECIFICATION DEFECT** — see §1.1 |
| **H-C2** | ✅ partition membership complete: **discovery 354, validation 99** |
| **H-C3** | ✅ S3 membership == frozen **n = 453** |
| **H-C4** | ✅ complete-case exclusion is **exactly GSM972262** |

### 1.1 ⚠️ H-C1 is mis-specified — THE issue requiring approval

**H-C1 as written demands set equality between the deposited matrix and the frozen Phase 3 bundle.** On execution it failed and halted.

**The failure is in the criterion, not the data.** The Phase 3 bundle is **S0-scoped (461) by design**; the deposited matrix holds all **585** arrays. **Equality was never achievable.**

The integrity property that actually matters was then evaluated:

| Check | Result |
|---|---|
| frozen bundle ⊆ deposited matrix | **TRUE** |
| **bundle samples missing from the matrix** | **0** |
| matrix-only samples | 124 |

The 124 extras reconcile exactly:

| Category | n |
|---|---:|
| non-tumoural mucosa | 19 |
| stage 0 / I / IV / missing-stage primary tumours | 97 |
| stage II–III primary tumours **lacking documented RFS** | **8** |

469 stage II–III primary tumours − 461 in S0 = **8**, and those 8 are precisely GSM971967, GSM972350, GSM972356, GSM972446, GSM972448, GSM972463, GSM972500, GSM972514 — the same records documented in Phase 1. **Every extra is legitimately outside S0.**

> **I did not silently amend the frozen criterion.** Both evaluations are recorded in `qc_summary.json` under `H-C1`, flagged `requires_approval: true`. Execution continued only because **no frozen-bundle sample is missing** — the data-integrity property is satisfied.
>
> **Requested decision:** approve amending H-C1 from *set equality* to *the frozen bundle is a subset of the matrix with zero missing IDs*. This narrows nothing and broadens nothing; it corrects a comparison between artifacts of different scope.

---

## 2. Sample-level hard criteria — 4 criteria, 0 failures

| ID | Metric | Result |
|---|---|---|
| **H-S1** | non-finite values per sample | ✅ **0 samples** |
| **H-S2** | samples absent from the frozen manifest | ✅ **0** |
| **H-S3** | exactly-duplicated sample columns (SHA-256 of column bytes) | ✅ **0 groups** |
| **H-S4** | non-tumoural samples in the analysis cohort | ✅ **0** |

---

## 3. Soft criteria — all three flagged; **no sample excluded**

Tukey 1.5 × IQR fences, computed on each metric's own distribution across the 453 S3 arrays, exactly as frozen.

| ID | Metric | Fence | Flagged |
|---|---|---|---:|
| **S-S1** | per-array median | [4.6907, 4.8532] | **9** |
| **S-S2** | per-array IQR | [2.8033, 3.4216] | **13** |
| **S-S3** | Σ(1 − Spearman ρ) to all other arrays | [14.701, 32.098] | **15** |

**29 distinct samples** carry at least one flag (6.4% of the cohort). Per the frozen policy, **every one produced a seven-field investigation record and none was excluded.**

### 3.1 Sample IDs requiring investigation (29)

| Sample | Flags | Batch | Subtype | Independent technical evidence |
|---|---|---|---|---|
| **GSM972059** | S-S1+S-S2+S-S3 | batch03 | C2 | **PARTIAL — extreme PCA \|z\| > 4** |
| **GSM972077** | S-S3 | batch03 | C6 | **PARTIAL — extreme PCA \|z\| > 4** |
| **GSM972105** | S-S1+S-S2+S-S3 | batch03 | C5 | **PARTIAL — degradation ratio > 1 log2 from unflagged median** |
| GSM971986 | S-S1+S-S3 | batch01 | C5 | none found |
| GSM971993 | S-S1 | batch01 | C3 | none found |
| GSM972061 | S-S2 | batch03 | C1 | none found |
| GSM972079 | S-S3 | batch03 | C2 | none found |
| GSM972089 | S-S3 | batch03 | C4 | none found |
| GSM972097 | S-S1+S-S3 | batch03 | C4 | none found |
| GSM972113 | S-S3 | batch03 | C4 | none found |
| GSM972149 | S-S1 | batch04 | C1 | none found |
| GSM972177 | S-S2 | batch04 | C6 | none found |
| GSM972220 | S-S1 | batch05 | C2 | none found |
| GSM972221 | S-S2 | batch05 | C2 | none found |
| GSM972257 | S-S2 | batch05 | C1 | none found |
| GSM972264 | S-S2 | batch05 | C2 | none found |
| GSM972274 | S-S2 | batch05 | C2 | none found |
| GSM972275 | S-S2+S-S3 | batch05 | C4 | none found |
| GSM972283 | S-S2 | batch05 | C2 | none found |
| GSM972295 | S-S2 | batch05 | C3 | none found |
| GSM972305 | S-S1 | batch07 | C2 | none found |
| **GSM972311** | S-S3 | batch07 | C5 | none found — **but see §3.3** |
| GSM972357 | S-S3 | batch08 | C6 | none found |
| GSM972361 | S-S3 | batch08 | C4 | none found |
| GSM972444 | S-S2+S-S3 | batch05 | C4 | none found |
| GSM972445 | S-S2 | batch05 | C2 | none found |
| GSM972467 | S-S3 | batch07 | C4 | none found |
| GSM972479 | S-S3 | batch07 | C4 | none found |
| GSM972499 | S-S1 | batch08 | C2 | none found |

**26 of 29 have NO independent technical evidence of invalidity.** Three have partial evidence, and in all three it is another *statistical* flag rather than a demonstrated technical fault.

### 3.2 Direct technical evidence vs biological hypothesis

**Direct technical evidence — and the key negative result.** The only independent technical measurement available on this matrix is the AFFX 3′/5′ degradation contrast. It **does not discriminate flagged from unflagged samples**:

| | flagged (n=29) | unflagged (n=424) |
|---|---:|---:|
| ACTB 3′/5′ log2 diff, median | 0.792 | 0.772 |
| GAPDH 3′/5′ log2 diff, median | 0.259 | 0.188 |

Moreover, of the **7 arrays** with the largest ACTB contrast (> 2 log2), only **one** (GSM972105) is flagged at all.

> **Conclusion: the one available piece of direct technical evidence does not corroborate the flags.**

**Batch association — evidence, but moderate.** Flag rate by processing batch:

| batch | n | flagged | % |
|---|---:|---:|---:|
| batch05 | 49 | 10 | **20.4** |
| batch03 | 102 | 8 | 7.8 |
| batch08 | 41 | 3 | 7.3 |
| batch01 | 33 | 2 | 6.1 |
| batch07 | 81 | 4 | 4.9 |
| batch04 | 103 | 2 | 1.9 |
| batch02 / batch09 / batch10 | 14 / 15 / 15 | 0 | **0.0** |

χ² = 23.2, df = 8, **p = 0.003**, **Cramér's V = 0.226**. Flags are **batch-associated but not batch-determined** — a moderate association, far weaker than the near-deterministic batch patterning seen for molecular-assay missingness in Phase 6D.

**Balanced across partitions:** discovery 6.8% vs validation 5.1% — flags are not a partition artefact.

**Biological hypothesis — labelled as hypothesis.** Flagged samples span CIT subtypes C1–C6 with no single subtype dominating. CMS4/stromal phenotypes are known to shift global expression distributions for genuinely biological reasons, so a biological contribution is plausible. **Neither the technical nor the biological explanation is established by these data**, and each investigation record says so.

### 3.3 GSM972311 — the clustering singleton

The 2-cluster dendrogram cut splits **452 / 1**, isolating **GSM972311**.

| | |
|---|---|
| Σ(1 − ρ) | **57.72** — vs fence 32.098, and the highest in the cohort |
| per-array median | 4.8096 — **within** the S-S1 fence |
| per-array IQR | 2.8450 — **within** the S-S2 fence |
| PCA z-scores | PC1 −0.48, PC2 −2.19 — **unremarkable** |
| degradation contrast | **no deviation** from the unflagged median |

**It is correlation-distant without being distributionally or degradation-abnormal.** No independent technical evidence of invalidity was found. **It is not excluded**, and the pattern — distant in correlation, normal in every direct technical measure — is noted for the record.

---

## 4. Descriptive findings — 8, none triggering exclusion

| ID | Finding |
|---|---|
| **D-S1** | PCA variance explained PC1–PC5: **9.7%, 8.1%, 6.8%, 5.0%, 3.8%**. No dominant component; no outlier rule applied |
| **D-S2** | Hierarchical clustering 2-cluster split **452 / 1** (GSM972311). No cut-height rule |
| **D-S3** | 9 batches; SD of batch-mean PC1 = **0.041**. **No batch correction applied** |
| **D-S4** | 3′/5′ log2 diff — ACTB median **0.773** (0.15–2.63), GAPDH median **0.190**. **Affymetrix < 3 threshold NOT applied** (scale mismatch) |
| **H-G1** | ✅ non-finite values in the gene matrix: **0** |
| **D-G1** | zero-variance genes in S3: **0**. Frozen per-fold M2 rule governs |
| **D-G2** | low expression: mean < 4.0 → 6,370; < 5.0 → 10,509; < 6.0 → 14,224; sd < 0.3 → 8,756. **No filter applied** |
| **D-G3** | **10 AFFX control probe sets** across ACTB, GAPDH, STAT1 in the frozen median representation — **NOT removed** |
| **D-G4** | mapping: 43,163 unique / 2,619 ambiguous / 8,893 unmapped (frozen Phase 2) |

**On D-S1:** PC1 explains only 9.7%. There is no single dominant technical axis — consistent with a matrix that has already been quantile-normalised and batch-adjusted, and an argument against reading any PCA position as a technical defect.

---

## 5. Outputs

| File | Contents |
|---|---|
| `qc_summary.json` | all 20 criterion results, counts, flags, restriction attestations |
| `qc_criteria_results.csv` | per-criterion table |
| `sample_qc_metrics.csv` | 453 × (median, IQR, distance, flags, non-finite count) |
| `investigation_records.json` | **29 seven-field records** |
| `investigation_summary.csv` | one-line-per-sample view |
| `pca_scores.csv` | PC1–PC5 with dataset and batch |
| `batch_composition.csv` | batch × partition |
| `degradation_ratios.csv` | ACTB / GAPDH 3′/5′ per sample |

---

## 6. Unresolved technical concerns

| # | Concern | Status |
|---|---|---|
| 1 | **H-C1 mis-specification** | **REQUIRES APPROVAL** (§1.1) |
| 2 | 26 of 29 flagged samples have **no independent technical evidence** of invalidity | no exclusion proposed |
| 3 | The only direct technical measure (3′/5′) **does not discriminate** flagged from unflagged | limits what QC on this matrix can establish |
| 4 | Flags are batch-associated (V = 0.226) but **batch correction is prohibited** | documented; not actioned |
| 5 | **GSM972311** is correlation-distant with no supporting technical evidence | documented; not excluded |
| 6 | QC on an already RMA+ComBat matrix is **partly self-fulfilling** (Phase 7B §1.2) | inherent limitation of the approved input |

---

## 7. Decisions — all three APPROVED

All three requested decisions were approved and are recorded in §0 and in
`config/analysis/QC_DECISIONS_RECORD.json`:

1. ✅ **H-C1 amended** — defect, corrected criterion, rationale, approval and the original execution result all recorded (spec revision 3).
2. ✅ **All 453 S3 samples retained** — 0 excluded, 29 investigation records preserved, with the explicit statement that no independent evidence demonstrated technical invalidity.
3. ✅ **batch05 and GSM972311 recorded as QC observations**, explicitly **not** exclusion justifications.

**Next gate:** EDA requires separate approval. Not started.
