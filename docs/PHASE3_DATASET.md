# Phase 3 — Analysis Dataset Documentation

**Date:** 2026-11-02 · **Status:** Built and verified. 57/57 tests pass; deterministic rebuild confirmed.
**Build:** `.venv/bin/python scripts/23_build_analysis_bundle.py`
**Verify:** `.venv/bin/python scripts/24_verify_analysis_bundle.py` · `.venv/bin/python -m pytest`

---

## 1. Processing provenance — read before use

> All 585 samples are annotated with the same deposited processing description,
> `RMA normalisation (Bioconductor R affy package) + ComBat (R sva package)`.
> **The deposit does not document the fitting population of either procedure** —
> which arrays entered a given RMA or ComBat run is not stated. The 566 tumour
> samples and the 19 non-tumoural samples were submitted ~3 years apart
> (2012-07-23 and 2015-05-12) in separate accession blocks, and the ArrayExpress
> record carrying the ComBat batch assignments covers only the 566 tumour
> samples. It must therefore **not** be assumed that a single joint
> normalisation or batch-correction was fitted across every array. Expression
> values are used as deposited; this project performs no normalisation and no
> batch correction.

> Discovery and CIT validation have **distinct sample partitions but associated
> processing history**. Batch composition differs sharply between them
> (batch01 33/0, batch02 14/0, batch09 1/14). **Distinct partitions do NOT imply
> independence of processing-related information.** CIT validation is internal
> validation only. Batch labels are preserved; no batch correction is applied and
> no batch covariate is automatically included.

---

## 2. Bundle contents

`data/processed/analysis/` (gitignored; reproducible from frozen inputs)

| File | Shape | SHA-256 (16) |
|---|---|---|
| `expression_median_S0.parquet` | 21,815 × 461 float64 | `a1ca41807c534584` |
| `expression_repiqr_disc443_S0.parquet` | 21,815 × 461 float64 | `e636e259b4feabfa` |
| `cohort_membership.csv` | 461 rows | `bd3671f8df64c594` |
| `endpoint.csv` | 461 rows | `802a68575a34a3b6` |
| `covariates.csv` | 461 × 16 | `47b2d9829a4fccec` |
| `covariates_excluded.csv` | 4 rows | `2b25a1765b89adf2` |
| `partition_and_batch.csv` | 461 rows | `57a45bab8efb81f9` |
| `exclusions.csv` | 8 rows | `41f7f2919416e1d2` |
| `analysis_bundle_manifest.json` | — | provenance record |

**Storage design:** S3 ⊂ S1 ⊂ S0 and S3 ⊂ S2 ⊂ S0, so matrices are stored **S0-scoped** and each cohort materialises as an exact sample-ID subset. One superset cannot drift the way four near-duplicate copies can.

---

## 3. Data dictionary

### `endpoint.csv`

| Column | Status | Definition |
|---|---|---|
| `rfs_event_source` | **verbatim** | `rfs.event` as deposited, never altered |
| `rfs_time_source` | **verbatim** | `rfs.delay` as deposited, never altered |
| `os_event_source` | **verbatim, separate** | `os.event` — preserved, **not** combined with RFS |
| `os_time_source` | **verbatim, separate** | `os.delay (months)` |
| `rfs_event_full` | **PRIMARY** | `= rfs_event_source` |
| `rfs_time_full` | **PRIMARY** | `= rfs_time_source` |
| `rfs_event_60units` | **CONDITIONAL** | `0 if rfs_time_source > 60 else source` |
| `rfs_time_60units` | **CONDITIONAL** | `60 if rfs_time_source > 60 else source` |
| `reclassified_by_horizon` | derived | true for the 6 discovery events beyond 60 units |

**Why `60units` and not `5y`.** Naming these columns `5y` would assert the months assumption in the artifact's own vocabulary. No source documents the unit of `rfs.delay`. A reader who later learns the unit is not months finds a column that never claimed otherwise.

### `cohort_membership.csv`

`sample_id` · `in_S0` / `in_S1` / `in_S2` / `in_S3` · `exclusion_criterion`

### `covariates.csv`

Carried (**not selected** — Phase 3 makes no modelling decision): age at diagnosis, sex, TNM stage/T/N/M, tumour location, adjuvant chemotherapy, MMR, CIMP, CIN, KRAS, BRAF, TP53, plus `flag_stage_m_conflict`.

`flag_stage_m_conflict` is true for **GSM972402 (CIT446)** only — `tnm.m = M1` with `tnm.stage = 2`. **Both fields are carried verbatim; neither is recoded.** The conflict is unresolved.

### `covariates_excluded.csv`

| Variable | Reason |
|---|---|
| `cit.molecularsubtype` | derived **from** the expression data |
| `rna.extraction.method` | technical; a batch diagnostic, never a predictor |
| `chemotherapy.adjuvant.type` | 69% missing, 5 levels on 51 events |
| `dependancy_sample` | tumour/normal pairing key |

### `partition_and_batch.csv`

`sample_id` · `dataset` (discovery / validation) · `combat_batch` (9 levels, 461/461 coverage) · `combat_mod`

---

## 4. Cohorts

| Cohort | Role | n | events (full) | events (60u)* | disc n/ev | val n/ev |
|---|---|---:|---:|---:|---:|---:|
| **S3** | **PRIMARY** | **453** | **136** | **130** | 354 / 107 | 99 / 29 |
| S0 | sensitivity — retain all 8 | 461 | 139 | 133 | 359 / 109 | 102 / 30 |
| S1 | sensitivity — C1 only | 458 | 136 | 130 | 357 / 107 | 101 / 29 |
| S2 | sensitivity — C2 only | 456 | 139 | 133 | 356 / 109 | 100 / 30 |

\* conditional on the provisional months assumption

**Exclusions** (`exclusions.csv`), each traceable to a named criterion:

| Criterion | Rule | Samples |
|---|---|---|
| **C1** — recurrence-free at origin | `rfs.delay == 0 AND rfs.event == 1` | GSM972062, GSM972206, GSM972402 |
| **C2** — >0 units of follow-up | `rfs.delay == 0 AND rfs.event == 0` | GSM972229, GSM972246, GSM972258, GSM972401, GSM972456 |

---

## 5. Representations

| Role | File | Frozen spec | Selection population |
|---|---|---|---|
| **PRIMARY** | `expression_median_S0.parquet` | `frozen_gene_to_probeset_median.csv` | none — median has no selection step |
| **SENSITIVITY** | `expression_repiqr_disc443_S0.parquet` | `frozen_gene_to_probe_discovery443.csv` | **discovery_443**; validation never used |

`median` is a **provisional, uniform aggregation rule — not claimed biologically optimal.** `mean` is historical evidence only; the α-stratified hybrid is **not approved and not implemented**.

---

## 6. S3 batch composition by partition

| batch | discovery | validation | |
|---|---:|---:|---|
| batch01 | 33 | **0** | discovery-only |
| batch02 | 14 | **0** | discovery-only |
| batch03 | 94 | 8 | |
| batch04 | 76 | 27 | |
| batch05 | 45 | 4 | |
| batch07 | 44 | 37 | most balanced |
| batch08 | 37 | 4 | |
| batch09 | **1** | **14** | validation-weighted 14:1 |
| batch10 | 10 | 5 | |

---

## 7. Unresolved — carried, not resolved

U1 `rfs.delay` time unit (PROVISIONAL: months) · U2 death-without-recurrence censoring (PROVISIONAL) · U3 zero-time records · U4 CIT446 stage/M · U5 `os.delay = 0` sentinel · U6 R0 status · U7 which aggregation is biologically better (**unanswerable** — no reference measurement) · U8 cross-cohort batch behaviour · U9 external probe availability / annotation version · U10 era effects.

---

## 8. What Phase 3 establishes, and what it does not

**Establishes:** computational correctness, exact cohort/event counts, label-based alignment, leakage prevention, deterministic reproducibility.

**Does not establish:** biological validity. The processing provenance is incompletely documented (§1), no normalization or batch correction has been applied, and no aggregation strategy has been shown more accurate than another.

**Not performed:** survival modelling · feature selection · differential expression · normalization · batch correction · external validation.
