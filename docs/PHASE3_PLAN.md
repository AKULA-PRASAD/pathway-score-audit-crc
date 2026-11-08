# Phase 3 — Implementation Plan (PROPOSED, not executed)

**Date:** 2026-11-02 · **Revision 2** — approved with four corrections, now incorporated.
**Status:** APPROVED FOR IMPLEMENTATION.

**Objective:** connect the approved expression representation to the approved cohort and endpoint decisions, producing a reproducible analysis-ready dataset **without silently resolving any outstanding scientific uncertainty.**

---

## 1. Design principles

| # | Principle | Consequence for the design |
|---|---|---|
| 1 | **Never recompute an approved decision** | Phase 3 *reads* the frozen Phase 1 manifest and the frozen Phase 2B representations. It never re-derives cohort membership or re-selects probes. |
| 2 | **Uncertainty stays visible in the artifact** | The conditional endpoint columns are named `rfs_*_60units`, not `_5y`. The column name itself refuses to assert the unit. |
| 3 | **Alignment by identifier, never by position** | Every join is an explicit `sample_id` reindex with a coverage assertion. No `.values`, no positional slicing. |
| 4 | **Store the superset, derive the subsets** | S3 ⊂ S1 ⊂ S0 and S3 ⊂ S2 ⊂ S0. Storing S0-scoped matrices plus a membership table avoids 8 near-duplicate matrices and removes drift risk between them. |
| 5 | **Every exclusion traceable to a named criterion** | Each excluded sample carries the criterion that removed it. |
| 6 | **No transformation is fitted in Phase 3** | Nothing is centred, scaled, filtered or corrected, so there is no fitted object that *could* leak. Leakage tests still assert this, because the guarantee must survive future edits. |

---

## 2. Proposed workflow

### Step 3.1 — Verify frozen inputs
Re-derive both approved representations from their frozen specifications and assert bit-identity against the existing Phase 2 candidate matrices. Confirms the frozen specs fully determine the matrices.

### Step 3.2 — Assemble cohort membership
Read `reports/phase1/cohort_manifest.csv` verbatim. Emit a membership table: one row per S0 sample, with `in_S0/S1/S2/S3`, `dataset`, and for excluded samples the criterion that removed them. **No eligibility logic is re-implemented.**

### Step 3.3 — Construct endpoint table
Carry the original `rfs.event` / `rfs.delay` through **unchanged**, and add derived columns as *new* fields (§4).

### Step 3.4 — Assemble covariate table
Carry the clinical covariates audited in Phase 1 §4.4. **Carried, not selected** — Phase 3 makes no modelling decision. Variables Phase 1 ruled out of a clinical baseline are carried in a separate `excluded_covariates` table with the documented reason, so nothing disappears silently.

### Step 3.5 — Attach batch labels
Join the nine recovered ComBat batch assignments by `sample_id`. **No correction performed.**

### Step 3.6 — Build the analysis bundle
Write S0-scoped matrices for both representations plus the tables above, aligned on a single canonical sample order.

### Step 3.7 — Verify and freeze
Run the full test suite (§6), write a provenance manifest with input hashes, and freeze.

---

## 3. Cohort construction

Counts below are **expected values**, taken from the approved Phase 1 manifest. Phase 3 asserts them; it does not produce them.

| Cohort | Role | n | events (full) | events (60u)* | disc n/ev | val n/ev |
|---|---|---:|---:|---:|---:|---:|
| **S3** | **PRIMARY** — both criteria applied | **453** | **136** | **130** | 354 / 107 | 99 / 29 |
| S0 | sensitivity — retain all 8 | 461 | 139 | 133 | 359 / 109 | 102 / 30 |
| S1 | sensitivity — Criterion 1 only | 458 | 136 | 130 | 357 / 107 | 101 / 29 |
| S2 | sensitivity — Criterion 2 only | 456 | 139 | 133 | 356 / 109 | 100 / 30 |

\* conditional on the months assumption (§4.2)

**Exclusion traceability** — every excluded sample named, with its criterion:

| Criterion | Rule | n | Samples |
|---|---|---:|---|
| **C1** — recurrence-free at time origin | `rfs.delay == 0 AND rfs.event == 1` | 3 | GSM972062, GSM972206, GSM972402 |
| **C2** — > 0 units of recorded follow-up | `rfs.delay == 0 AND rfs.event == 0` | 5 | GSM972229, GSM972246, GSM972258, GSM972401, GSM972456 |

Base before exclusions: 461 (stage II–III primary tumours with documented RFS).

**Discovery / CIT validation membership is preserved as an explicit column**, never by file partitioning. Source metadata is not modified.

---

## 4. Endpoint construction

### 4.1 Fields

| Column | Source | Rule |
|---|---|---|
| `rfs_event_source` | `rfs.event` | **verbatim copy**, never altered |
| `rfs_time_source` | `rfs.delay` | **verbatim copy**, never altered |
| `rfs_event_full` | derived | `= rfs_event_source` |
| `rfs_time_full` | derived | `= rfs_time_source` |
| `rfs_event_60units` | derived, **conditional** | `0 if rfs_time_source > 60 else rfs_event_source` |
| `rfs_time_60units` | derived, **conditional** | `60 if rfs_time_source > 60 else rfs_time_source` |
| `time_unit_assumption` | constant | `"months — PROVISIONAL, undocumented in any source"` |
| `censoring_interpretation` | constant | `"death without recurrence = censoring — PROVISIONAL, inferred from 49/49 data pattern, not a published rule"` |

**The conditional columns are named `60units`, not `5y`.** Calling them `5y` would encode the unverified months assumption into the artifact's vocabulary. A reader who later discovers the unit is not months finds a column that never claimed otherwise.

### 4.2 Uncertainties carried, not resolved

| Uncertainty | How Phase 3 carries it |
|---|---|
| **Time unit** | Both endpoint definitions stored side by side; conditional columns named `60units`; assumption recorded as a data field in every bundle manifest. **A 60-unit cutoff is conditional on the months assumption, not a confirmed fact.** |
| **Death-without-recurrence censoring** | Recorded as a data field. `os.event` / `os.delay (months)` are **carried into the bundle** so an alternative DFS-style endpoint remains constructible without rebuilding. |
| **Zero-time records** | Not resolved — S0/S1/S2/S3 all materialisable, each sample tagged with its criterion. |
| **CIT446 stage/M conflict** | `tnm.stage` and `tnm.m` both carried verbatim; a `flag_stage_m_conflict` boolean marks it. **Not recoded.** |

### 4.3 Endpoint hierarchy (Correction 3)

| Rank | Endpoint | Status |
|---|---|---|
| **PRIMARY** | `rfs_event_full` / `rfs_time_full` on cohort **S3** | **The primary stored endpoint.** Unconditional — it does not depend on the time-unit assumption. |
| **CONDITIONAL** | `rfs_event_60units` / `rfs_time_60units` | **Conditional on the provisional months assumption.** Stored separately, never overwriting the full endpoint, and labelled conditional in the manifest. Not a confirmed fact. |
| **SEPARATE** | `os_event_source` / `os_time_source` | Overall-survival fields **preserved separately and verbatim**. Carried so an alternative DFS-style endpoint remains constructible without a rebuild. **Not combined with RFS by Phase 3.** |
| **COHORT SENSITIVITY** | S0 / S1 / S2 | Materialisable for every endpoint variant |

**Neither the death-censoring interpretation nor the zero-time question is
resolved by Phase 3.** Both are recorded as data fields; the OS columns are
preserved separately so that an alternative censoring interpretation can be
implemented later without re-deriving anything.

---

## 5. Expression matrices

### 5.1 Representations

| Role | Source matrix | Frozen spec |
|---|---|---|
| **PRIMARY** | `gene_matrix_median.parquet` | `frozen_gene_to_probeset_median.csv` (21,815 genes / 43,163 probes) |
| **SENSITIVITY** | `gene_matrix_representative_iqr_discovery_443.parquet` | `frozen_gene_to_probe_discovery443.csv` (21,815 genes / 21,815 probes) |

Both already exist and are **not rebuilt**; Step 3.1 verifies they reproduce bit-identically from their frozen specs.

### 5.2 Alignment contract

Every alignment is an explicit `sample_id` reindex:
1. canonical sample order = the S0 sample list sorted by the probe cache's column order;
2. `matrix.reindex(columns=cohort_sample_ids)` — label-based, raises on any missing id;
3. assert the resulting column list **equals** the cohort's manifest sample list, in order;
4. assert no `NaN` introduced by the reindex.

**No positional indexing anywhere.** A test asserts that shuffling the matrix columns before alignment yields an identical result (§6, T-07).

### 5.3 Storage

S0-scoped matrices (21,815 × 461, float64, ≈ 80 MB each) plus a membership table. Each cohort materialises as an exact sample-ID subset, verified against its manifest.

**Rationale:** S3/S1/S2 are strict subsets of S0. Storing 8 matrices would duplicate ~630 MB and create four opportunities for them to drift apart. One superset plus a verified accessor cannot drift.

**No outcome-based gene selection, differential-expression filtering or model-driven feature selection occurs.** All 21,815 genes are carried.

---

## 6. Data integrity and provenance

### 6.1 Recorded in every bundle manifest

Input SHA-256 for each file; sample IDs and canonical order; gene IDs and order; dimensions; dtypes; missing-value counts; environment report; and the two provisional-assumption fields from §4.1.

**Approved input hashes (as of this plan):**

| File | SHA-256 (first 16) | Size |
|---|---|---:|
| `data/processed/probe_matrix.parquet` | `693bb18c05b7cb1f` | 319.3 MB |
| `frozen_gene_to_probeset_median.csv` | `08f8715fb7eca5a0` | 1.07 MB |
| `frozen_gene_to_probe_discovery443.csv` | `c386d5f11cc511dc` | 0.79 MB |
| `selection_population_discovery443.txt` | `9490d2264ebd32cf` | 7 KB |
| `reports/phase1/cohort_manifest.csv` | `283e764835045d0c` | 0.11 MB |
| `metadata/GSE39582_sample_metadata.csv` | `37cef7feb767766f` | 0.15 MB |
| `data/external/E-GEOD-39582.sdrf.txt` | `3f373cf1801df6aa` | 0.54 MB |

A mismatch at build time **aborts**; it does not warn.

### 6.2 Preserved untouched

Probe cache · both frozen specs · all six Phase 2 candidate matrices · Phase 1 manifest · scripts 01–14 · `metadata/` · `data/raw/`.

---

## 7. Processing provenance and batch limitations

### 7.1 Processing provenance — qualified wording (Correction 1)

**What is documented.** Every one of the 585 samples carries the identical
`!Sample_data_processing` value:

```
RMA normalisation (Bioconductor R affy package) + ComBat (R sva package )
```

**What is NOT documented: the fitting population.** The deposit nowhere states
which arrays were processed together in any given RMA or ComBat run. The field
names the *methods*; it does not define the *scope*.

**Two observations caution specifically against assuming a single joint run
over all 585 arrays:**

| Observation | Evidence |
|---|---|
| The two sample groups were deposited nearly three years apart | 566 tumour samples submitted **2012-07-23**; 19 non-tumoural samples submitted **2015-05-12** |
| They occupy separate accession blocks | tumours `GSM971957–GSM972522`; non-tumoural `GSM1681353–GSM1681371` |
| The record carrying the ComBat batch assignments covers only the tumours | ArrayExpress E-GEOD-39582 contains **566** samples, not 585 |

**Approved wording, to be used verbatim in every bundle manifest and in all
downstream documentation:**

> All 585 samples are annotated with the same deposited processing
> description, `RMA normalisation (Bioconductor R affy package) + ComBat
> (R sva package)`. **The deposit does not document the fitting population of
> either procedure** — which arrays entered a given RMA or ComBat run is not
> stated. The 566 tumour samples and the 19 non-tumoural samples were
> submitted ~3 years apart in separate accession blocks, and the ArrayExpress
> record carrying the ComBat batch assignments covers only the 566 tumour
> samples. It must therefore **not** be assumed that a single joint
> normalisation or batch-correction was fitted across all 585 arrays.
> Expression values are used as deposited; this project performs no
> normalisation and no batch correction.

Statements of the form "jointly RMA-normalised and ComBat-corrected across all
585 arrays" are **withdrawn** from this plan and must not reappear.

### 7.2 Cohort–batch imbalance (Correction 2)

Nine ComBat batch labels were recovered from the ArrayExpress SDRF. Coverage
of the S3 cohort is **453/453**. **No batch correction is performed, and no
batch covariate is automatically included in any model.** All labels are
preserved and carried into every bundle for later investigation.

**S3 batch composition by partition — the imbalance is severe:**

| batch | discovery | validation | note |
|---|---:|---:|---|
| **batch01** | **33** | **0** | discovery-only |
| **batch02** | **14** | **0** | discovery-only |
| batch03 | 94 | 8 | |
| batch04 | 76 | 27 | |
| batch05 | 45 | 4 | |
| batch07 | 44 | 37 | most balanced |
| batch08 | 37 | 4 | |
| **batch09** | **1** | **14** | validation-weighted 14:1 |
| batch10 | 10 | 5 | |

**Two batches (01, 02) contain discovery samples only; batch09 is weighted
14:1 toward validation.** Partition membership and processing batch are
therefore strongly associated.

**Approved wording for every bundle manifest:**

> Discovery and CIT validation have **distinct sample partitions but
> associated processing history**. Batch composition differs sharply between
> them (batch01 33/0, batch02 14/0, batch09 1/14). **Distinct partitions do
> NOT imply independence of processing-related information.** CIT validation
> is internal validation only. Batch labels are preserved; no batch correction
> is applied and no batch covariate is automatically included.

## 8. Unresolved scientific decisions

Carried into Phase 3 artifacts, **not** resolved by them:

| # | Question | Status | Carried how |
|---|---|---|---|
| U1 | `rfs.delay` time unit | **PROVISIONAL: months.** Six documentary sources exhausted. | Both endpoint variants stored; columns named `60units`; assumption a data field |
| U2 | Death without recurrence = censoring | **PROVISIONAL**, from 49/49 pattern; no published rule | Recorded as a data field; `os.*` carried so DFS stays constructible |
| U3 | Zero-time records (8) | **UNRESOLVED** | All four cohorts materialisable; per-sample criterion tags |
| U4 | CIT446 stage II vs M1 | **UNRESOLVED** | Both fields verbatim + conflict flag; not recoded |
| U5 | `os.delay == 0`: value or sentinel | **UNRESOLVED** | Values carried verbatim |
| U6 | R0 / curative resection status | **UNRESOLVED**, not derivable | Documented limitation |
| U7 | Which aggregation is biologically better | **UNANSWERABLE** — no reference measurement | Both representations shipped in every bundle |
| U8 | Cross-cohort batch behaviour | **UNMEASURED** | Batch labels carried; Phase 4 obligation |
| U9 | External probe availability / annotation version | **UNVERIFIED** for any cohort | Documented limitation |
| U10 | Era effects (1987–2007, 7 centres) | **UNRESOLVED**, no per-sample date | Documented limitation |

---

## 9. Leakage and integrity test plan

Proposed `tests/test_phase3.py`, **36 tests in seven groups**.

### Group A — cohort membership
| ID | Test |
|---|---|
| T-01 | S3 has **exactly 453** samples, **136** events (full), **130** (60units) |
| T-02 | S0/S1/S2 have exactly 461/458/456 samples and 139/136/139 events |
| T-03 | Per-partition counts match §3 exactly (S3: 354/107 disc, 99/29 val) |
| T-04 | S3 ⊂ S1 ⊂ S0 and S3 ⊂ S2 ⊂ S0 as sample-ID sets |
| T-05 | The 8 excluded samples are exactly the named IDs, each tagged with its criterion |
| T-06 | Cohort membership is byte-identical to the Phase 1 manifest (not recomputed) |

### Group B — alignment
| ID | Test |
|---|---|
| T-07 | **Shuffling matrix columns before alignment yields an identical result** — proves label-based, not positional |
| T-08 | Every cohort matrix's column list **equals** its manifest sample list, in order |
| T-09 | Reindex introduces zero NaN |
| T-10 | Both representations share identical gene index and sample order |
| T-11 | A deliberately missing sample id raises, never silently drops |

### Group C — endpoint
| ID | Test |
|---|---|
| T-12 | `rfs_*_source` byte-identical to `metadata/` |
| T-13 | `*_full` equals source exactly |
| T-14 | `*_60units` rule holds for every row; the 6 known reclassified discovery samples are exactly those reclassified |
| T-15 | Validation reclassifies **0** events (its latest recurrence is at 45) |
| T-16 | `time_unit_assumption` and `censoring_interpretation` present and non-empty in every bundle |
| T-17 | No column named `*_5y` or `*_years` exists anywhere |

### Group D — leakage prevention
| ID | Test |
|---|---|
| T-18 | **No validation sample id appears in `selection_population_discovery443.txt`** |
| T-19 | The frozen rep-IQR map reproduces from discovery_443 alone, bit-identically |
| T-20 | Rebuilding the rep-IQR map from discovery_443 **after deleting all validation columns** gives the identical map — proves validation cannot influence it |
| T-21 | Phase 3 fits **no** transformation: no scaler/centring/filter object is produced or persisted |
| T-22 | No endpoint column is referenced by any expression-side module (static grep) |
| T-23 | Gene set is identical across all four cohorts — no cohort-dependent gene filtering |

### Group E — integrity
| ID | Test |
|---|---|
| T-24 | All recorded input hashes match the live files |
| T-25 | Zero non-finite values in every cohort matrix; dtypes are float64 |
| T-26 | Probe cache and all six Phase 2 candidates unchanged |

### Group F — reproducibility
| ID | Test |
|---|---|
| T-27 | Full bundle rebuild is bit-identical (hash of every output) |
| T-28 | Bundle construction is CWD-independent |

### Group G — additional integrity checks (Correction 4)

| ID | Test |
|---|---|
| **T-29** | **The rep-IQR selection procedure cannot access validation columns.** The selection runs to completion on a matrix containing **only the 443 discovery columns**, and produces the frozen map bit-identically. A procedure that required validation data would raise. |
| **T-30** | **Provenance guard.** Attempting to label a selection as `discovery_443` while passing any validation sample id in `selection_columns` raises. The provenance label cannot be mis-attached. |
| **T-31** | All nine batch labels preserved for every covered sample; **no batch covariate appears in `covariates.csv`** and no batch-derived column appears in the expression matrices. |
| **T-32** | **No batch-correction artifact exists anywhere in the bundle** — no fitted ComBat object, no adjusted matrix, no residualised values. Asserted by scanning bundle contents and manifest keys. |
| **T-33** | OS fields (`os_event_source`, `os_time_source`) present and **byte-identical** to `metadata/`, and stored **separately** from the RFS endpoint columns. |
| **T-34** | The manifest records the endpoint hierarchy explicitly: full RFS = `primary`, 60-unit = `conditional`, OS = `separate_preserved`. |
| **T-35** | **Processing-provenance wording check.** The manifest contains the approved qualified statement and contains **no** phrase asserting a joint fit across all 585 arrays (regex-tested against `jointly`, `all 585 arrays`, `across all 585`). |
| **T-36** | The cohort–batch imbalance table is present in the manifest and records batch01 = 33/0, batch02 = 14/0, batch09 = 1/14. |

**On T-20 and T-29** — together these are the strongest available guarantee. Because Phase 3 fits nothing, the only discovery-fitted object in the pipeline is the rep-IQR probe map. T-20 proves the map is unchanged when validation data is deleted; **T-29 proves the selection procedure itself never needs validation columns to run.**

---

## 10. Files

### Created

| Path | Purpose |
|---|---|
| `src/cohort.py` | Read the frozen manifest; materialise a cohort by sample ID |
| `src/endpoint.py` | Endpoint derivation (§4) |
| `src/bundle.py` | Bundle assembly, manifest writing, hash recording |
| `scripts/23_build_analysis_bundle.py` | Entry point |
| `scripts/24_verify_analysis_bundle.py` | Standalone verification |
| `tests/test_phase3.py` | 36 tests (§9) |
| `data/processed/analysis/` | Bundle output (gitignored) |
| `docs/PHASE3_DATASET.md` | Dataset documentation + data dictionary |

### Modified

| Path | Change |
|---|---|
| `.gitignore` | add `data/processed/analysis/` |
| `docs/PHASE3_PLAN.md` | this file, marked executed |

### Explicitly NOT modified

Scripts 01–14 · `metadata/` · `data/raw/` · `reports/phase1/cohort_manifest.csv` · `data/processed/probe_matrix.parquet` · `data/processed/frozen/*` · all six Phase 2 candidate matrices · `requirements.txt`.

### Proposed output layout

```
data/processed/analysis/
├── expression_median_S0.parquet              21,815 x 461  float64  ~80 MB
├── expression_repiqr_disc443_S0.parquet      21,815 x 461  float64  ~80 MB
├── cohort_membership.csv                     461 rows, in_S0/S1/S2/S3 + criterion
├── endpoint.csv                              461 rows, source + full + 60units
├── covariates.csv                            461 rows, carried clinical variables
├── covariates_excluded.csv                   variables ruled out + reason
├── partition_and_batch.csv                   sample_id, dataset, combat_batch
└── analysis_bundle_manifest.json             hashes, dims, dtypes, assumptions
```

---

## 11. Execution estimate and dependencies

| Step | Estimate |
|---|---|
| 3.1 verify frozen representations (rebuild both from spec) | ~3–4 min |
| 3.2–3.5 assemble tables | < 10 s |
| 3.6 write bundle (~160 MB) | ~30 s |
| 3.7 verification + full test suite | ~2–3 min |
| **Total** | **≈ 7–9 minutes** |

**Dependencies: none new.** pandas 2.3.3, numpy 2.5.4, pyarrow 26.0.0, pytest 9.2.0 are installed. `statsmodels` and `scikit-survival` remain deferred — Phase 3 performs no modelling and does not need them.

**Disk:** ≈ 160 MB added; 132 GB free.

---

## 12. Boundaries

Phase 3 as planned performs **no** survival modelling, **no** outcome-based feature selection, **no** differential expression, **no** normalization, **no** batch correction, **no** external validation, and resolves **none** of U1–U10. It assembles, aligns, documents and verifies.
