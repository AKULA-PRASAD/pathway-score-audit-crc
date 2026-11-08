# Phase 4 — Pre-Acquisition Investigation and Decision Record

**Date:** 2026-11-02 · **Revision 2** — conditional approval; five decisions applied.
**Status:** Investigation complete; **CEL acquisition NOT authorized**. **No external expression data retrieved. No CEL files downloaded. No packages installed.** Phase 3 bundle untouched.

## 0. DECISION RECORD

Five decisions applied on review. They govern every section below.

### D1 — Patient overlap

| | |
|---|---|
| **GSE14333** | Retained as **external EXPLORATORY candidate** |
| **GSE17536** | **EXCLUDED from independent validation** — established patient overlap with GSE14333 |
| Overlap evidence | **129 unambiguous 1:1 matches** and **2 ambiguous matches per side** preserved **separately** (§1.1) |
| Subtraction | **Prohibited.** The 159-sample remainder is **not** a newly independent cohort and must never be presented as one |
| GSE17536 ↔ GSE17537 | Non-overlap finding **preserved with its limitations** (§1.4) |

### D2 — Endpoint and stage

| | |
|---|---|
| GSE14333 `DFS_Cens` | **0 = event, 1 = censored** — supported by cross-deposit comparison (§2.1) |
| Death without recurrence | **UNRESOLVED** — not stated in any external cohort |
| Absolute time unit | **UNRESOLVED** — agreement between cohorts established; absolute unit is not |
| The 185 Dukes B/C cases | **"candidate stage II–III proxy population"** — **NOT individually confirmed AJCC stage II–III** |
| Dukes/AJCC concordance | **demonstrated in the 129 matched patients only**; not a property of the remaining 156 |
| GSE33113 | **separate stage-II-only, composite-endpoint cohort** |

### D3 — Event-count policy (pre-specified, before any result exists)

Operational **screening** thresholds for confirmatory external validation — **not universal statistical guarantees**:

| Events | Classification |
|---:|---|
| **< 100** | **Exploratory only** |
| **100–199** | **Conditional** — formal precision assessment required |
| **≥ 200** | **Preferred** — still subject to formal precision assessment |

The literature requires **precision-based sample-size assessment**. For a time-to-event prediction model this means a **simulation-based assessment targeting calibration-slope precision** and the other relevant performance estimates — not a threshold lookup.

**Current consequence:**

| Cohort | Events | Classification |
|---|---:|---|
| GSE14333 | **48** | **Exploratory only** |
| GSE33113 | **18** | **Exploratory only, separate estimand** |

**Their event counts must NOT be combined.** 48 + 18 = 66 is not a meaningful quantity: the two cohorts carry different endpoint constructs, different time units and different stage coverage.

> **The current evidence does not support confirmatory external validation.**

### D4 — fRMA pilot

**Preparation approved; execution NOT authorized.** The protocol is a separate document: **`docs/PHASE4_FRMA_PILOT.md`**. Proposed environment R 4.5.2 + Bioconductor 3.22. **No packages installed. No CEL files downloaded — including for the hard gate.**

### D5 — Acquisition

> **The ~971-array CEL acquisition is NOT authorized.** No external expression data will be downloaded.

Full acquisition remains unauthorized until **all four** hold:

1. fRMA **or an alternative uniform preprocessing strategy** is technically validated;
2. a suitable external cohort with **compatible endpoint definitions and adequate event information** is identified;
3. a **formal precision-based sample-size assessment** supports the intended validation claims;
4. the **exact external estimand and target population** are defined.

GSE14333 and GSE33113 remain **exploratory candidates**.

---

**Method:** GEO sample-level metadata, Europe PMC bibliographic records, Bioconductor package indices, and cross-cohort metadata matching. Evidence tables at `reports/phase4/`.

---

## 1. GSE14333 ↔ GSE17536 overlap — RESOLVED

### 1.1 Sample-level overlap table

Matching key: exact `(DFS_time, age, sex)`. GSE14333 records specimen identifiers (`T5266A1`); GSE17536 records `MCC Patient N`.

| Quantity | Value |
|---|---:|
| Matching rows | 133 |
| Distinct GSE14333 GSMs matched | **131** |
| Distinct GSE17536 GSMs matched | **131** |
| **Unambiguous 1:1 pairs** | **129** |
| Ambiguous (one-to-many) | 2 on each side |

**Preserved separately, per D1:**

| File | Contents |
|---|---|
| `reports/phase4/overlap_GSE14333_GSE17536_unambiguous.csv` | **129** one-to-one pairs, both GSM accessions |
| `reports/phase4/overlap_GSE14333_GSE17536_ambiguous.csv` | **4 rows** — a 2×2 ambiguity, preserved apart |
| `reports/phase4/overlap_GSE14333_GSE17536_all.csv` | all 133 rows |

**The 2 ambiguous matches per side.** GSE14333 `GSM358396` (T2511B) and `GSM358397` (T5162A1) both match GSE17536 `GSM437139` (MCC Patient 47) and `GSM437169` (MCC Patient 77) — all four share `(22.38, 65, F)` and all are Dukes B / AJCC 2. **Many-to-many on the key; not resolvable to individual patients from metadata alone.** They are counted in the 131 matched GSMs but excluded from the 129 confirmed pairs.

Sample of the recovered accession mapping:

| GSE14333 GSM | specimen | Dukes | GSE17536 GSM | title | AJCC | DFS_time | age | sex |
|---|---|---|---|---|---|---:|---:|---|
| GSM358341 | T5266A1 | A | GSM437114 | MCC Patient 22 | 1 | 3.64 | 78 | M |
| GSM358343 | T1856A | A | GSM437103 | MCC Patient 11 | 1 | 16.47 | 80 | F |
| GSM358344 | T4919A1 | A | GSM437113 | MCC Patient 21 | 1 | 19.75 | 58 | M |
| GSM358347 | T866D | A | GSM437098 | MCC Patient 6 | 1 | 27.02 | 63 | F |
| … | | | | | | | | |

### 1.2 Three independent corroborations

The match is not a single-key coincidence. Three further quantities agree across the 129 pairs:

| Corroboration | Result |
|---|---|
| **Stage** | Dukes A↔AJCC 1 (23), B↔2 (51), C↔3 (55). **Zero off-diagonal** |
| **Time values** | max \|ΔDFS_time\| across 129 pairs = **0.000000** |
| **Recurrence status** | **127/129 (98.4%) concordant** with GSE17536's explicitly labelled `dfs_event` |

### 1.3 Literature context

| | |
|---|---|
| **GSE14333** | Jorissen et al., *Clin Cancer Res* 2009;15:7642–51, DOI 10.1158/1078-0432.CCR-09-1431, PMC2920750. Abstract: *"A multisite single-platform microarray study was done on **553 colorectal cancers**"* — GSE14333 (290) is a subset of a larger multisite collection |
| **GSE17536 / GSE17537** | Smith et al., *Gastroenterology* 2010;138:958–68, DOI 10.1053/j.gastro.2009.11.005, PMC3388775. Abstract: *"55 patients with colorectal cancer from Vanderbilt Medical Center (VMC) as the training dataset and **177 patients from the Moffitt Cancer Center** as the independent dataset"* |

**Full text could not be retrieved.** PubMed, PMC and the Europe PMC `fullTextXML` endpoint were all blocked (CAPTCHA / HTTP 500); both papers are non-open-access. Abstracts were retrieved via Europe PMC. **The Methods sections remain unread**, so the *reason* for the overlap (shared multisite recruitment vs. direct reuse) is not established. This is recorded as an open item, not glossed.

### 1.4 Classification

| Pair | Evidence | Category |
|---|---|---|
| **GSE14333 ↔ GSE17536** | 131 GSMs, 129 unambiguous, 3 independent corroborations | **CONFIRMED OVERLAP** |
| **GSE17536 ↔ GSE17537** | 4 raw matches — **all with `DFS_time = 0.0`**, a degenerate key; stage *disagrees* (2 vs 4); one GSE17537 sample matched two GSE17536 samples | **NOT OVERLAP.** Artifacts of a degenerate key. Previously "suspected"; now **resolved as coincidence** |
| GSE14333 ↔ GSE17537 | 0 matches | **No overlap identified** — but GSE17537 is small and the test has limited power. **Not proof of independence** |

### 1.5 Resolution — one cohort whole, no subtraction

**Subtraction is prohibited (D1).** Removing the 131 matched samples from GSE14333 would leave 159 samples / 77 Dukes B+C / **21 events**. That remainder:

- is **not a newly independent cohort** — it is an arbitrary residue of a defined study population;
- has **no publishable provenance** — no publication describes it;
- is **worse powered** than either source cohort.

**It must never be presented as independent.** The figures are recorded here only to document why the approach was rejected.

**Resolution applied:** use **GSE14333 whole** as the external exploratory candidate; **exclude GSE17536 entirely** from independent validation. GSE14333 is a single defined study population, independent of GSE39582, and carries the larger event count.

| | Eligible | Events | Status |
|---|---:|---:|---|
| **GSE14333 whole** | 185 (Dukes B+C) | **48** | **exploratory candidate** |
| ~~GSE17536~~ | 114 | 31 | **excluded — patient overlap** |
| ~~GSE14333 minus overlap~~ | 77 | 21 | **rejected — not an independent cohort** |

## 2. Endpoint compatibility

### 2.1 GSE14333 — `DFS_Cens` polarity CONFIRMED

The field is undocumented in GEO, and the publication Methods could not be retrieved. It was instead resolved by **cross-cohort label comparison against the 129 matched patients**, using GSE17536's explicitly labelled field:

| GSE14333 `DFS_Cens` | GSE17536 `no recurrence` | GSE17536 `recurrence` |
|---|---:|---:|
| **0** | 0 | **26** |
| **1** | **101** | 2 |

**127/129 concordant (98.4%). `DFS_Cens = 0` is EVENT; `1` is CENSORED.** This is a cross-deposit label validation, materially stronger than the data-pattern argument used in rev. 2 — but it is still **not a statement from the publication**, and two pairs disagree (same patient, different recurrence status in the two deposits), which is itself a data-quality signal.

### 2.2 Time units

| Test | Result |
|---|---|
| \|DFS_time(GSE14333) − dfs_time(GSE17536)\| over 129 pairs | **max 0.000000** |

Both cohorts use the **identical unit**. Observed maximum 142.55 is consistent with months. **This establishes agreement between the two, not an absolute unit** — neither deposit labels it. Our own `rfs.delay` unit also remains provisional, so a shared horizon is still not interpretable.

### 2.3 Endpoint constructs — not harmonised

| Cohort | Construct | Death w/o recurrence | Status |
|---|---|---|---|
| **GSE39582 (ours)** | first recurrence from surgery | censoring (**inferred**, 49/49) | reference |
| **GSE14333** | "disease-free survival" | **NOT STATED** | Methods unread |
| **GSE17536** | "cancer recurrence"; abstract says *"recurrence-free survival"* | **NOT STATED** | Methods unread |
| **GSE33113** | **"meta or recurrence within 3 years"** — explicit composite | **NOT STATED** | — |

**No external cohort documents its death-without-recurrence rule.** GSE17536's abstract uses "recurrence-free survival", which is encouraging for compatibility with our estimand but is an abstract phrase, not a Methods definition.

**GSE33113 is preserved as a separate stage-II-only estimand.** Its endpoint is a metastasis-or-recurrence composite with a 3-year framing, its time is in **days** (max 3599 ≈ 9.9 y), and it is 100% AJCC stage II. **No harmonisation with our RFS is proposed.**

---

## 3. Stage compatibility — Dukes mapping RESOLVED

The 129 matched patients carry **both** staging systems, giving a direct empirical mapping:

| | AJCC 1 | AJCC 2 | AJCC 3 |
|---|---:|---:|---:|
| **Dukes A** | **23** | 0 | 0 |
| **Dukes B** | 0 | **51** | 0 |
| **Dukes C** | 0 | 0 | **55** |

**Zero off-diagonal.** Dukes A = AJCC I, **B = II, C = III**, exactly, for every matched patient.

> **D2 constraint.** This concordance is **demonstrated in the 129 matched patients ONLY**. It is not a property of the other 156 GSE14333 samples, for which no AJCC stage exists in any source.
>
> Accordingly the 185 Dukes B/C cases are designated the **"candidate stage II–III proxy population"** — **NOT individually confirmed AJCC stage II–III**. Every downstream artifact, table and claim must use that wording.

**Scope of this result, stated honestly:** it is established on 129 patients who appear in both deposits. It is strong evidence that the two coding systems were applied equivalently *in this collection*, and it removes the rev. 2 "unresolved" status for GSE14333's B/C → II/III mapping. It is **not** a general claim that Dukes and AJCC are universally equivalent (they are not — Dukes B does not distinguish T3 from T4, and Dukes C does not distinguish N1 from N2). Per-patient T and N are **not** in GSE14333's GEO record, so sub-stage stratification remains impossible.

---

## 4. fRMA feasibility

**Moved to a separate document per D4: `docs/PHASE4_FRMA_PILOT.md`.**

Summary of what was established here:

| Finding | Status |
|---|---|
| **Bioconductor 3.23 requires R 4.6**; local R is **4.5.2** → correct pairing is **Bioconductor 3.22** | established |
| On Bioc 3.22: `frma` **1.62.0** (not 1.64.0 as cited in the plan rev. 2), `hgu133plus2frmavecs` **1.5.0**, `hgu133plus2cdf` 2.18.0, `affy` 1.88.0, `oligo` 1.74.0, `arrayQualityMetrics` 3.66.0 | established |
| The vectors package names **`hgu133ahsentrezgcdf` 19.0.0** — an hgu133**a** Entrez-gene custom CDF inside an hgu133**plus2** package | **BLOCKING, unresolved** |
| fRMA's documented claim is *"reduced susceptibility to batch effects"* — **not elimination** | established |
| Reference-population composition of the frozen vectors | **undocumented** |

**Execution is NOT authorized.** The hard gate requires **no CEL files**.

## 5. Updated candidate ranking

| Cohort | Independence | Eligible II–III | Events | Endpoint | Stage | CEL | Class |
|---|---|---:|---:|---|---|---|---|
| **GSE14333** | independent of GSE39582; **supersedes GSE17536** | **185** | **48** | DFS; polarity **confirmed**; unit = GSE17536's | **Dukes↔AJCC resolved** | listed, unverified | ⚠️ **Conditional — lead candidate** |
| **GSE33113** | no overlap identified | **96 (II only)** | **18** | composite, **days** — separate estimand | AJCC II | listed, unverified | ⚠️ **Conditional — stage II stratum only** |
| ~~GSE17536~~ | **131 patients shared with GSE14333** | 114 | 31 | — | — | — | ❌ **Excluded — redundant with GSE14333** |
| GSE38832 | no overlap identified (weak) | 74 | **≤9** | DFS, months | AJCC | listed, unverified | ❌ **Underpowered** |
| GSE17537 | not overlapping GSE17536 (resolved) | 34 | **3** | DFS | AJCC | listed, unverified | ❌ **Unusable** |
| GSE92921 | unknown | unknown | unknown | DFS, days | AJCC | unknown | ⚠️ Unassessed |
| GSE31595 | unknown | unknown | unknown | **RFS, months labelled** | TNM | unknown | ⚠️ Unassessed, n=37 |
| GSE29621 | unresolved vs GSE17536 | unknown | unknown | unknown | unknown | unknown | ⚠️ Unassessed |
| GSE37892, GSE39084 | **UNRESOLVED** (CIT investigators) | — | — | — | — | — | ❌ Not usable |
| GSE41258, GSE12945 | — | — | — | — | — | — | ❌ GPL96 |

### 5.1 Defensible counts after overlap exclusion

| | Eligible | Events |
|---|---:|---:|
| **GSE14333** (stage II+III) | **185** | **48** |
| **GSE33113** (stage II only) | **96** | **18** |
| **Total** | **281** | **66** |
| **Spanning stage II *and* III** | **185** | **48** |

---

## 6. GO / NO-GO — decided

### 6.1 Revised gate table

| Gate | Criterion | Status |
|---|---|---|
| **G1** | Independence established | ✅ **MET (metadata-level)** for GSE14333 and GSE33113. Genotype concordance outstanding |
| **G2** | CEL availability confirmed | ⏳ **Not started** — not required while acquisition is unauthorized |
| **G3** | Endpoint construct documented | ⚠️ **PARTIAL.** Polarity resolved; **death-handling UNRESOLVED**; **absolute time unit UNRESOLVED**; Methods unread |
| **G4** | Sufficient events | ❌ **FAILED.** GSE14333 **48**, GSE33113 **18** — both **< 100 → exploratory only** (D3). Counts must not be combined |
| **G5** | Stage mapping resolved | ⚠️ **PARTIAL.** Concordance shown in the **129 matched patients only**; the 185 are a **proxy population**, not confirmed AJCC II–III |
| **G6** | Uniform processing validated | ⏳ **Blocked** on the fRMA hard gate (`PHASE4_FRMA_PILOT.md`) |

### 6.2 Decision

> ## ❌ NO-GO for CEL acquisition
>
> **The ~971-array acquisition is NOT authorized. No external expression data will be downloaded.**
>
> **G4 has failed on pre-specified criteria.** Both candidate cohorts fall below 100 events and are **exploratory only**. The current evidence **does not support confirmatory external validation.**

### 6.3 Status of the candidates

| Cohort | Status |
|---|---|
| **GSE14333** | **exploratory candidate** — 185 candidate stage II–III proxy population, 48 events |
| **GSE33113** | **exploratory candidate, separate estimand** — 96 stage II, composite endpoint, days, 18 events |
| GSE17536 | **excluded** — patient overlap with GSE14333 |
| GSE17537, GSE38832 | not pursued — 3 and ≤9 events |
| GSE37892, GSE39084 | not pursued — independence unresolved |
| GSE41258, GSE12945 | not pursued — GPL96 |

### 6.4 Four conditions to unlock acquisition (D5)

| # | Condition | Current status |
|---|---|---|
| 1 | Uniform preprocessing strategy **technically validated** (fRMA or alternative) | **blocked** on the hard gate |
| 2 | Suitable external cohort with **compatible endpoint definitions and adequate event information** | **not identified** — best candidate is exploratory-only |
| 3 | **Formal precision-based sample-size assessment** supporting the intended claims | **not performed** |
| 4 | **Exact external estimand and target population defined** | **not defined** |

**All four must hold. None currently does.**

### 6.5 What this means for the project

A NO-GO on confirmatory external validation is **a result, not a failure**. The defensible paths from here:

1. **Single-cohort internally-validated study** with explicitly stated limits — the honest default.
2. **Exploratory external assessment** on GSE14333 and/or GSE33113, reported *as exploratory*, never as confirmatory validation, and with each cohort's estimand attached.
3. **Expand the candidate search** — GSE92921, GSE31595 and GSE29621 remain unassessed, and other platforms/consortia were not searched. Unlikely to reach 100+ independent events on GPL570, but not excluded.

A precision-based assessment (D3) should be run **before** committing to any of these, since it determines what claim the available events can support.

## 7. Remaining unresolved questions

| # | Question | Status |
|---|---|---|
| Q1 | GSE39084, GSE37892 independence | **Unresolved** — no comparable metadata for a patient-level test |
| Q2 | *Why* GSE14333 and GSE17536 overlap (shared recruitment vs reuse) | **Unresolved** — Methods unread; all full-text routes blocked |
| Q3 | 2 of 129 matched pairs disagree on recurrence status | **Unresolved** — same patient, different label across deposits |
| Q4 | Death-without-recurrence rule in **any** external cohort | **Unresolved** — none states it; ours is an inference |
| Q5 | Absolute time unit in GSE14333/GSE17536 | **Agreement established, absolute unit not** — neither deposit labels it |
| Q6 | `hgu133plus2frmavecs` CDF keying | **BLOCKING** — pilot step 2 |
| Q7 | `hgu133plus2frmavecs` reference population composition | **Undocumented** |
| Q8 | GNUSE exclusion threshold | **Must be pre-specified with a documented source** |
| Q9 | GSE92921, GSE31595, GSE29621 eligible counts | **Unknown — unassessed** |
| Q10 | Our own `rfs.delay` unit | **PROVISIONAL** |
| Q11 | T/N sub-staging for GSE14333 | **Unavailable** in the GEO record |

---

## 8. Boundaries observed

No external expression data retrieved · no CEL files downloaded · no packages installed · no batch correction · no normalization · no modelling · Phase 3 bundle unmodified · frozen representations unchanged.
