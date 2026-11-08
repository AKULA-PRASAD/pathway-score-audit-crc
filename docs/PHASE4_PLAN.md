# Phase 4 — External Validation Feasibility (rev. 3)

**Date:** 2026-11-02 · **Revision 3** — superseded in part by the decision record.
**Status:** PLANNING ONLY; **acquisition NO-GO**. No external expression data retrieved. No CEL files downloaded. No R/Bioconductor installed. Phase 3 bundle untouched.

**Method note.** Everything below rests on **GEO sample-level metadata, Bioconductor package indices and published literature**. No expression matrix or CEL file was retrieved. Working files are in a scratchpad outside the project.

---

> ## ⚠️ SUPERSEDED IN PART — see `PHASE4_PREACQUISITION.md` §0 and §6
>
> The pre-acquisition investigation resolved several items this plan left open
> and **changed the outcome**:
>
> | This plan (rev. 2) said | Now established |
> |---|---|
> | GSE14333 ↔ GSE17536 overlap "confirmed, 130 triples" | **131 GSMs; 129 unambiguous pairs + a 2×2 ambiguity**, accessions recorded |
> | GSE17536 ↔ GSE17537 "suspected overlap" | **NOT overlap** — degenerate `t=0.0` key artifacts |
> | Dukes mapping "unresolved" | **resolved in the 129 matched patients only**; the 185 are a **proxy population**, not confirmed AJCC II–III |
> | `DFS_Cens` polarity `[INFERRED]` | **confirmed by cross-deposit comparison** (127/129) |
> | `frma 1.64.0` | **1.62.0** — Bioc 3.23 needs R 4.6; local R is 4.5.2 → **Bioc 3.22** |
> | "Conditional GO, reduced scope" | **❌ NO-GO.** G4 failed on pre-specified criteria |
> | fRMA protocol (§6) | moved to **`PHASE4_FRMA_PILOT.md`** |
>
> **Event-count policy (D3):** <100 events = exploratory only · 100–199 =
> conditional + precision assessment · ≥200 = preferred + precision assessment.
> GSE14333 (**48**) and GSE33113 (**18**) are **exploratory only**, and their
> counts **must not be combined**.
>
> §§1–11 below are retained unaltered as the record of the reasoning at the time.

---

## 0. What changed in this revision

| # | Correction | Outcome |
|---|---|---|
| A | "CEL availability confirmed" was unjustified | **Downgraded to `listed, unverified`** for every cohort, with the required inspection specified (§3) |
| B | Shared investigators treated as near-proof | **GSE39084 reclassified `independence unresolved`.** A patient-level test was then run and found a **confirmed-level overlap of 130 patients between GSE14333 and GSE17536** (§2) |
| C | Total GEO n used as validation n | **Eligible stage II–III counts and event counts computed per cohort** — they are far smaller (§4) |
| D | fRMA recommended prematurely | **Technical comparison prepared; no commitment** (§6) |
| E | Pooling implied | **Cohort-specific validation; pooling replaced by meta-analysis of per-cohort estimates** (§7) |
| F | Ranking implied certainty | **Reissued with `unknown` preserved as unknown; no composite score** (§8) |

---

## 1. Candidate cohorts

| Accession | GEO n | Platform | Deposited processing | Endpoint field |
|---|---:|---|---|---|
| GSE17536 | 177 | GPL570 | RMA | `dfs_event` (cancer recurrence) |
| GSE17537 | 55 | GPL570 | RMA | `dfs_event` |
| GSE14333 | 290 | GPL570 | MAS5.0 | `DFS_Time` / `DFS_Cens` |
| GSE33113 | 96 | GPL570 | MAS5.0 (GCOS) | "meta or recurrence within 3 years" |
| GSE38832 | 122 | GPL570 | MAS5 | `dfs_event` |
| GSE37892 | 130 | GPL570 | RMA-variant | dates: surgery / distant metastasis / last contact |
| GSE39084 | 70 | GPL570 | — | — |
| GSE31595 | 37 | GPL570 | — | `relapse free survival time` (months, labelled) |
| GSE92921 | 59 | GPL570 | MAS5 | `dfs_event` / `dfs_time` (days) |
| GSE29621 | 65 | GPL570 | — | not inspected |
| GSE12945 | 62 | **GPL96** | — | — |
| GSE41258 | 390 | **GPL96** | — | — |

---

## 2. Independence and overlap — evidence-based

**Principle applied: shared investigators or institutions are overlap *warnings*, not proof of duplicate patients.** Where possible, a patient-level test was run: matching exact `(DFS_time, age, sex)` triples across cohorts. Identical continuous follow-up times to 2 decimal places, paired with identical age and sex, are not plausibly coincidental at scale.

### 2.1 Confirmed overlap

| Pair | Evidence | Verdict |
|---|---|---|
| **GSE14333 ↔ GSE17536** | **130 matching `(DFS_time, age, sex)` triples.** 130 of GSE17536's 143 distinct `DFS_Time` values (**90.9%**) also occur in GSE14333. Example matches: `(0.92, 54, M)`, `(2.26, 40, F)`, `(3.32, 83, F)`, `(3.64, 78, M)` | **CONFIRMED OVERLAP** (metadata-level). GSE17536 appears to be largely contained in GSE14333. **They cannot both be used.** |

This was found despite the two cohorts sharing **no authors and no institution** — which is precisely why "no shared authors" cannot be used to infer independence.

### 2.2 Suspected overlap

| Pair | Evidence | Verdict |
|---|---|---|
| GSE17536 ↔ GSE17537 | **3** matching triples (of 171 vs 42 keys). Different institutions (MCC vs VMC) | **SUSPECTED, low.** 3 matches is within coincidence range for age+sex+time. Needs patient-level resolution |
| **GSE39084** ↔ GSE39582 | Contributors **L. Marisa, A. De Reynies**; `@ligue-cancer.net` contact. **No patient-level test possible** — GEO record exposes no comparable endpoint fields | **INDEPENDENCE UNRESOLVED.** Shared investigators are a warning, not proof. Must not be used until resolved, and must not be labelled proven-overlapping either |
| **GSE37892** ↔ GSE39582 | Contributor **Olschwang is a Marisa 2013 co-author**; samples carry a French `center` field (e.g. "CAL (Nice Lacassagne)") consistent with the CIT collection. Records **dates**, not durations, so no triple test was possible | **INDEPENDENCE UNRESOLVED** |

### 2.3 No overlap identified

| Pair | Evidence | Verdict |
|---|---|---|
| GSE38832 ↔ GSE17536 / GSE17537 | **0** shared `dfs_time` values | **No overlap identified** — but **weak**: GSE38832 exposes no age or sex, so only a time-only test was possible |
| GSE38832 ↔ GSE14333 | 1 shared time value | **No overlap identified** (weak, same limitation) |
| GSE14333 ↔ GSE17537 | **0** matching triples | **No overlap identified** |

**"No overlap identified" is not the same as "independent."** For every cohort the definitive test remains genotype concordance from CEL-derived SNP-sensitive probe signals, which has not been performed.

---

## 3. Raw archive vs confirmed CEL availability

**Corrected classification.** GEO listing a `RAW.tar` does **not** confirm usable CEL files.

What was actually verified, from the per-sample GEO records:

| Cohort | Series RAW.tar listed | Per-sample supplementary file | Mapping | **CEL availability** |
|---|---|---|---|---|
| GSE17536 | yes | 177 × `.CEL.gz` | 1:1 GSM | **listed, unverified** |
| GSE17537 | yes | 55 × `.CEL.gz` | 1:1 GSM | **listed, unverified** |
| GSE38832 | yes | 122 × `.CEL.gz` | 1:1 GSM | **listed, unverified** |
| GSE14333 | yes | 290 × `.CEL.gz` | 1:1 GSM | **listed, unverified** |
| GSE33113 | yes | 96 × `.CEL.gz` | 1:1 GSM | **listed, unverified** |
| all others | yes | not inspected | — | **unknown** |

Each inspected cohort lists exactly **one `.CEL.gz` per sample, named by GSM** — stronger than a bare RAW.tar, since it establishes a sample↔file mapping. **It still does not establish that the files are readable, complete, or GPL570 CELs.**

### 3.1 Required archive-manifest inspection (Step 4A2)

Before any cohort is accepted, and **without processing expression values**:

1. **Manifest listing** — enumerate archive contents without extracting payloads; record file name, size, and count.
2. **File count vs sample count** — assert exactly one CEL per GSM, no extras, no missing.
3. **Sample mapping** — every CEL resolves to exactly one GSM; no duplicates.
4. **Header-only integrity read** — read each CEL *header* to confirm: valid CEL magic/version, array type = `HG-U133_Plus_2`, expected row/column dimensions, scan date present.
5. **Checksums** — SHA-256 per file, recorded before use.
6. **Completeness vs publication** — compare array count against the numbers reported in the source publication.

A cohort is reclassified from `listed, unverified` to **`CEL confirmed`** only after all six pass.

---

## 4. Eligible sample and event counts

**Total GEO n is not the eligible validation n.** Computed from sample-level metadata:

| Cohort | GEO n | Stage field | **Eligible for stage II–III** | **Recurrence events** | Non-events | Unknown |
|---|---:|---|---:|---:|---:|---:|
| **GSE14333** | 290 | Dukes | **185** (B+C) | **48** | 137 | 0 |
| **GSE17536** | 177 | AJCC | **114** (II+III) | **31** | 83 | 0 |
| **GSE33113** | 96 | AJCC | **96** (**II only**) | **18** | 72 | 6 |
| **GSE38832** | 122 | AJCC | **74** (II+III) | **≤9** | 65 | 0 |
| **GSE17537** | 55 | AJCC | **34** (II+III) | **3** | 31 | 0 |
| GSE31595 | 37 | TNM | unknown | unknown | — | — |
| GSE92921 | 59 | AJCC | unknown | unknown | — | — |
| GSE29621 | 65 | — | unknown | unknown | — | — |

**This is the most consequential finding of the revision.** Event counts in the *eligible* population are far smaller than the headline sample sizes:

- **GSE17537 (3 events) and GSE38832 (≤9 events) cannot validate anything on their own.**
- Only **GSE14333 (48)** and **GSE17536 (31)** offer meaningful event counts — **and they overlap by 130 patients, so at most one may be used.**
- GSE33113 contributes **18 events, stage II only**.

**Realistic maximum independent external event count: ~48 + 18 = 66**, from one of {GSE14333, GSE17536} plus GSE33113 — and GSE33113 covers only half the target stage range.

### 4.1 GSE14333 `DFS_Cens` polarity — resolved

The field is undocumented. Resolved empirically:

| Evidence | Finding |
|---|---|
| Dukes A with `DFS_Cens = 1` | **89%** |
| Dukes D with `DFS_Cens = 1` | **0%** (all NA) |
| median `DFS_Time`, `Cens = 0` | **19.1** |
| median `DFS_Time`, `Cens = 1` | **47.6** |

If `1` meant *event*, the best-prognosis stage would carry the most events and events would have *longer* follow-up than censored observations. Both are backwards.

**Therefore `DFS_Cens = 1` is CENSORED and `0` is EVENT.** This is an **inference from data pattern, not documentation** — tier `[INFERRED]`, and it must be re-verified against the source publication before use.

---

## 5. Endpoint and stage compatibility

### 5.1 Endpoint constructs are NOT interchangeable

| Construct | Definition | Cohorts |
|---|---|---|
| **First recurrence (RFS)** | time from surgery to first recurrence, any site | **GSE39582 (ours)**, GSE31595 |
| **Disease-free survival (DFS)** | varies; usually recurrence *or* death, but **none of these cohorts state which** | GSE17536, GSE17537, GSE38832, GSE14333 |
| **Metastasis-or-recurrence composite** | explicit composite | GSE33113 |
| **Distant-metastasis-free survival** | distant metastasis only; local recurrence excluded | GSE37892 |

**These are different estimands.** DFS that counts death as an event is not RFS. A metastasis-only endpoint omits local recurrence entirely. **No cohort states its death-without-recurrence rule** — including ours, where it is itself an inference.

### 5.2 Time units

| Cohort | Unit | Basis |
|---|---|---|
| GSE39582 (ours) | **months — PROVISIONAL** | undocumented; six sources exhausted |
| GSE17536, GSE17537, GSE14333 | months (observed max 142.55) | **inferred from range** |
| GSE38832, GSE31595 | months | **explicitly labelled** |
| **GSE33113** | **days** | **established**: max 3599 → 9.9 y; months would give 300 y |
| GSE92921 | days | inferred (2101) |
| GSE37892 | **dates, not durations** | explicit |

### 5.3 Dukes staging — mapping NOT assumed

GSE14333 uses Dukes A/B/C/D. **The common Dukes B → AJCC II and C → AJCC III equivalence is not applied here.** Dukes B spans T3–T4 N0 without the T-substaging AJCC uses, and Dukes C does not distinguish N1 from N2. Without the underlying T and N values — **which GSE14333's GEO record does not provide** — the mapping cannot be verified.

**Reported as unresolved.** If GSE14333 is used, the mapping must be established from the source publication's clinical tables, or the cohort analysed on its own Dukes scale with the difference stated.

### 5.4 GSE33113 scope

**Retained as a possible stage-II-only validation cohort.** It is 100% AJCC stage II and **cannot validate the stage II–III target population.** It may validate the stage II stratum only, and any result must be reported as such.

---

## 6. fRMA vs conventional RMA — technical comparison, no commitment

### 6.1 Verified availability

| Fact | Status |
|---|---|
| `frma` Bioconductor package | **Confirmed** — version 1.64.0, Bioconductor 3.23 |
| **`hgu133plus2frmavecs`** (GPL570 frozen vectors) | **CONFIRMED PRESENT** in the Bioconductor annotation repository, **version 1.5.0** |
| Platform match | `hgu133plus2` **is** HG-U133 Plus 2.0 = GPL570 ✅ |
| Vector provenance | Package page states only: *"created by frmaTools version 1.19.3 and hgu133ahsentrezgcdf version 19.0.0"*. **It does NOT state how many arrays were used, their source, or how batches were chosen** |
| Dependencies | `frma`, `hgu133plus2frmavecs`, `affy`, `oligo`, `preprocessCore`, `frmaTools` |

**Original method reference:** McCall MN, Bolstad BM, Irizarry RA. *Frozen robust multiarray analysis (fRMA).* **Biostatistics** 2010;11(2):242–253. PMID 20097884.

### 6.1.1 Local toolchain status — nothing installed

| Component | Status |
|---|---|
| **R** | **4.5.2 present on this machine, PRE-EXISTING.** Binary dated 2025-10-31, symlink 2026-02-24 — both predate this session. **Not installed by this project.** |
| `frma` | **NOT installed** |
| `hgu133plus2frmavecs` | **NOT installed** |
| `affy` / `oligo` | not checked (not needed until 4C) |

**No R or Bioconductor package was installed during Phase 4 planning**, per the review restriction. That R happens to be available is convenient for Step 4C but changes nothing now; installing the Bioconductor stack remains a Step 4C action requiring approval.

### 6.2 Comparison

| | **Option 1: fRMA** | **Option 2: conventional RMA, one run over all eligible CELs** |
|---|---|---|
| **Principle** | Probe-effect and variance estimates are **precomputed and frozen** from large public databases; each array is normalised to that fixed reference | All arrays quantile-normalised and summarised **together** in a single run |
| **Documented claim** | *"comparable to RMA when the data are analyzed as a single batch and outperforms RMA when analyzing multiple batches"*; *"reduced susceptibility to batch effects"* (McCall 2010) | Standard, universally understood, no external reference needed |
| **Scale** | Common **external** reference — a new cohort can be added later without reprocessing anything | Common **internal** reference — the scale is defined by whichever arrays are in the run |
| **Adding a cohort later** | No reprocessing required | **Every cohort must be reprocessed**; all prior results shift |
| **Probe-level output** | Yes — probe-set summaries on the standard `hgu133plus2` CDF | Yes — same probe-set identifiers |
| **Compatible with our frozen representations?** | **Yes** — both emit the 54,675 GPL570 probe-set IDs our representations are keyed on | **Yes** |
| **Key limitation** | Frozen vectors derive from a reference population whose composition the package page does not document. Colorectal tumour tissue may be under-represented. **Reduced susceptibility to batch effects is not elimination of study effects** | The scale depends on run composition, so it is **not reproducible by a third party** without the identical array set. Mixing cohorts in one quantile normalisation can itself induce cross-cohort structure |
| **Risk** | Inherits an undocumented reference | Fragile to run composition; not portable |

### 6.3 Handling GSE39582 and the 19 non-tumoural samples

**GSE39582 tumour CELs must go through the identical pipeline.** Its deposited values are RMA + ComBat with an **undocumented fitting population** (Phase 3 §1); leaving discovery on that scale while reprocessing external cohorts would guarantee incomparability. Its `GSE39582_RAW.tar` is listed, so this appears feasible — **subject to the same §3.1 inspection**.

**The 19 non-tumoural samples must be processed and reported separately**, never pooled into a tumour normalisation run. They were submitted ~3 years later in a separate accession block, are absent from the ArrayExpress record carrying the ComBat batches, and are not part of any analysis cohort. Including them in a quantile-normalisation run would let non-tumour tissue influence the tumour scale.

### 6.4 What neither option achieves

**Neither fRMA nor uniform RMA guarantees elimination of study effects.** Scanner generation, RNA extraction protocol, tissue handling, collection era (1987–2007 for ours) and laboratory practice leave signal that no probe-level normalisation removes. fRMA's documented claim is *reduced susceptibility*, not removal. **Any cross-cohort performance claim must be accompanied by distributional comparison evidence**, and a validation failure must not be attributed to biology without first excluding processing as the cause.

**No commitment is made.** The decision belongs to Step 4C after §3.1 inspection establishes what is actually available.

---

## 7. Cohort-specific validation strategy

**Expression matrices from different cohorts will NOT be pooled before validation.** Pooling would require a shared normalisation whose validity is exactly what is in question, and would make a single bad cohort contaminate all results.

**Procedure:**

1. Each eligible cohort is assessed **separately**, start to finish, with its own coverage audit, endpoint documentation and performance estimate.
2. Results are reported **per cohort**, with that cohort's endpoint construct, time unit, stage coding and processing provenance attached to the estimate.
3. **Only if** two or more cohorts qualify **and** their endpoint construct and estimand are established as compatible, a **meta-analysis of per-cohort validation estimates** may follow — combining estimates, never combining expression matrices.
4. Cohorts with different endpoint constructs (metastasis-only vs any recurrence) are **not** meta-analysed together. Where constructs differ but are judged comparable, that judgement is stated explicitly with its justification.

---

## 8. Reissued ranking

**No composite score is computed.** Unknown values are shown as `unknown`; a score built from missing evidence would manufacture false precision.

| Cohort | Independence | Eligible II–III | Events | Endpoint | CEL | GPL570 | Processing | Covariates | Class |
|---|---|---:|---:|---|---|---|---|---|---|
| **GSE14333** | no overlap w/ GSE39582 identified; **overlaps GSE17536 (130 pts)** | **185** (Dukes B+C) | **48** | DFS, rule unstated; polarity `[INFERRED]`; **Dukes mapping unresolved** | listed, unverified | ✅ | MAS5 → reprocess | stage, age, sex, location, AdjCTX, AdjXRT | ⚠️ **Conditional** |
| **GSE17536** | no overlap w/ GSE39582 identified; **overlaps GSE14333 (130 pts)** | **114** | **31** | DFS, rule unstated; months inferred | listed, unverified | ✅ | RMA → reprocess | stage, age, sex, ethnicity, grade | ⚠️ **Conditional** |
| **GSE33113** | no overlap identified | **96** (**II only**) | **18** | metastasis-or-recurrence composite; **days** | listed, unverified | ✅ | MAS5 → reprocess | stage, age, sex | ⚠️ **Conditional — stage II only** |
| GSE38832 | no overlap identified (**weak** — no age/sex) | 74 | **≤9** | DFS, months labelled | listed, unverified | ✅ | MAS5 → reprocess | stage only | ⚠️ Conditional — **underpowered alone** |
| GSE17537 | 3 triples vs GSE17536 (suspected, low) | 34 | **3** | DFS, months inferred | listed, unverified | ✅ | RMA → reprocess | stage, age, sex | ⚠️ Conditional — **unusable alone** |
| GSE92921 | unknown | unknown | unknown | DFS, days | unknown | ✅ | MAS5 | stage, mutations, MSI | ⚠️ Conditional — **unassessed** |
| GSE31595 | unknown | unknown | unknown | **RFS, months labelled** — closest match to ours | unknown | ✅ | unknown | stage, age, sex, nodes, chemo | ⚠️ Conditional — n=37 |
| GSE29621 | Moffitt/Yeatman — **unresolved** vs GSE17536 | unknown | unknown | unknown | unknown | ✅ | unknown | unknown | ⚠️ Conditional — **unassessed** |
| **GSE37892** | **UNRESOLVED** (Marisa co-author; French centres) | unknown | unknown | **distant metastasis only** — different construct | listed, unverified | ✅ | RMA-variant | stage, age, sex, centre, dates | ❌ **Not usable until independence resolved** |
| **GSE39084** | **UNRESOLVED** (Marisa, De Reynies contributors) | unknown | unknown | unknown | listed, unverified | ✅ | unknown | unknown | ❌ **Not usable until independence resolved** |
| GSE41258 | unknown | unknown | unknown | unknown | unknown | ❌ **GPL96** | unknown | unknown | ❌ **Unsuitable — platform** |
| GSE12945 | unknown | unknown | unknown | unknown | unknown | ❌ **GPL96** | unknown | unknown | ❌ **Unsuitable — platform** |

GPL96 loss, measured: **10,379/21,815 representative probes present (47.6%)**; median representation — 7,495 genes complete (34.4%), **5,156 partial (23.6%)**, 9,164 absent (42.0%).

### 8.1 Revised candidate set

**GSE14333 and GSE17536 are mutually exclusive.** Choose one:

| Choice | For | Against |
|---|---|---|
| **GSE14333** | 185 eligible, **48 events** — the larger | Dukes mapping unresolved; `DFS_Cens` polarity inferred; MAS5 |
| **GSE17536** | AJCC staging, RMA deposit, cleaner metadata | 114 eligible, **31 events**; likely a subset of GSE14333 |

**Proposed, subject to §3.1 and 4A:** one of {GSE14333, GSE17536} **plus GSE33113 for the stage II stratum only**, with GSE92921 and GSE31595 assessed as possible additions once their eligible counts are known.

**Maximum realistic independent external event count: ~66** — and only ~48 of those span stage II *and* III.

---

## 9. Revised go/no-go criteria

> **SUPERSEDED.** The gate criteria below were adopted, evaluated, and
> **G4 failed**. The decided outcome is **NO-GO for CEL acquisition** —
> see `PHASE4_PREACQUISITION.md` §6. The criteria are retained here as the
> record of what was agreed before the evidence was in.


Phase 4 proceeds to validation **only if all** of the following hold. Each is a gate, not a preference.

| # | Criterion | Threshold |
|---|---|---|
| G1 | **Independence established** for at least one cohort | No confirmed overlap with GSE39582; genotype-concordance check passed |
| G2 | **CEL availability confirmed** | All six §3.1 checks pass for that cohort **and** for GSE39582 |
| G3 | **Endpoint construct documented** | Definition, time origin, unit and death-handling established from the source publication — not inferred |
| G4 | **Eligible events sufficient** | A pre-specified minimum event count, justified by a Riley-style calculation for *validation* precision, agreed **before** results are seen |
| G5 | **Stage mapping resolved** | AJCC II–III identifiable, or the cohort analysed on its own scale with the difference stated |
| G6 | **Uniform processing achieved** | All cohorts incl. GSE39582 through one pipeline, with distributional comparison evidence |

**If G1–G6 cannot all be met, the correct outcome is NO-GO**, and the project is reported as a single-cohort internally-validated study with its limits stated. **That is a legitimate scientific result, not a failure.** Given §4, G4 is the most likely gate to fail.

---

## 10. Unresolved questions

| # | Question | Status |
|---|---|---|
| Q1 | Is GSE39084 genuinely overlapping with GSE39582, or merely same-investigator? | **Unresolved** — no comparable metadata fields for a patient-level test |
| Q2 | Same for GSE37892 | **Unresolved** — records dates, not durations |
| Q3 | Direction and completeness of the GSE14333 ↔ GSE17536 overlap | 130 triples matched; which is the superset not established |
| Q4 | Do the 3 GSE17536 ↔ GSE17537 matches represent real shared patients? | **Unresolved** — within coincidence range |
| Q5 | GSE14333 `DFS_Cens` polarity | `[INFERRED]` from data pattern; must be confirmed against the publication |
| Q6 | Dukes ↔ AJCC mapping for GSE14333 | **Unresolved** — T/N not in the GEO record |
| Q7 | Death-without-recurrence rule in **every** external cohort | **Unresolved** — none states it; ours is itself an inference |
| Q8 | Our own `rfs.delay` time unit | **PROVISIONAL** — makes any shared horizon uninterpretable |
| Q9 | Composition of the `hgu133plus2frmavecs` reference population | **Undocumented** on the package page |
| Q10 | Eligible counts for GSE92921, GSE31595, GSE29621 | **Unknown** — not yet assessed |
| Q11 | Whether uniform reprocessing actually makes cohorts comparable | **Unmeasurable until attempted**; neither method guarantees it |

---

## 11. Boundaries observed

No external expression data retrieved · no CEL files downloaded · no R/Bioconductor packages installed · no batch correction · no normalization · no survival modelling · no outcome-based gene selection · frozen representations unchanged · Phase 3 bundle unmodified · all existing evidence preserved.
