# Phase 7B — QC Protocol (APPROVED WITH CLARIFICATIONS)

**Date:** 2026-11-02 · **Revision 2** — seven review clarifications applied.
**Status:** ✅ **APPROVED WITH CLARIFICATIONS** (2026-11-02). ⚠️ **QC NOT EXECUTED.**
**Spec:** `config/analysis/QC_SPECIFICATION.json` · **Tests:** `tests/test_phase7b_qc.py`

**Restrictions observed:** no QC execution, no exclusions, no normalization, no ComBat, no packages installed, no CEL files, no modelling, no outcome analysis, no modification of Phase 2 / Phase 3 / frozen M2 / A-19 artifacts.

---

## 0. APPROVAL AND BINDING CLARIFICATIONS

> ✅ **Specification approved** with seven clarifications, all applied below and enforced by tests.

### C1 — Thresholds are now immutable

**All 20 criteria and their current thresholds are PRESERVED.**

> **Thresholds must NOT be changed after observing QC distributions** without a **separate documented decision and explicit approval**, including a written record of what was observed before the change was proposed.
>
> **A threshold changed after inspection is no longer pre-specified.** That is the whole point of fixing them now.

### C2 — Hard exclusions are narrowed, not broadened

A hard criterion may trigger automatic halt or rejection **only when both** hold:

1. the failure is **deterministic** — no threshold judgement is involved; **and**
2. the sample/data integrity is **demonstrably invalid** — not inferred.

**The 9 hard criteria are fixed and may not be broadened:** H-S1, H-S2, H-S3, H-S4, H-G1, H-C1, H-C2, H-C3, H-C4.

Note that two of the nine do **not** act automatically even so: **H-S3** (duplicate columns) investigates first, and **H-G1** (non-finite gene row) **escalates** rather than removing a gene. Automatic action is confined to cases where invalidity is self-evident.

### C3 — Soft flags cannot exclude. Ever, on their own.

> **No sample may be excluded SOLELY because it crosses a Tukey 1.5 × IQR fence, or appears unusual in PCA, clustering, per-array median, per-array IQR, or between-array distance.**
>
> **Crossing a fence is a reason to look, never a reason to exclude.**

Every soft and descriptive criterion now carries `automatic_exclusion: PROHIBITED` in the machine-readable spec.

### C4 — The exclusion pathway, and the record it requires

```
soft flag
   ↓
individual investigation record  (7 required fields, one per sample)
   ↓
STOP — request explicit approval
   ↓
exclusion only if approved
```

**Required fields in every investigation record:**

| # | Field | Content |
|---|---|---|
| 1 | **sample_id** | the GSM accession |
| 2 | **flagged metric and value** | which criterion fired, the observed value, the fence crossed |
| 3 | **metadata and batch** | partition, ComBat batch, stage, other recorded characteristics |
| 4 | **independent technical evidence** | evidence of invalidity **independent of the flag itself** — **state explicitly if NONE exists** |
| 5 | **biological plausibility** | could this reflect biology rather than technique? (e.g. a CMS4 stromal phenotype) |
| 6 | **effect on cohort composition and event counts** | reported for **transparency only** — **must never inform the decision** |
| 7 | **proposed action and scientific justification** | |

> **Field 6 is the dangerous one.** Reporting the effect on event counts is required for transparency, but using it to decide whether to exclude would be outcome-guided selection. The spec records this prohibition explicitly.

### C5 — AFFX sensitivity analysis: approved, conditional, and non-rescuing

**Approved**, for **M3 only**, triggered **only** if M3 materially implicates ACTB, GAPDH or STAT1. The **frozen primary representation is preserved**, and **both analyses are reported**, clearly labelled.

> **Anti-rescue clause.** The sensitivity analysis must **not** be used to rescue an unfavourable result. If the primary is unfavourable and the sensitivity favourable, **both are reported and the primary remains primary.** Running it selectively, or reporting only the more favourable of the two, is prohibited.

### C6 — Preprocessing decision preserved

**The deposited matrix is used as-is.** No additional normalization, no ComBat or other batch correction, no fRMA execution, no CEL acquisition.

---

## 0. The methodological posture that shapes this protocol

> **I deliberately did NOT compute the QC metric distributions before writing this specification.**

Choosing a numeric cutoff by inspecting the distribution it is meant to screen is circular: the threshold then encodes what the data happens to look like, and "outlier" becomes "whatever I decided to call one." §3 of the brief asks for exactly this discipline.

Consequently, every criterion below is one of three things:

| | |
|---|---|
| a **deterministic** integrity check needing no threshold | e.g. non-finite values, ID matching |
| a **data-adaptive** rule with a documented source | Tukey boxplot-whisker, the arrayQualityMetrics convention |
| **descriptive only**, with no rule, because none is justifiable | PCA, clustering, 3′/5′ ratios |

**Where a defensible constant does not exist, I have not invented one.** Three criteria are explicitly investigation-only, with exclusion requiring separate written approval.

---

## 1. QC scope under the approved input

**Approved input (Phase 7A Rank 1):** the deposited matrix — 54,675 probe sets × 585 samples, log2, RMA + ComBat with an **undocumented fitting population**.

### 1.1 Metrics that CANNOT be computed — and why that is a limitation, not an oversight

| Standard Affymetrix QC metric | Why unavailable |
|---|---|
| **GNUSE** | requires fRMA — **not used under Rank 1** |
| **NUSE / RLE** | require probe-level PLM fits from CEL files |
| **RNA degradation slope** (`affy::AffyRNAdeg`) | requires CEL probe-level data |
| **Percent present / MAS5 calls** | require CEL data |
| **Background, spike-in controls** | not recoverable post-RMA |

> Per the brief: **GNUSE is not applied.** It is meaningless without fRMA, and applying a metric outside the pipeline it was defined for would be worse than omitting it.

### 1.2 What an already-transformed matrix permits

RMA quantile-normalises across arrays and ComBat then adjusts for batch. Two consequences:

1. **Distributional QC is partly self-fulfilling.** Per-array medians are near-constant *because normalisation made them so*. A departure is still informative — it means an array did not normalise as expected — but the metric is less sensitive than on raw data.
2. **Batch-related signal has already been adjusted.** Residual batch structure understates what was originally present (Phase 4 §5). QC can describe what remains; it cannot recover what was removed.

---

## 2. Criteria — 20 across three tiers

| Tier | n | Meaning |
|---|---:|---|
| **HARD EXCLUSION** | **9** | deterministic failure where integrity is **demonstrably invalid**; fixed list, **may not be broadened** (§C2) |
| **SOFT INVESTIGATION** | **3** | pre-specified data-adaptive flag; triggers **investigation only**. **Automatic exclusion PROHIBITED** — requires an investigation record (§C4) and explicit written approval |
| **DESCRIPTIVE** | **8** | reported only; **never independently triggers exclusion** |

### 2.1 Sample-level

| ID | Metric | Rule | Tier | Action |
|---|---|---|---|---|
| **H-S1** | non-finite values in the sample's column | any > 0 | **HARD** | exclude + report *(currently 0)* |
| **S-S1** | per-array **median** log2 | **Tukey**: beyond Q3 + 1.5·IQR of the metric's own cross-array distribution | **SOFT** | **investigate only — exclusion prohibited on this flag alone** |
| **S-S2** | per-array **IQR** | Tukey | **SOFT** | **investigate only — exclusion prohibited on this flag alone** |
| **S-S3** | sum of (1 − Spearman ρ) to every other array | Tukey | **SOFT** | **investigate only — exclusion prohibited on this flag alone** |
| **D-S1** | PCA (PC1–PC5, variance explained) | **no rule** | DESCRIPTIVE | report only |
| **D-S2** | hierarchical clustering (correlation, average linkage) | no cut-height | DESCRIPTIVE | report only |
| **D-S3** | batch composition, PC-vs-batch association | no threshold | DESCRIPTIVE | report only — **batch correction PROHIBITED** |
| **D-S4** | 3′/5′ log2 difference, ACTB and GAPDH | **threshold NOT applied** | DESCRIPTIVE | report only |

**Why Tukey, and why it is legitimate to pre-specify.** The rule flags values beyond Q3 + 1.5·IQR of the metric's own distribution. It is the outlier convention used by **arrayQualityMetrics** (Kauffmann, Gentleman & Huber, *Bioinformatics* 2009;25(3):415–416), which applies boxplot whiskers to per-array summaries and to the row-sums of the between-array distance matrix. Because it is **data-adaptive**, it can be fixed now without my having seen the data — which is precisely the property a pre-specified rule needs.

**S-S2 justification.** IQR captures dispersion the median does not, and IQR-based statistics are reported to perform well on *normalised* Affymetrix data — the relevant case here (Mieczkowski et al., *Bioinformatics* 2014;30(10):1392–1398, IQRray).

**Why PCA and clustering get no rule (D-S1, D-S2).** No literature-justified PCA outlier threshold exists for this setting, and defining one after seeing the scores would be circular. More importantly, **PCA structure in this cohort is expected to reflect biology** — the CMS4 stromal axis dominates (Phase 6C §3.4) — so a PCA outlier is not presumptively technical. They are diagnostics, not criteria.

### 2.2 The 3′/5′ ratio — a documented threshold that does NOT transfer (D-S4)

Phase 7A found the AFFX degradation controls present in the matrix, so the ratio is computable:

| | |
|---|---|
| ACTB | `AFFX-HSAC07/X00351_3_at`, `_5_at`, `_M_at` |
| GAPDH | `AFFX-HUMGAPDH/M33197_3_at`, `_5_at`, `_M_at` |

Affymetrix guidance recommends a **3′/5′ signal ratio < 3** for GAPDH and ACTB as an RNA-degradation indicator.

> **That threshold is explicitly NOT applied here.** It is defined on **MAS5 signal ratios on the linear scale**. Our values are **log2, RMA-summarised and ComBat-adjusted** — the ComBat step having moved each probe by a batch-specific amount that does not cancel in a 3′/5′ contrast. **The quantity is not the one the constant was calibrated for, so the constant does not transfer.**
>
> Computing the ratio and comparing it to 3 would look rigorous and be wrong. It is reported as a **relative** diagnostic across arrays only.

### 2.3 Probe / gene-level

| ID | Metric | Rule | Tier | Effect on frozen artifacts |
|---|---|---|---|---|
| **H-G1** | non-finite values in a gene row | any > 0 | **HARD** | **ESCALATE, do not drop** — removing a gene would reopen a frozen spec |
| **D-G1** | zero-variance gene count | **no cohort-level exclusion** | DESCRIPTIVE | none — the **frozen per-fold M2 rule** governs |
| **D-G2** | per-gene mean/SD distribution; counts below candidate cutoffs | **no filter applied** | DESCRIPTIVE | filtering **PROHIBITED** without approval |
| **D-G3** | AFFX control inventory in the frozen representations | **no removal** | DESCRIPTIVE | **none — removal would reopen a frozen decision** |
| **D-G4** | probe-to-gene mapping summary | no action | DESCRIPTIVE | none |

**Critical design rule: no gene-level QC criterion may remove a gene from a frozen representation.** Zero-variance handling is **already frozen** in the M2 scoring spec — σ = 0 genes are excluded **per fold**, with the exclusion recorded. A cohort-level removal would both duplicate and contradict that. H-G1 therefore **escalates** rather than acting.

### 2.4 Cohort integrity

| ID | Check | Tier | Expected |
|---|---|---|---|
| **H-C1** | sample-ID set equality: matrix ↔ frozen Phase 3 bundle | **HARD — halts QC** | 0 mismatches |
| **H-C2** | discovery/validation membership matches `partition_and_batch.csv` | **HARD — halts** | 0 |
| **H-C3** | S3 membership equals frozen `cohort_membership.csv` (n = 453) | **HARD — halts** | 0 |
| **H-C4** | **GSM972262 excluded** under the frozen complete-case strategy | **HARD — verify** | **exactly 1 exclusion, documented** |
| **H-S2** | sample absent from the frozen manifest | **HARD** | 0 |
| **H-S3** | exactly-duplicated sample columns | **HARD** | 0 expected; investigate if found |
| **H-S4** | non-tumoural sample in the analysis cohort | **HARD** | 0 — S3 excludes the 19 by construction |

**QC cannot add, reassign or redefine samples.** Every cohort check compares against a frozen artifact and **halts** on mismatch rather than reconciling.

---

## 3. Protection against analytical bias

| Prohibited | Status |
|---|---|
| outcome-guided sample exclusion | **prohibited** |
| outcome-guided threshold adjustment | **prohibited** |
| changing S3 based on model performance | **prohibited** |
| modifying frozen Phase 2 representations | **prohibited** |
| modifying either frozen M2 specification | **prohibited** |
| removing genes for weak outcome association | **prohibited** |
| additional batch correction | **prohibited** |
| applying an expression filter without approval | **prohibited** |
| silently removing or replacing AFFX probes | **prohibited** |

**No outcome column is opened by any criterion in this protocol.** QC runs entirely on the expression matrix, the frozen membership tables and the frozen specifications. Any deviation from the frozen cohort or representations **requires explicit approval**.

**Structural guarantee:** the only criteria that can exclude a sample automatically are deterministic integrity failures (H-S1, H-S2, H-S4) and cohort-integrity halts. **Every judgement-based flag is investigation-only.** There is no path from a QC number to an exclusion without a human decision.

---

## 4. AFFX limitation and the proposed sensitivity analysis

**Findings preserved from Phase 7A, not corrected:**

| | |
|---|---|
| AFFX control probe sets in the frozen **median** representation | **10**, across **ACTB, GAPDH, STAT1** |
| Selected as representatives by **rep-IQR** | **`AFFX-HSAC07/X00351_5_at`** (ACTB), **`AFFX-HUMGAPDH/M33197_5_at`** (GAPDH) |
| What those probes are | **5′ RNA-degradation controls** — high IQR because 5′/3′ ratios track degradation, a **technical** quantity |
| Present in the frozen **M2 4-set primary panel**? | **NO** |
| Frozen specifications | **UNCHANGED** |

### Proposed sensitivity analysis — **conditional, and pre-specified now**

> **Warranted conditionally, M3 only.** If an M3 exploratory result materially implicates **ACTB, GAPDH or STAT1**, repeat that analysis with the three AFFX-derived genes removed and report **both results side by side, clearly labelled**.
>
> The removal happens **inside the sensitivity analysis only**. The frozen representations are **not** modified.

**Why conditional rather than unconditional:** removing three genes from a 21,815-gene exploratory model unconditionally would be a change to the analysis set made on technical grounds alone, and M3 is already labelled exploratory with mandatory null-simulation reporting. Pre-specifying the trigger now — before any M3 result exists — prevents it from becoming a post-hoc rescue.

**Why not fix the frozen representation:** it is frozen. Retrofitting would reopen a settled decision, invalidate the Phase 2/3 hashes, and set a precedent that frozen artifacts can be edited when something inconvenient is found later.

---

## 5. Proposed thresholds — summary

| Criterion | Threshold | Source |
|---|---|---|
| non-finite values | **> 0** | deterministic |
| ID / membership mismatches | **any** | frozen artifacts |
| duplicate columns | **exact equality** | deterministic |
| per-array median, IQR, distance | **Tukey: Q3 + 1.5·IQR** | arrayQualityMetrics (Kauffmann 2009) |
| PCA / clustering | **none** | no justifiable constant |
| 3′/5′ ratio | **none** — Affymetrix `< 3` **not applied** | scale mismatch (§2.2) |
| zero / low variance | **none at cohort level** | frozen per-fold M2 rule |
| expression filter | **none** | never approved |

---

## 6. Unresolved decisions requiring approval

| # | Decision |
|---|---|
| 1 | **May any SOFT_INVESTIGATION flag lead to exclusion**, and on what written justification? |
| 2 | **Should a numeric threshold be adopted for S-S1/S-S2/S-S3** once the distributions are seen? This would need a **second approval** and must not be set to achieve a desired exclusion |
| 3 | Is the **AFFX sensitivity analysis** triggered? (conditional on M3 results) |
| 4 | Approve this specification so QC may execute |

---

## 7. Status

| Item | Status |
|---|---|
| QC specification written | ✅ 20 criteria |
| Thresholds justified or explicitly withheld | ✅ |
| QC distributions inspected | ❌ **deliberately not** |
| **QC executed** | ❌ **NO — Phase 7C pending approval** |
| Samples excluded | ❌ none |
| Frozen artifacts modified | ❌ none |
| **Specification approved** | ✅ **APPROVED WITH CLARIFICATIONS** (2026-11-02) — thresholds immutable |
