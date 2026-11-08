# Phase 4 — fRMA Technical Pilot Protocol

**Date:** 2026-11-02 · **Status:** **PREPARED, NOT EXECUTED.** Execution requires separate authorization (D4).

> **Standing restrictions for this protocol**
> - **No packages installed** until execution is authorized.
> - **No CEL files downloaded — including for the hard gate.** The gate deliberately requires none.
> - No normalization, no batch correction, no modelling.
> - If identifiers prove incompatible: **STOP. Do not remap or substitute identifiers to force compatibility.**

---

## 1. Purpose

Determine whether frozen RMA can place GSE39582 and any external GPL570 cohort on a common scale **while preserving the exact probe-set identifiers that the two frozen gene-level representations depend on**.

The two representations are keyed on GPL570 probe-set IDs:

| Representation | Requires |
|---|---:|
| median (primary) | **43,163** probe-set IDs grouped into 21,815 genes |
| representative-IQR (sensitivity) | **21,815** specific probe-set IDs |

If fRMA output does not carry these identifiers, fRMA is unusable for this project.

---

## 2. Proposed environment — not installed

| Component | Version | Status |
|---|---|---|
| R | **4.5.2** | **pre-existing** on this machine (binary 2025-10-31); not installed by this project |
| Bioconductor | **3.22** (released 2025-10-30) | **not installed** |
| `frma` | **1.62.0** | not installed |
| `hgu133plus2frmavecs` | **1.5.0** | not installed |
| `hgu133plus2cdf` | 2.18.0 | not installed |
| `affy` | 1.88.0 | not installed |
| `oligo` | 1.74.0 | not installed |
| `preprocessCore` | 1.72.0 | not installed |
| `affyPLM` | 1.86.0 | not installed |
| `arrayQualityMetrics` | 3.66.0 | not installed |

**Version correction on record:** Bioconductor **3.23 requires R 4.6**. The plan rev. 2 cited `frma 1.64.0`, which belongs to 3.23. On R 4.5.2 the installable version is **1.62.0**.

Installation, when authorized, targets a **project-local R library** so the system R is untouched, with every version pinned and hashed.

---

## 3. HARD GATE — identifier compatibility (no CEL files required)

**This gate runs on the frozen-vector package alone.** It is cheap, requires no expression data, and decides whether the rest of the protocol happens at all.

### 3.1 The uncertainty being resolved

The `hgu133plus2frmavecs` package page states, in full:

> *"This package was created by frmaTools version 1.19.3 and **hgu133ahsentrezgcdf version 19.0.0**."*

This names an **hgu133a Entrez-gene custom CDF** inside an **hgu133plus2** package. Two readings:

| Reading | Consequence |
|---|---|
| **(a)** Boilerplate artifact in the DESCRIPTION; vectors are keyed on standard `hgu133plus2cdf` probe sets | fRMA is compatible — proceed |
| **(b)** Vectors genuinely built on a custom Entrez-gene CDF | fRMA emits **custom gene-level identifiers**, not the 54,675 Affymetrix probe sets. **Incompatible — stop** |

Neither the package page nor the frma vignette resolves this. **It must be checked against the actual object, not inferred.**

### 3.2 Gate checks

Inspect the **actual frozen-vector object**, not documentation:

| # | Check | Pass criterion |
|---|---|---|
| G-1 | Load `hgu133plus2frmavecs`; enumerate the object's components | object loads; components enumerable |
| G-2 | **Vector lengths** | the probe-level vectors have the length implied by the standard GPL570 CDF |
| G-3 | **Vector names / feature identifiers** | names resolve to standard Affymetrix probe-set IDs (`1007_s_at` form), **not** custom Entrez IDs (`100_at` / `..._at` ENTREZG form) |
| G-4 | **Probe-set count after summarisation** | **exactly 54,675** |
| G-5 | **Set equality against our frozen probe list** | identifier set is **equal** to the 54,675 in `data/processed/probe_matrix.parquet` — not a subset, not a superset |
| G-6 | **Does the custom CDF change feature definitions?** Compare the CDF the vectors reference against `hgu133plus2cdf` 2.18.0: probe-set membership, probe counts per set | **identical feature definitions**; any difference means a probe set with the same name measures different probes |
| G-7 | **Coverage of both frozen representations** | all **43,163** median-grouping probes and all **21,815** representative-IQR probes present |

### 3.3 Gate outcome

| Result | Action |
|---|---|
| **All of G-1…G-7 pass** | fRMA is identifier-compatible. Proceed to §4 **only with separate authorization** |
| **Any check fails** | **STOP.** Record which check failed and the observed values. **Do NOT remap, substitute, subset or coerce identifiers.** The processing strategy reverts to conventional RMA (§5), which uses `hgu133plus2cdf` directly and has no such ambiguity |

> **Explicit prohibition.** Forcing compatibility by remapping identifiers would silently redefine what every gene in our frozen representations measures. That is exactly the failure mode the frozen representations exist to prevent.

---

## 4. Pilot steps — only if the gate passes AND execution is authorized

**Requires a small CEL sample (~20 arrays). Not authorized.**

| # | Step | Pass criterion |
|---|---|---|
| P-1 | Acquire ~20 CELs: 10 GSE39582 tumour + 10 GSE14333, with per-file SHA-256 and header-only validation (array type `HG-U133_Plus_2`, dimensions, scan date) | all headers valid |
| P-2 | Run `frma()`; take `exprs(object)` | `rownames(exprs)` set-equal to our 54,675 probe-set IDs |
| P-3 | **Array QC**: `GNUSE(object, type="stats")` → median, IQR, 95th, 99th per array | values recorded per array |
| P-4 | Run conventional RMA (`affy::rma`) on the **identical** 20 CELs | same 54,675 identifiers |
| P-5 | **Compare** fRMA vs RMA: per-probe correlation, per-array distribution, rank agreement | quantified; no claim of equivalence asserted |
| P-6 | Apply **both frozen representations** to each output | both apply with **zero substitutions**; gene counts = 21,815 |
| P-7 | Cross-cohort distributional comparison on the pilot arrays | recorded as evidence, not as proof of comparability |

### 4.1 Array-level failure criteria — threshold NOT yet set

`GNUSE` is the documented array-quality metric; the vignette's worked example shows median values around 1.08–1.16.

> **No exclusion threshold is adopted here.** A GNUSE cut-off must be **pre-specified with a documented source** before any array is excluded. I do not have such a source, and inventing a cut-off post hoc would let array exclusion be tuned to the result.

Until a sourced threshold exists, QC values are **recorded and reported, not acted on**.

---

## 5. Comparator — conventional RMA

Uniform `affy::rma` across all eligible CELs in a single run.

| | fRMA | conventional RMA |
|---|---|---|
| Scale reference | **external**, frozen | **internal**, defined by the run's array set |
| Adding a cohort later | no reprocessing | **everything reprocessed**; prior results shift |
| Third-party reproducible | yes — frozen reference is a published artifact | **no** — requires the identical array set |
| CDF ambiguity | **the §3 blocking question** | none — uses `hgu133plus2cdf` directly |
| Non-tumoural arrays | structurally safe to separate (each array independent of others) | **unsafe to mix** — quantile normalisation spans the run |

**No preference is committed.** If the gate fails, conventional RMA is the fallback by default.

---

## 6. Non-tumoural arrays

The 19 GSE39582 non-tumoural arrays are **processed separately and never included in any tumour normalisation run**.

- Under **fRMA** each array normalises to the frozen external reference, so separation is structurally guaranteed — no array can influence another.
- Under **conventional RMA** it is not: quantile normalisation is computed across the run, so including normal mucosa would shift the tumour scale.

They were submitted ~3 years after the tumours in a separate accession block, are absent from the ArrayExpress record carrying the ComBat batches, and belong to no analysis cohort.

---

## 7. Freezing the procedure for future cohorts

On success, the pipeline is frozen as a reusable specification recording: R version · Bioconductor release · every package version with hash · frozen-vector package version · the CDF actually used · the exact call and arguments · the identifier set emitted · per-array QC values.

A third party must be able to apply the identical procedure to a future independent cohort **without reprocessing anything that already exists** — which is the property fRMA offers and conventional RMA does not.

---

## 8. What this pilot cannot establish

| | |
|---|---|
| **Batch-effect elimination** | fRMA's documented claim is *"reduced susceptibility to batch effects"* and down-weighting of probes that "appear to be batchy" — **not elimination** (McCall, Bolstad & Irizarry, *Biostatistics* 2010;11:242–253) |
| **Reference adequacy** | the composition of the `hgu133plus2frmavecs` reference population is **undocumented**; whether colorectal tumour tissue is adequately represented cannot be checked |
| **Cross-cohort comparability** | a pilot on ~20 arrays produces evidence, not proof. Any validation claim must carry distributional comparison evidence alongside it |
| **Whether validation is worth doing** | that is settled by the event-count policy (D3), not by this protocol. Both candidate cohorts are currently **exploratory only** |

---

## 9. Authorization status

| Item | Status |
|---|---|
| Protocol prepared | ✅ |
| Package installation | ❌ **not authorized** |
| Hard gate execution (§3, no CELs) | ❌ **not authorized** |
| Pilot execution (§4, ~20 CELs) | ❌ **not authorized** |
| Full acquisition (~971 arrays) | ❌ **NO-GO** (D5) |
