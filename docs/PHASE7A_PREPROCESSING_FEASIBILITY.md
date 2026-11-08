# Phase 7A — Expression Preprocessing Feasibility Audit

**Date:** 2026-11-02 · **Status:** Audit complete. ✅ **RANK 1 APPROVED 2026-11-02.** **No normalization, no ComBat, no fRMA, no CEL downloads, no package installation, no modelling.**
**Reproduce:** `.venv/bin/python scripts/30_preprocessing_compatibility_audit.py` → `reports/phase7/`

---

## 0. APPROVED DECISION

> ✅ **APPROVED — Rank 1: use the deposited GSE39582 expression matrix as-is** for the current internally-validated study.
>
> | | |
> |---|---|
> | Approved | 2026-11-02T22:05:42.198598+00:00 |
> | Consequence | **no CEL acquisition, no fRMA, no re-normalisation, no re-ComBat** |
> | Frozen artifacts | **unchanged** — the deposited matrix already underlies every Phase 2/3 artifact |
> | Binding limitation | the undocumented ComBat fitting population (§1.3) **forecloses cross-cohort claims** and must appear in the study's limitations |
> | Phase 6B prerequisite #4 | **satisfied** |
> | Phase 6B prerequisite #5 (QC thresholds) | **addressed in Phase 7B** |
>
> **Ranks 2–4 are not adopted.** fRMA/CEL reprocessing remains deferred, not rejected (§4 Rank 2).

---

## 0.1 Headline

| Question | Verdict |
|---|---|
| **fRMA identifier hard gate** | **CONDITIONAL GO** — the CDF ambiguity is **resolved**; the residual loaded-object check needs installation and is still outstanding |
| **Is the deposited matrix usable?** | **Yes for internal analysis** — it is the basis of every frozen Phase 2/3 artifact — **but it carries a material comparability limitation** that forecloses cross-cohort claims |
| **Recommended strategy** | **Rank 1: keep the deposited matrix** for the internally-validated study now authorised. CEL reprocessing is not justified by the current scope |

---

## 1. Deposited expression matrix

### 1.1 Confirmed properties

| Property | Value |
|---|---|
| Dimensions | **54,675 probe sets × 585 samples** |
| Sample coverage | 566 primary tumours + 19 non-tumoural mucosa |
| Scale | **log2** (range −0.399 … 16.112; 2^max ≈ 7.1 × 10⁴, consistent with RMA output) |
| Non-finite values | **0** |
| Documented processing | `RMA normalisation (Bioconductor R affy package) + ComBat (R sva package)` |
| **Fitting population** | **NOT DOCUMENTED** |

### 1.2 Forensic evidence — what can be inferred

The deposit names the methods but not their scope. Three measurements narrow it:

| Observation | Value | Inference |
|---|---:|---|
| SD across arrays of the per-array **median** | **0.034** | quantile normalisation was applied — per-array distributions are near-identical |
| SD of per-array 10th / 90th percentile | 0.037 / 0.102 | same |
| Sorted value vectors **exactly** identical across arrays? | **No** | a post-normalisation adjustment was applied **on top** — consistent with ComBat following RMA |

**Conclusion (inference, not documentation):** the deposited values are **quantile-normalised and then batch-adjusted**. The ordering in the deposit's own description — "RMA … + ComBat" — is corroborated by the data's structure.

### 1.3 The ComBat fitting population — partially narrowed, still unresolved

Phase 4 established that the deposit does not state which arrays entered which run, and that the 566 tumours and 19 non-tumoural samples were submitted ~3 years apart in separate accession blocks, with ArrayExpress's ComBat batch labels covering **only the 566 tumours**.

This audit adds a measurement:

| Comparison | max \|difference\| between sorted value vectors |
|---|---:|
| two **tumour** arrays | 0.525 |
| two **non-tumoural** arrays | **0.190** |
| **tumour vs non-tumoural** | **2.152** |

The cross-group difference is **4.1×** the within-tumour difference. Under a shared quantile normalisation, sorted vectors would be identical *before* ComBat, so residual differences reflect only the subsequent adjustment.

**Inference (not proof):** the 19 non-tumoural arrays were **probably not in the same normalisation/ComBat run** as the 566 tumours. Their much *smaller* within-group spread (0.190 vs 0.525) is consistent with their not having been ComBat-adjusted across 9 batches, as the tumours were.

**What remains unresolved:**
1. Whether RMA was fitted once over the 566 tumours or in sub-batches.
2. Which covariates entered the ComBat model — ArrayExpress exposes `normalizationcombatmod` with values `modcan1`/`modcan2`, but their meaning is undocumented.
3. Whether the 19 normals were processed alone, with the tumours, or in a third run.
4. Whether any array was excluded before deposit.

**None of these is resolvable from the deposit.** Author correspondence is the only route.

### 1.4 Is this a material comparability limitation?

**Yes — and its materiality is scope-dependent.**

| Use | Material? | Why |
|---|---|---|
| **Internal analysis within GSE39582** | **No** | Every sample is on one internally consistent scale. Relative-risk modelling, discrimination and internal validation are unaffected |
| **Discovery → CIT validation** | **Partly** | Both partitions share the scale, but batch and partition are confounded (Cramér's V 0.488, two discovery-only batches), so CIT is an internal holdout at most — as already frozen |
| **Cross-cohort transport** | **YES — disqualifying** | The values sit on a scale defined by an undocumented ComBat fit over this cohort. No external cohort shares it, and no transformation can recover it without knowing the fit |

**This is why Phase 4 concluded NO-GO on external validation and Phase 5 confirmed it on precision grounds. The deposited matrix is adequate for exactly the study that is now authorised, and inadequate for the one that is not.**

---

## 2. fRMA identifier hard gate

### 2.1 Our side — verified against actual identifiers

| Check | Result |
|---|---|
| Expression probe IDs | **54,675** |
| GPL570 annotation probe IDs | **54,675** |
| **Identical sets** | **TRUE** |
| Match standard Affymetrix probe-set grammar | **54,662 / 54,675 (99.98%)** |
| The 13 non-matching | **all `AFFX-` controls** with `/` in the name (e.g. `AFFX-HSAC07/X00351_3_at`) — standard CDF members |
| **Standard-CDF vs ENTREZG discriminator** | see below |

> **Note on a discarded metric.** An earlier draft counted `<digits>_at` IDs as evidence of an Entrez CDF. **That is not a discriminator** — standard Affymetrix IDs have the same form (`202763_at`). The check was replaced with the valid one:

| Discriminator | Our set | A Brainarray ENTREZG CDF |
|---|---:|---|
| probe-set count | **54,675** | ~20,000 |
| non-plain suffixes (`_s_at`, `_x_at`, `_a_at`, …) | **17,261** | **0** — all plain `_at` |
| `AFFX-` controls | **62** | **0** |

Suffix breakdown: `_at` 37,414 · `_s_at` 11,295 · `_x_at` 4,207 · `_a_at` 1,720 · `_g_at` 13 · `_r_at` 12 · `_i_at` 9 · `_f_at` 5.

**Our identifiers are unambiguously standard CDF.**
| Frozen **median** probe set (43,163) ⊂ platform | **TRUE** |
| Frozen **rep-IQR** probe set (21,815) ⊂ platform | **TRUE** |

**Our frozen representations are keyed entirely on standard Affymetrix probe-set identifiers from the GPL570 platform record. There is no Entrez-CDF contamination on our side.**

### 2.2 The CDF ambiguity — RESOLVED

Phase 4 flagged a blocking uncertainty: `hgu133plus2frmavecs`'s DESCRIPTION says it *"was created by frmaTools version 1.19.3 and **hgu133ahsentrezgcdf** version 19.0.0"* — an Entrez custom CDF named inside a plus2 package.

**The package reference manual resolves it.** Verbatim from the manual:

> **`hgu133plus2frmavecs` — Vectors used by fRMA for hgu133plus2.**
> *Usage:* `data(hgu133plus2frmavecs)` · `data(hgu133plus2hsentrezgfrmavecs)`
> *Format:* A list with **6 or 7** elements — `normVec`, `probeVec`, `probeVarWithin`, `probeVarBetween`, `probesetSD`, `medianSE`, **`version` — "the Entrez Gene alternative CDF version *if applicable*"**

**The package ships TWO separate vector objects:**

| Object | CDF |
|---|---|
| **`hgu133plus2frmavecs`** | **standard** |
| `hgu133plus2hsentrezgfrmavecs` | Entrez alternative |

The DESCRIPTION string refers to how the **alternative** object was built. The `version` element exists only *"if applicable"* — i.e. for the altCDF object, which is why the list has "6 **or** 7" elements. **The standard object is standard-CDF keyed**, and with `hgu133plus2cdf` loaded, `frma()` emits the standard 54,675 probe sets.

### 2.3 Verdict: **CONDITIONAL GO**

| | |
|---|---|
| ✅ Our identifiers are 100% standard-CDF | verified against actual IDs |
| ✅ Both frozen probe sets are strict subsets of the platform | verified |
| ✅ The Entrez-CDF ambiguity is explained | reference manual, verbatim |
| ⏳ **Outstanding** | that the **loaded** standard object yields exactly our 54,675 identifiers |

**The residual check requires installing `hgu133plus2frmavecs`, which this phase prohibits.** It is a five-minute check once installation is authorised, and it is the last item on the gate.

**This is not a NO-GO.** A NO-GO would be warranted if our identifiers were non-standard, or if the only available vectors were Entrez-keyed. Neither holds. **No mapping has been forced**: nothing was remapped, substituted or coerced.

### 2.4 NEW FINDING — AFFX control probesets inside the frozen representations

The identifier check surfaced something not previously noticed:

| Gene | AFFX probes in frozen **median** spec |
|---|---|
| **ACTB** | `AFFX-HSAC07/X00351_3_at`, `_5_at`, `_M_at` |
| **GAPDH** | `AFFX-HUMGAPDH/M33197_3_at`, `_5_at`, `_M_at` |
| **STAT1** | `AFFX-HUMISGF3A/M97935_3_at`, `_5_at`, `_MA_at`, `_MB_at` |

**10 AFFX control probe sets across 3 genes** are in the median representation, because GPL570 annotates these controls with gene symbols.

**More consequentially, the rep-IQR criterion selected two of them as representatives:**

| Gene | selected probe |
|---|---|
| ACTB | **`AFFX-HSAC07/X00351_5_at`** |
| GAPDH | **`AFFX-HUMGAPDH/M33197_5_at`** |

These are **5′ RNA-degradation control probe sets**. Their IQR is high precisely because 5′/3′ ratios track RNA degradation — a **technical** quantity. This is a concrete instance of the limitation stated in Phase 2B §3: *"variability is not biological quality… an IQR criterion cannot distinguish technical variability from real biology."*

**Scope of the problem — bounded:**

| | |
|---|---|
| AFFX-derived genes in the frozen **M2 4-set panel** | **NONE** |
| Affected genes | ACTB, GAPDH, STAT1 — available to **M3 only** |
| Frozen specifications | **unchanged** — this is documented, not corrected |

**Recommendation:** record as a known limitation of the M3 exploratory model. Do **not** retrofit the frozen representation — that would reopen a frozen decision. If M3 results implicate ACTB or GAPDH, the AFFX provenance must be reported.

---

## 3. Alternative preprocessing strategies

### 3.1 Requires raw CEL files

| Strategy | Data needed | Reproducible by a third party? | Key limitation |
|---|---|---|---|
| **Uniform RMA** over all eligible CELs | **all** CELs (~566 + any external) | **No** — the scale depends on the exact array set in the run | Adding a cohort later changes every prior value |
| **fRMA** with frozen vectors | CELs, per array | **Yes** — the frozen reference is a published artifact | Reference population composition **undocumented**; *reduced* batch susceptibility, not elimination |
| **GCRMA / MAS5** | CELs | yes/no respectively | Not comparable with the deposited RMA scale; no advantage here |

**Acquisition requirement for any of these: `GSE39582_RAW.tar`, ~566 arrays.** Phase 4's §3.1 six-step inspection (manifest, count, GSM mapping, header-only array-type check, SHA-256, publication reconciliation) would have to pass first. **Currently classified `listed, unverified`.**

### 3.2 Applicable to the deposited matrix

| Strategy | Validity | Assessment |
|---|---|---|
| **Use as deposited** | **Valid for internal analysis** | One internally consistent scale; it is what every frozen Phase 2/3 artifact is built on |
| Re-normalise the deposited matrix | **Not valid** | Values are already quantile-normalised and batch-adjusted. Re-normalising a normalised matrix compounds transformations without a defined target |
| Apply further ComBat | **Not valid, and out of scope** | ComBat has already been applied by an undocumented fit; a second pass on top is uninterpretable. Explicitly prohibited by the current protocol |
| Within-sample rank transformation | **Technically valid, no benefit here** | Would make values robust to between-cohort scale — but there is no external cohort, and it would invalidate the frozen Phase 2 representations |

> **Technical availability is not validity.** Re-normalising or re-ComBat-ing the deposited matrix is trivially executable and scientifically indefensible, because the starting point is already transformed by a procedure whose parameters are unknown.

---

## 4. Ranked recommendation

### Rank 1 — **Use the deposited matrix as-is** ✅ RECOMMENDED

| Criterion | Assessment |
|---|---|
| **Scientific rationale** | The authorised study is a **single-cohort, internally-validated** analysis (Phase 5). Every sample sits on one consistent scale; relative-risk estimation, discrimination and internal validation are valid on it. Reprocessing would buy comparability with cohorts we have decided not to use |
| **Compatibility with frozen Phase 2** | **Perfect — it IS the basis of them.** The probe cache, both frozen representations and the Phase 3 bundle all derive from it. Zero migration |
| **Compatibility with discovery/internal-validation design** | **Adequate with the already-frozen caveat**: CIT is an internal holdout, never external validation |
| **Batch-leakage risk** | **Pre-existing and documented, not introduced.** Batch is confounded with partition (V = 0.488); the 9 batch labels are carried in the Phase 3 bundle for diagnostics. **No batch correction is performed** |
| **Distribution-mismatch risk** | **None internally.** Material only for cross-cohort transport, which is out of scope |
| **Data acquisition** | **None.** Already in hand and hash-verified |
| **Remaining uncertainties** | The ComBat fitting population stays unresolved (§1.3); cross-cohort claims remain foreclosed |

### Rank 2 — **fRMA from CELs** (deferred, not rejected)

Right strategy for a *different* study — one with external validation. It gives a third-party-reproducible common scale and would resolve §1.3 by replacing the undocumented pipeline entirely.

**Why not now:** requires ~566-array acquisition plus Phase 4 §3.1 verification; requires installing the Bioconductor stack; **would invalidate every frozen Phase 2/3 artifact** and force re-freezing the representation, the cohort bundle and M2's resolved memberships; and buys comparability with cohorts Phase 5 ruled out on precision grounds (**48 external events against ~338 required**).

**Revisit if** a cohort with ~340+ compatible stage II–III events appears.

### Rank 3 — **Uniform RMA from CELs**

Same costs as Rank 2 without the third-party reproducibility: the scale depends on the exact array set, so it is not portable. **Dominated by Rank 2** whenever CEL reprocessing is justified at all.

### Rank 4 — **Re-normalising or re-ComBat-ing the deposited matrix** ❌ NOT RECOMMENDED

Technically trivial, scientifically indefensible (§3.2). Listed only to record that it was considered and rejected.

---

## 5. Unresolved questions

| # | Question | Status |
|---|---|---|
| Q1 | ComBat fitting population — which arrays, which covariates | **UNRESOLVED.** `modcan1`/`modcan2` undocumented. Author correspondence only |
| Q2 | Were the 19 non-tumoural arrays in the same run? | **INFERRED: probably not** (4.1× cross-group difference). Not proven |
| Q3 | Was RMA fitted once over 566 tumours or in sub-batches? | **UNRESOLVED** |
| Q4 | Were any arrays excluded before deposit? | **UNRESOLVED** |
| Q5 | Does the **loaded** `hgu133plus2frmavecs` yield exactly 54,675 standard IDs? | **OUTSTANDING** — needs installation |
| Q6 | Composition of the fRMA reference population | **UNDOCUMENTED** on the package page |
| Q7 | Should AFFX-derived ACTB/GAPDH/STAT1 be excluded from M3? | **OPEN** — documented as a limitation; frozen specs unchanged |
| Q8 | Is `GSE39582_RAW.tar` complete and valid? | **UNVERIFIED** — Phase 4 §3.1 not run |

---

## 6. Recommended next decision

> **Adopt Rank 1 — proceed with the deposited matrix — and record the comparability limitation in the study's limitations section.**

This requires **no acquisition, no installation, no reprocessing, and no change to any frozen artifact.**

The two remaining Phase 6B prerequisites then become:

| Prerequisite | Status under Rank 1 |
|---|---|
| Preprocessing strategy decided (#4) | **satisfied by this decision** |
| **QC thresholds pre-specified (#5)** | **still outstanding** — but under Rank 1 the relevant QC is on the *deposited* matrix (distributional outliers, array-level summaries), **not** GNUSE, which requires fRMA |

**If instead Rank 2 is preferred**, the next step is the single cheapest one: authorise installing `hgu133plus2frmavecs` **alone** to discharge Q5 — no CEL files needed — before committing to acquisition.

---

## 7. Restrictions observed

No normalization · no ComBat · no fRMA execution · no CEL downloads · no package installation · no survival or recurrence analysis · no modelling or feature selection · no modification of frozen Phase 2, Phase 3 or M2 artifacts · no commits.
