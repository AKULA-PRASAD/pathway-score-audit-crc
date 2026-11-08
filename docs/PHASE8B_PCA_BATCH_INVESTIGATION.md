# Phase 8B — Outcome-Blind PCA and Batch-Structure Investigation

**Version:** 1.0 · **Status:** executed · **Date:** 2026-11-02
**Script:** `scripts/33_phase8b_pca_batch.py` · **Tests:** `tests/test_phase8b_pca_batch.py`
**Outputs:** `reports/phase8b/` — 15 CSV + 1 JSON + 3 figures
**Cohort:** S3, n = 453, **all retained** (Phase 7C decision D2, unchanged)
**Matrix:** deposited GSE39582 values **as-is** — no normalization, no ComBat, no batch correction

---

## 1. Methods and exact inputs

### 1.1 Input artifacts

| Artifact | SHA-256 (16) | Role |
|---|---|---|
| `data/processed/analysis/expression_median_S0.parquet` | `a1ca41807c534584` | frozen median gene-level representation |
| `data/processed/analysis/cohort_membership.csv` | `bd3671f8df64c594` | S3 definition (453) |
| `data/processed/analysis/partition_and_batch.csv` | `57a45bab8efb81f9` | `combat_batch`, `combat_mod`, partition |
| `data/processed/analysis/covariates.csv` | `47b2d9829a4fccec` | clinical / molecular annotation |
| `metadata/GSE39582_sample_metadata.csv` | `37cef7feb767766f` | provenance fields (**outcome columns dropped at read**) |
| `config/pathways/m2_panel_4set.json` | `cb4de071e2f64a6b` | M2 primary (rev. 2, read-only) |
| `config/pathways/m2_panel_5set.json` | `10dae4874a918160` | M2 sensitivity (rev. 2, read-only) |
| `reports/phase7c/*` | — | QC metrics, degradation ratios (read-only) |
| `reports/phase8a/*` | — | PCA scores, sex control (read-only, for reproduction) |

`data/processed/analysis/endpoint.csv` exists and is **never opened**.

### 1.2 PCA method

Gene-wise centring, **no scaling**. Phase 8B decomposes the 453 × 453 **sample covariance**
(`numpy.linalg.eigh`) where Phase 8A took a thin **SVD** of the 21,815 × 453 centred matrix —
a deliberately different algorithm, so that an error in one path cannot cancel in the other.
Scores are rescaled by `√(ev·(n_genes−1))` onto Phase 8A's SVD scale, so every PC value in this
report is directly comparable with the Phase 8A report. Sign convention: largest-magnitude
loading forced positive. No RNG anywhere in the decomposition.

### 1.3 Statistical approach

- **Effect sizes are primary.** With n = 453, small p-values are cheap. Categorical associations use Kruskal–Wallis with ε²; continuous use Spearman ρ. All 360 PC × annotation tests carry **Benjamini–Hochberg q-values**.
- **Stratified permutation test (new in 8B).** To ask whether batch structure exists *beyond* MMR composition, batch labels are permuted **within MMR strata** (10,000 permutations, seed 20261102). The null preserves the batch × MMR contingency, so rejection means structure that composition cannot explain. This keeps all 453 samples, unlike a single-stratum analysis.
- **Detection sensitivity (new in 8B).** A known shift is injected into batch05's PC3 and the stratified test re-run (300 simulated datasets × 500 permutations per effect size) to establish what the test *could* have found.

### 1.4 Outcome-blind enforcement — two layers

The sample metadata file itself carries `rfs.event`, `rfs.delay`, `os.event`,
`os.delay_(months)`. Phase 8A loaded that file whole (using only `mmr.status`); Phase 8B
tightens this:

| Layer | Mechanism |
|---|---|
| **File** | `_guard_outcome_files()` patches `pandas.read_csv/read_parquet/read_table` to raise `PermissionError` on any path containing `endpoint.csv`, for the whole run |
| **Column** | `load_metadata_outcome_blind()` drops all four outcome columns at read time and **asserts** they are absent before returning |

Tests prove both (§7): the file guard must actually raise; the metadata loader must return
identical results when the underlying outcome columns are **corrupted**; and a static scan
confirms no outcome column name appears anywhere in the script outside the ban-list constant.

---

## 2. Reproduction of Phase 8A

| Check | Result |
|---|---|
| Max \|difference\| in variance explained, PC1–20 | **4.9 × 10⁻⁹** |
| Min **signed** correlation of PC scores, PC1–20 | **0.999999999999999** |
| Max relative difference in PC scores, PC1–20 | **5.3 × 10⁻⁹** |

Variance explained PC1–5: **10.32 %, 7.68 %, 6.36 %, 5.92 %, 3.21 %** — identical to Phase 8A.

The residual ~5 × 10⁻⁹ is the float64 floor: forming the sample covariance squares the condition
number. The signed correlations (exact to 10⁻¹²) are the strict check, and they also confirm the
sign convention reproduced.

**Phase 8A is reproduced.** The `combat_mod` ↔ MMR concordance also reproduces exactly: **0
discordant** among the 410 S3 samples with MMR recorded, and **0 among 519** SDRF tumour arrays.
*(Phase 8A's 519 is the SDRF denominator; the 461-row S0 bundle yields 417, also 0 discordant.
Both are reported here so the two figures are never conflated.)*

---

## 3. Variance structure

| | |
|---|---|
| PCs for 50 % / 80 % / 90 % of variance | **17 / 127 / 237** |
| Horn parallel analysis (25 gene-permuted replicates) | **48** PCs exceed the permuted 95th percentile |
| Broken-stick criterion | 97 PCs |

The spectrum is shallow and **high-dimensional**: even by the most conservative of the three
criteria, tens of components carry joint structure. No single technical gradient dominates.

Horn's criterion bounds *how many* components carry joint structure. It says **nothing** about
whether that structure is biological or technical.

---

## 4. PC associations with technical and provenance variables

360 tests (PC1–20 × 18 variables), BH-adjusted. Strongest by effect size:

| Annotation | PC | ε² / \|ρ\| | q (BH) | New in 8B |
|---|---|---|---|---|
| `cit.molecularsubtype` | PC3 | **0.689** | 6.0e-63 | ✱ |
| `cit.molecularsubtype` | PC1 | 0.591 | 7.8e-54 | ✱ |
| `cit.molecularsubtype` | PC4 | 0.287 | 5.8e-25 | ✱ |
| CIMP status | PC3 | 0.207 | 1.9e-17 | |
| MMR status | PC3 | 0.199 | 9.3e-18 | |
| CIN status | PC3 | 0.198 | 5.6e-16 | |
| `combat_mod` | PC3 | 0.195 | 3.4e-19 | |
| MMR status | PC1 | 0.183 | 1.6e-16 | |
| `combat_mod` | PC1 | 0.164 | 2.1e-16 | |
| **ComBat batch** | **PC3** | **0.061** | 2.8e-05 | |
| **RNA extraction method** | **PC3** | **0.052** | **3.0e-04** | **✱** |
| **RNA extraction method** | **PC1** | **0.032** | **1.0e-02** | **✱** |

> ⚠ **`cit.molecularsubtype` is circular and is not evidence.** C1–C6 are the depositors' own
> **expression-derived** clustering of this same array data. Large PC associations are therefore
> near-tautological. It is included only to show where the published subtypes sit in this PCA.

---

## 5. Batch structure within MMR strata

### 5.1 Cell sizes — the limiting factor

| batch | 01 | 02 | 03 | 04 | 05 | 07 | 08 | 09 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| `modcan1` (dMMR) | 11 | **0** | 5 | **0** | 31 | 13 | 2 | 1 | 2 |
| `modcan2` (pMMR) | 22 | 14 | 97 | 103 | 18 | 68 | 39 | 14 | 13 |

The dMMR stratum has **64 samples across 6 usable batches, minimum cell 2**, and is entirely
absent from two batches. Any within-dMMR test is severely underpowered.

### 5.2 Per-stratum tests — **this extends Phase 8A, and qualifies it**

Phase 8A tested the **pMMR stratum only** and reported the association vanishing. Phase 8B tests
both strata:

| Stratifier | Stratum | PC | n | batches | min cell | ε² | p |
|---|---|---|---|---:|---:|---|---|
| `combat_mod` | modcan2 (pMMR) | PC1 | 388 | 9 | 13 | −0.018 | 0.996 |
| `combat_mod` | modcan2 (pMMR) | PC3 | 388 | 9 | 13 | −0.012 | 0.894 |
| `combat_mod` | **modcan1 (dMMR)** | **PC3** | 64 | 6 | **2** | **+0.130** | **0.028** |
| `combat_mod` | modcan1 (dMMR) | PC1 | 64 | 6 | 2 | −0.009 | 0.485 |
| `mmr.status` recorded | pMMR | PC3 | 350 | 9 | 13 | −0.010 | 0.814 |
| `mmr.status` recorded | **dMMR** | **PC3** | 59 | 5 | **2** | **+0.102** | **0.050** |

**The dMMR stratum shows a nominal PC3 × batch association that Phase 8A did not test.** It is
uncorrected for multiplicity, rests on minimum cell size 2, and is not confirmed by the
better-powered conditional test below. It is reported because it exists, **not** because it is
established.

### 5.3 Stratified permutation test — all 453 retained

Batch labels permuted **within** MMR strata, 10,000 permutations:

| Stratifier | PC | n | observed H | null p95 | p |
|---|---|---:|---|---|---|
| `combat_mod` | PC1 | 453 | 20.86 | 42.11 | 0.862 |
| `combat_mod` | PC2 | 453 | 17.43 | 27.05 | 0.393 |
| `combat_mod` | PC3 | 453 | 34.91 | 47.20 | 0.381 |
| `mmr.status` recorded | PC3 | 410 | 39.96 | 48.21 | 0.226 |

Conditioning on MMR composition while keeping every sample, **no PC1–PC5 location association
with batch survives**.

### 5.4 Detection sensitivity — what this test *could* have found

| Injected batch05 PC3 shift | 0.0 SD | 0.1 | 0.2 | **0.3** | 0.5 | 0.75 | 1.0 |
|---|---|---|---|---|---|---|---|
| Power at α = 0.05 | 0.000 | 0.000 | 0.000 | **0.973** | 1.000 | 1.000 | 1.000 |

The test reaches 80 % power at a shift of about **0.3 SD of PC3 (≈ 5.5 PC3 units)**. A batch
location effect smaller than roughly 0.2–0.3 SD would **not** have been detected.

> **This is not evidence that batch effects are absent.** ComBat was applied to this deposit
> before release, so between-batch variance was removed *by construction* and cannot be recovered
> from these data. The null results above are bounded statements about what remains detectable in
> an already-corrected matrix — nothing more.

---

## 6. batch05 PC3 decomposition

Bootstrap, 5,000 resamples:

| `combat_mod` | group | n | mean PC3 | 95 % CI |
|---|---|---:|---|---|
| modcan1 (dMMR) | **batch05** | 31 | **+22.39** | [+18.21, +26.44] |
| modcan1 (dMMR) | other batches | 34 | **+18.15** | [+13.27, +22.83] |
| modcan2 (pMMR) | **batch05** | 18 | **−2.55** | [−9.25, +4.84] |
| modcan2 (pMMR) | other batches | 370 | **−3.42** | [−5.07, −1.74] |

batch05 overall mean PC3 **+13.23** against a cohort mean of 0.00; batch05 is **63.3 % dMMR**
against **14.3 %** cohort-wide.

**Direct observation.** Within each MMR stratum, batch05's mean PC3 is close to that stratum's
other batches, and the confidence intervals overlap substantially (dMMR: +22.4 vs +18.2;
pMMR: −2.55 vs −3.42). The elevation of batch05's *overall* mean is carried by its dMMR share.

**Explanation — not established.** These data **cannot** distinguish "batch05 contains more dMMR
tumours" from "something about batch05 processing affects dMMR arrays specifically", because the
two are perfectly confounded in this deposit. No causal claim is made.

---

## 7. Sex-axis discordance — reassessment of the eight

Marker panel widened from Phase 8A's 6 Y-linked genes to **11** (adding ZFY, NLGN4Y, TXLNGY,
TMSB4Y, PRKY), with per-gene male–female separation 0.16–4.86 log₂.

| sample | recorded sex | batch | Y (6-gene) | XIST | Y markers consistent | 7C soft flag | ACTB 3′/5′ | mean r |
|---|---|---|---|---|---|---|---|---|
| GSM972036 | Female | batch03 | 6.28 | 4.78 | 10/11 | — | 1.436 | 0.946 |
| GSM972111 | Male | batch03 | 3.64 | 6.85 | 10/11 | — | 0.668 | 0.942 |
| GSM972113 | Female | batch03 | 6.53 | 3.27 | 10/11 | **flagged** | 0.969 | 0.927 |
| GSM972174 | Male | batch04 | 3.60 | 8.21 | **11/11** | — | 0.989 | 0.943 |
| GSM972248 | Male | batch05 | 4.32 | 8.44 | **11/11** | — | 0.574 | 0.959 |
| GSM972283 | Male | batch05 | 3.26 | 8.84 | **11/11** | **flagged** | 1.111 | 0.930 |
| GSM972479 | Male | batch07 | 3.31 | 9.56 | **11/11** | **flagged** | 0.587 | 0.897 |
| GSM972495 | Female | batch08 | 6.20 | 2.54 | 10/11 | — | 1.082 | 0.933 |

**Observations.** Discordance is consistent across the *whole* widened panel (10–11 of 11
markers), not an artefact of one probe. The eight span **five** batches, arguing against a single
plate-level handling error. Degradation ratios (ACTB 0.574–1.436) sit inside the cohort range
(0.152–2.629).

**Statistical association.** 3/8 carry a Phase 7C soft flag against 29/453 cohort-wide —
Fisher **OR = 9.67, p = 0.011**.

**Explanation — unresolved.** A Phase 7C soft flag reports array **quality**, not array
**identity**. The enrichment is consistent with degraded arrays producing unreliable marker calls,
and equally with coincidence in eight samples; it is evidence for neither a mix-up nor its
absence. Five of the eight carry no technical flag at all.

> **No independent technical evidence establishes a sample swap or a sex-label error.** Expression
> markers alone cannot separate a mix-up from loss of Y or XIST silencing; genotype-level data
> would be required. **These eight remain an unresolved QC observation. No sample is excluded or
> relabelled. The Phase 7C decision to retain all 453 stands.**

---

## 8. Additional technical structure — new in Phase 8B

### 8.1 RNA extraction method — associated, but confounded with batch

Six protocols in S3 (Rneasy Mini 143, TriReagent+Rneasy Mini 118, Trizol 82, Rneasy Micro 81,
RNA NOW 15, Cesium Chloride 14). **Phase 8A never examined this variable.**

PC3 ε² = **0.052** (q = 3.0e-04), PC1 ε² = 0.032 (q = 1.0e-02) — comparable to batch's own PC3
ε² of 0.061. But extraction method is **almost completely confounded with batch**
(χ² = 2265.0, df = 40, p ≈ 0), so the two cannot be separated in this deposit. Under MMR
stratification the association does **not** survive (PC3 p = 0.427, PC1 p = 0.725), the same
pattern as batch.

### 8.2 Two signals that **do** survive MMR stratification

Applying §5.3's discipline to the new findings:

| Test | observed H | null p95 | p (stratified) |
|---|---|---|---|
| RNA extraction method ~ PC3 | 28.35 | 41.06 | 0.427 |
| RNA extraction method ~ PC1 | 19.26 | 36.70 | 0.725 |
| **AFFX control PC1 ~ batch** | **24.85** | 15.50 | **0.0016** |
| **expression PC2 dispersion ~ batch** | **41.94** | 15.41 | **0.0001** |

**This is the principal new finding of Phase 8B.** Phase 8A's conclusion — that the batch–PC
association is composition — holds for **location on PC1 and PC3**, but does **not** generalise:

- **AFFX control probes.** The 62 Affymetrix control probes have their own PC1 (52.0 % of control-probe variance) which is associated with batch (ε² = 0.038, p = 0.0017) **and survives MMR stratification** (p = 0.0016). Max \|ρ\| with expression PC1–5 is 0.440.
- **Dispersion, not location.** Brown–Forsythe tests across batches: **PC2 W = 4.95, p = 1e-05**; PC5 W = 2.43, p = 0.014; PC1/PC3/PC4 null. Batch-associated **variance** heterogeneity on PC2 survives MMR stratification (p = 0.0001). Phase 8A compared batch **means** only and could not have seen this.

**Interpretation — plausible, not established.** Location-and-scale correction removes mean
differences more completely than higher-moment structure, and control probes were subject to the
same deposited processing as the gene index, so they are **not** an uncorrected technical
reference. What can be said: **residual batch-associated structure exists in this matrix that
MMR composition does not explain.** Its origin is not determined here, and no claim is made that
ComBat succeeded or failed.

### 8.3 Processing-order proxy — no signal

Numeric CIT accession index vs PC1–20: strongest \|ρ\| = 0.093 (PC8, q = 0.148) — nothing
survives correction. Caveat: accession order is a **proxy** for processing order, not a record of
it; no scan or hybridisation date is deposited.

### 8.4 Dependent arrays — none

The deposit flags related arrays via `dependancy_sample`. **0 of the 453 S3 samples carry one**,
so a future resampling scheme's independence assumption is not violated by declared sample
dependency.

---

## 9. Reproduced vs new

| Finding | Status |
|---|---|
| PCA spectrum, PC1–20 scores, sign convention | **Reproduced** (independent algorithm, 5e-9) |
| `combat_mod` ↔ MMR, 0 discordant at 410 and 519 | **Reproduced** |
| batch05 is 63.3 % dMMR vs 14.3 % cohort | **Reproduced** |
| PC1/PC3 batch association absent in the pMMR stratum | **Reproduced** |
| Eight both-axes sex-discordant samples, five batches | **Reproduced** |
| Horn / broken-stick dimensionality (48 / 97 PCs) | **New** |
| PC associations extended to PC1–20 with BH q-values | **New** |
| RNA extraction method, CIT subtype, accession order examined | **New** |
| Stratified permutation test retaining all 453 | **New** |
| Detection-sensitivity curve (80 % power at ≈ 0.3 SD) | **New** |
| Nominal PC3 × batch association **within the dMMR stratum** | **New** (qualifies Phase 8A) |
| batch05 decomposition with bootstrap CIs | **New** |
| **AFFX control-probe batch structure surviving stratification** | **New** |
| **PC2 dispersion heterogeneity across batches** | **New** |
| Soft-flag enrichment among the eight (OR 9.67, p = 0.011) | **New** |
| Column-level outcome blindness on the metadata file | **New** |

---

## 10. Limitations and unresolved questions

1. **ComBat was already applied.** Between-batch variance was removed before deposit. Every null batch result here is a statement about an already-corrected matrix. **No claim is made that ComBat removed all batch effects** — §8.2 shows it did not remove everything.
2. **Detection bound.** The stratified test reaches 80 % power at ≈ 0.3 SD of PC3. Smaller effects are invisible to it.
3. **dMMR stratum is small and unbalanced** — 64 samples, 6 usable batches, min cell 2, absent from two batches entirely.
4. **Batch, extraction method and (through `combat_mod`) MMR are mutually confounded.** No analysis of this deposit can separate them.
5. **`combat_mod` is not proven to be MMR.** Perfect concordance across 519 arrays is consistent with the depositors having used `mmr.status`, and equally with any variable partitioning these arrays identically. **Unresolved — author correspondence is the only route** (Phase 7A Q1).
6. **The MLH1/EPM2AIP1 observation is hypothesis-generating only.** The 5 MMR-unrecorded `modcan1` samples sit at mean z = −0.754, between recorded dMMR (−1.343) and pMMR (+0.249). **This was not used to impute MMR status**; A-19 complete-case is frozen.
7. **The eight sex-discordant samples are unresolved** and require genotype-level evidence.
8. **AFFX controls are not an independent technical reference** — they went through the same deposited processing.
9. **`cit.molecularsubtype` is circular** and is not evidence about structure.
10. **U1 / U2 unchanged** — `rfs.delay` unit `[PROVISIONAL: months]`, death-without-recurrence handling. Untouchable in an outcome-blind phase.

---

## 11. Decision log

**No decision was changed by Phase 8B.** This phase was investigative only.

| Decision | Phase | Status after 8B |
|---|---|---|
| Cohort S3 = 453, all retained | 7C D2 | **UNCHANGED** — 0 excluded, 0 relabelled |
| H-C1 subset criterion | 7C D1 (rev. 3) | **UNCHANGED** |
| batch05 / GSM972311 are observations, not exclusion grounds | 7C D3 | **UNCHANGED** |
| Median gene-level representation | 2B | **UNCHANGED** (hash verified) |
| Analysis bundle | 3 | **UNCHANGED** (9 hashes verified) |
| M2 primary = 4-set, sensitivity = 5-set | 6C | **UNCHANGED** — membership and roles untouched |
| M2 σ̂ `ddof=1` | M2-A1 | **UNCHANGED** — pin verified in force, read from spec |
| Complete-case missing data | A-19 | **UNCHANGED** |
| Endpoint definition | 6A | **NOT ACCESSED** |

The eight sex-discordant samples **remain an unresolved QC observation**. No new exclusion
criterion was created; no exploratory pattern was converted into one.

---

## 12. Decisions requested

None are required to proceed. Three items are offered for your consideration:

1. **The AFFX / PC2-dispersion residual structure (§8.2)** is new evidence that batch-associated structure remains after MMR composition is accounted for. Whether this warrants a sensitivity analysis in the eventual modelling protocol (e.g. batch as a stratification factor) is a protocol decision, not an 8B finding.
2. **`combat_mod` ≡ MMR (§4)** can only be closed by author correspondence.
3. **The eight sex-discordant samples (§7)** remain unresolved; closing them requires genotype-level data that this project does not hold.

**Phase 8B is complete. No Phase 8C or Phase 9 activity has been started.**
