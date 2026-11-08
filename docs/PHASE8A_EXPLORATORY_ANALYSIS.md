# Phase 8A — Outcome-Blind Exploratory Data Analysis

**Status:** executed · **Date:** 2026-11-02 · **Script:** `scripts/32_phase8a_eda.py`
**Outputs:** `reports/phase8a/` (21 CSV + 1 JSON + 6 figures)
**Cohort:** S3, n = 453, all retained (Phase 7C decision D2)
**Representation:** frozen median gene-level aggregation (Phase 2B), 21,815 genes
**Matrix:** deposited GSE39582 values used **as-is** — no normalization, no ComBat, no batch correction

---

## 0. Compliance

| Restriction | Status | How it is enforced |
|---|---|---|
| No survival modelling | ✅ none | — |
| No recurrence association testing | ✅ none | — |
| No outcome column opened | ✅ **blocked in code** | `_guard_outcome_files()` patches `pandas.read_csv/read_parquet/read_table` to raise `PermissionError` on any path containing `endpoint.csv`, for the whole run |
| No feature selection | ✅ none | gene statistics written for **all 21,815** genes; no list is filtered, cut, or carried forward |
| No sample exclusion | ✅ 453/453 retained | asserted (`len(s3) == 453`) |
| No normalization / ComBat | ✅ none | — |
| No package installation | ✅ none | numpy/pandas/scipy/matplotlib, all pre-existing |
| Frozen artifacts unmodified | ✅ 14/14 hashes identical | hashed before **and** after the run; `assert not changed` |
| No commits | ✅ none | — |

Three analyses that are conventionally part of EDA were **deliberately not run**, because each
would have required outcome access or constituted selection: variance-based gene filtering,
supervised dimension reduction, and any stratification of figures by recurrence status.

---

## 1. Expression-level analysis

### 1.1 Per-sample distributions

| Statistic | Value |
|---|---|
| Matrix | 21,815 genes × 453 samples, **0 non-finite values**, **0 zero-variance genes** |
| Global range | 1.465 – 16.112 (log₂), mean 5.383, SD 2.007 |
| Per-sample median | 4.974 – 5.276 (total spread **0.302**; SD across samples **0.034**) |
| Per-sample IQR | 2.493 – 3.407 (SD across samples 0.132) |
| Per-sample mean | 5.383 ± 0.030 |

Per-sample central tendency is tightly clustered: the full range of sample medians, 0.302 log₂
units, is about one tenth of the median *gene's* own dynamic range. This is the expected
signature of quantile-normalised RMA output and is consistent with the Phase 7C hard-criteria
result (0 hard failures). Dispersion varies more than location — the IQR range (0.91 log₂) is
three times the median range.

### 1.2 Gene-level variability

Gene SD deciles (log₂ units): **0.079** (min), 0.181, 0.220, 0.255, 0.299, **0.351** (median),
0.408, 0.471, 0.561, 0.712 (p90), **3.110** (max).

| | |
|---|---|
| Genes with SD > 1.0 | 733 (3.4%) |
| Genes with SD > 0.5 | 5,784 (26.5%) |
| Genes with SD < 0.2 | 3,175 (14.6%) |
| Variance in the top 1% of genes | **13.5%** of total |
| Variance in the top 10% of genes | **45.5%** of total |
| Spearman(gene mean, gene SD) | **+0.633** |

Variance is concentrated but not pathologically so. The strong positive mean–SD rank correlation
is the residual intensity-dependent heteroscedasticity that log₂-RMA reduces but does not
eliminate; low-expressed genes are compressed toward the array background.

The gene-mean distribution (deciles 2.29 … 8.08, max 13.99) is right-skewed with a heavy low
shoulder — the usual "mostly not expressed / some expressed" structure of a whole-genome array
on one tissue.

**The 15 most variable genes** (descriptive only — no gene is selected or carried forward):
REG4, PRAC1, REG1A, OLFM4, CLCA1, REG1B, **RPS4Y1**, **UGT2B17**, ITLN1, DEFA5, SPINK4, LEFTY1,
MS4A12, CLCA4, **XIST**.

Twelve of the fifteen are intestinal differentiation / secretory-lineage genes (regenerating-islet
family, goblet and Paneth-cell markers, normal-colonocyte markers MS4A12 and CLCA4). The other
three are not tumour biology at all: RPS4Y1 and XIST are sex-linked, and UGT2B17 carries a common
germline deletion polymorphism. This is a useful calibration — **the highest-variance axis of this
matrix is partly germline and partly normal-epithelium content, not disease process**, which is
one argument for the pre-specified pathway model (M2) over a variance filter.

### 1.3 Sample-to-sample correlation

| | min | p1 | Q1 | median | Q3 | max |
|---|---|---|---|---|---|---|
| Pearson | 0.798 | 0.883 | 0.933 | **0.945** | 0.954 | 0.986 |
| Spearman | 0.818 | 0.890 | 0.937 | **0.947** | 0.955 | 0.985 |

No pair falls below 0.80. The distribution is unimodal with a short left tail.

### 1.4 Distributional heterogeneity

Statistic: sup-norm between each sample's quantile function (percentiles 1–99) and the
cohort-median quantile function, in log₂ units.

Range **0.023 – 0.516**, median **0.108**. The six largest: GSM972105 (0.516, batch03),
GSM972097 (0.449, batch03), GSM971986 (0.436, batch01), GSM972017 (0.430, batch02),
GSM972264 (0.362, batch05), GSM972267 (0.341, batch05).

Even the most deviant sample's quantile curve sits within ~0.52 log₂ of the cohort median across
the whole range — well under the 2.0 log₂ spread of a typical gene.

### 1.5 Positive control — sex-linked expression vs recorded sex

An outcome-blind identity check. The two axes (Y-linked mean over RPS4Y1/DDX3Y/KDM5D/EIF1AY/
USP9Y/UTY, and XIST) are classified **independently** rather than by a single Y-minus-XIST
margin, because the two common tumour-biological causes of apparent discordance move only one
axis while a sample mix-up moves both. Thresholds are the midpoints of each axis's 1-D 2-means
centres (Y > 4.93, XIST > 5.95), derived from the data and reported.

| Classification | n | Rate |
|---|---|---|
| Concordant on both axes | 376 / 453 | 83.0% |
| **Discordant on both axes** | **8** | **1.8% of cohort** |
| Y axis only (LOY-like) | 31 | all 31 are male — 12.2% of 255 males |
| XIST axis only | 38 | 37 are female — 18.7% of 198 females; 1 male |

The 8 both-axes cases — GSM972036, GSM972111, GSM972113, GSM972174, GSM972248, GSM972283,
GSM972479, GSM972495 — are spread across five ComBat batches (03, 03, 03, 04, 05, 05, 07, 08),
not concentrated in one, which argues against a single plate-level handling error.

**Not adjudicated.** Array expression alone cannot separate a mix-up from loss of Y or XIST
silencing; genotype-level data would be required. **No sample is excluded or relabelled**
(Phase 7C D2 is frozen, and exclusion is prohibited in this phase). Recorded as a descriptive
observation only.

---

## 2. PCA and sample structure

### 2.1 Method

Gene-wise centring, **no scaling** (matching the Phase 7B criterion D-S1 convention), then a thin
SVD of the centred 21,815 × 453 matrix via `numpy.linalg.svd(full_matrices=False)`. Fully
deterministic — no randomised solver, no RNG. Component signs are fixed by forcing the
largest-magnitude loading positive, so the output is reproducible run to run.

> **This full-cohort rotation is EXPLORATORY.** Phase 6B §3.2.2 requires any PCA rotation used as
> a model input to be re-estimated inside each training fold. It is not a model input here and
> must not become one.

### 2.2 Variance explained

| PC | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| Gene-level (primary) | **10.32%** | 7.68% | 6.36% | 5.92% | 3.21% |
| Gene-standardised (sensitivity) | 11.58% | 5.91% | 5.26% | 3.21% | 2.45% |
| Probe-level (Phase 7C D-S1) | 9.73% | 8.09% | 6.76% | 4.95% | 3.82% |

Cumulative: PC1–5 **33.5%**, PC1–10 42.6%, PC1–20 52.5%. **17** PCs reach 50% of variance;
**127** reach 80%.

The gene-level and the independently computed probe-level spectra agree closely, as they should —
median aggregation changes the row space but not the dominant sample-space structure.

The spectrum is **shallow**: no single component dominates, and the effective dimensionality is
high. There is no evidence of one overwhelming technical gradient.

### 2.3 PC associations with annotation (descriptive)

Kruskal–Wallis with ε² effect size for categorical annotation, Spearman ρ for age. These are
*descriptive associations in an outcome-blind annotation set*, not hypothesis tests about
recurrence. With n = 453, small p-values are cheap; effect sizes are the informative column.

| Annotation | PC | ε² | p |
|---|---|---|---|
| CIMP status | PC3 | **0.207** | 3.2e-19 |
| MMR status | PC3 | 0.199 | 1.3e-19 |
| CIN status | PC3 | 0.198 | 1.6e-17 |
| `combat_mod` | PC3 | 0.195 | 3.7e-21 |
| MMR status | PC1 | 0.183 | 3.6e-18 |
| `combat_mod` | PC1 | 0.164 | 5.3e-18 |
| CIMP status | PC1 | 0.138 | 2.1e-13 |
| BRAF mutation | PC3 | 0.135 | 1.4e-13 |
| Tumour location | PC1 | 0.100 | 1.1e-11 |
| **ComBat batch** | **PC3** | **0.061** | 2.8e-05 |
| **ComBat batch** | **PC1** | **0.029** | 7.5e-03 |
| TNM stage | PC3 | 0.009 | 0.027 |
| Sex | PC4 | 0.008 | 0.036 |
| Age (Spearman ρ) | PC5 | −0.195 | 3.0e-05 |
| Deposit partition (discovery/validation) | PC1 | 0.005 | 0.070 |

The leading components track the **MMR / CIMP / CIN / BRAF / proximal–distal** axis — the
canonical molecular axis of colorectal cancer. Batch associations are an order of magnitude
weaker. The discovery/validation partition shows essentially no expression-level separation
(ε² ≤ 0.005), which is consistent with Phase 1's finding that these are not independent external
cohorts.

### 2.4 Clustering

Average linkage on 1 − Pearson sample correlation. The tree is dominated by one large cluster with
a handful of peelable singletons at every cut:

| k | Cluster sizes |
|---|---|
| 2 | 452, 1 |
| 3 | 451, 1, 1 |
| 4 | 449, 2, 1, 1 |
| 5 | 446, 3, 2, 1, 1 |
| 6 | 445, 3, 2, 1, 1, 1 |

**GSM972311 is the singleton at every k** — independently reproducing the Phase 7C OBS-2
observation at gene level (it was found at probe level there). It also has the lowest mean
correlation to other samples (0.834). GSM972357 splits off at k = 3 (mean r 0.894); the pair
{GSM972275, GSM972361} at k = 4.

There is **no two-group structure** in this cohort at the whole-transcriptome level — no clean
split that would suggest two populations or a processing dichotomy.

### 2.5 Technical vs biological structure — `combat_mod` resolved, and the batch05 observation re-read

Phase 7A left open question **Q1**: the ArrayExpress SDRF exposes a ComBat model covariate
`normalizationcombatmod` with values `modcan1`/`modcan2` whose meaning is **undocumented**.

**Direct observation.** `combat_mod` is perfectly concordant with the deposited MMR status:

| | dMMR | pMMR |
|---|---|---|
| `modcan1` | 60 | 0 |
| `modcan2` | 0 | 350 |

**Zero discordant samples** among the 410 S3 samples with MMR recorded; **zero** among all **519**
tumour arrays in the SDRF with MMR recorded.

This matters for the batch05 observation. ComBat batches differ sharply in dMMR composition
(χ² = 136.0, df = 8, p = 1.6e-25): **batch05 is 63.3% `modcan1` against a cohort rate of 14.3%**,
while batch02 and batch04 contain none at all.

Stratifying the PC–batch association by `combat_mod` removes it entirely:

| PC | marginal ε² | p | pMMR-stratum ε² | p |
|---|---|---|---|---|
| PC1 | 0.029 | 0.008 | **−0.018** | 0.996 |
| PC2 | 0.021 | 0.026 | **−0.012** | 0.896 |
| PC3 | **0.061** | **2.8e-05** | **−0.012** | 0.894 |
| PC4 | −0.006 | 0.702 | −0.010 | 0.829 |
| PC5 | −0.001 | 0.463 | −0.012 | 0.901 |

batch05's PC3 mean of +13.2 decomposes into +22.4 for its 31 `modcan1` samples and **−2.55** for
its 18 `modcan2` samples — the latter indistinguishable from every other batch.

---

## 3. Frozen pathway-score analysis (M2)

Scores computed with the frozen rule, verbatim:
`score_P(i) = mean_{g ∈ P ∩ index, σ_g > 0} (x_gi − μ_g^train) / σ_g^train`.

The implementation makes `train_ids` a **mandatory argument with no default**, so a complete-cohort
fit cannot happen by omission — a caller must always state the estimation population.

### 3.1 Coverage against the frozen rule (≥80% of source members **and** ≥20 surviving genes)

| Gene set | source | resolved & present | coverage | zero-variance dropped | rule |
|---|---|---|---|---|---|
| EPITHELIAL_MESENCHYMAL_TRANSITION | 200 | 198 | 99.0% | 0 | **MET** |
| E2F_TARGETS | 200 | 195 | 97.5% | 0 | **MET** |
| WNT_BETA_CATENIN_SIGNALING | 42 | 41 | 97.6% | 0 | **MET** |
| P53_PATHWAY | 200 | 194 | 97.0% | 0 | **MET** |
| PI3K_AKT_MTOR_SIGNALING *(5-set only)* | 105 | 104 | 99.0% | 0 | **MET** |

No score was dropped, in any fold. Frozen membership and coverage rules are unchanged; no pathway
was added, removed, promoted or demoted.

### 3.2 Full-cohort scores — EXPLORATORY, NOT A MODEL INPUT

Files are named `pathway_scores_*_EXPLORATORY_full_cohort.csv` and carry
`full_cohort_scores_used_as_model_input: false` in the summary JSON.

| Set (4-set primary) | mean | SD | range |
|---|---|---|---|
| EMT | 0.000 | 0.495 | −1.105 … +2.156 |
| E2F_TARGETS | 0.000 | 0.554 | −2.479 … +1.306 |
| WNT_BETA_CATENIN | 0.000 | 0.253 | −0.653 … +0.854 |
| P53_PATHWAY | 0.000 | 0.195 | −0.555 … +0.761 |

Score SD is well below 1 because averaging ~200 imperfectly correlated standardised genes shrinks
the dispersion; this is a property of the frozen averaging rule, not a defect.

### 3.3 Collinearity — a check against the frozen record

| | 4-set (S3, n=453) | 4-set (Phase 6C, n=461) | 5-set (S3) | 5-set (Phase 6C) |
|---|---|---|---|---|
| max \|r\| | 0.316 | 0.318 | **0.631** | — |
| median \|r\| | 0.186 | 0.181 | 0.186 | — |
| n(\|r\| > 0.5) | **0** | 0 | **1** | 1 |
| condition number | 2.949 | 2.950 | 7.032 | 6.95 |
| min eigenvalue | 0.481 | 0.480 | 0.278 | 0.280 |
| PCs for 80% | 3 | 3 | 3 | — |

The Phase 6C figures are reproduced on the S3 cohort. The single |r| > 0.5 pair in the 5-set is
**P53_PATHWAY ↔ PI3K_AKT_MTOR (r = 0.631)** — identifying the specific pair behind the conditioning
penalty that made PI3K_AKT_MTOR a sensitivity-only addition. The immutability clause stands: it is
**not** promoted or removed on this basis, and this observation is outcome-blind, so it could not
serve as grounds either way.

### 3.4 Fold isolation — the leakage invariant, demonstrated

Five **outcome-blind deterministic splits** (seed 20261102, permutation of sorted sample IDs).

> These splits are **illustrative machinery verification only**. They are **not** the Phase 6B
> resampling scheme, which is defined at modelling time and may use outcome information for
> stratification — something this phase cannot access.

| Check | Result |
|---|---|
| μ/σ bit-identical when held-out samples are **deleted from the input matrix** | **5/5 folds** (asserted) |
| Coverage rule met in every fold | 5/5, both panels |
| Zero-variance exclusions | 0 in every fold |
| Max \|score shift\|, held-out samples, train-fold fit vs full-cohort fit | **0.0722 SD units** |

The last row is the quantity that matters: it is the magnitude of the leakage that a full-cohort
scoring step would have introduced. It is small but **not zero**, which is exactly why the frozen
rule puts standardisation inside the loop.

### 3.5 σ degrees of freedom — found unpinned here, **since pinned by amendment M2-A1**

> **Resolved 2026-11-02.** This section records the gap as Phase 8A found it. The approved
> amendment **M2-A1** pins σ̂ to **`ddof = 1`** in both frozen specifications and in the freeze
> record (spec revision 1 → 2); see `docs/PHASE6C_PATHWAY_PANEL_REVIEW.md` §0. Phase 8A was
> re-run afterwards and every pathway-score artifact is **byte-identical**, because this phase
> had already used `ddof = 1`. The scoring implementation now reads `scoring.sigma_ddof` from
> the specification and **raises** if it is absent, rather than inheriting a library default.

As found: the M2 freeze record specified "per-gene μ̂/σ̂ from training fold only" but did **not
pin the denominator** of σ̂. numpy defaults to ddof = 0, pandas to ddof = 1.

Measured effect: every score is rescaled by exactly √(n/(n−1)) — **1.0011055835** at n = 453,
matching the analytic constant to 10 decimal places. Because every gene in a fold shares the same
n_train, this is a single constant rescaling per fold: **rank order and all correlations are
unchanged, but the numeric score values are not.**

This phase used **ddof = 1** — now the pinned convention. The gap was real, and it was closed
before modelling rather than after seeing results, which is the only point at which such a
choice can be made without it becoming a tuning decision.

### 3.6 Pathway scores vs principal components (Spearman, descriptive)

| Pathway | PC | ρ |
|---|---|---|
| EMT | PC1 | **+0.828** |
| E2F_TARGETS | PC4 | −0.552 |
| E2F_TARGETS | PC2 | +0.511 |
| P53_PATHWAY | PC3 | +0.435 |
| WNT_BETA_CATENIN | PC4 | −0.418 |
| P53_PATHWAY | PC1 | +0.395 |
| WNT_BETA_CATENIN | PC2 | −0.352 |

Mean scores by annotation (descriptive):

| | EMT | E2F | WNT | P53 |
|---|---|---|---|---|
| `modcan1` (dMMR) | +0.188 | +0.244 | −0.089 | +0.201 |
| `modcan2` (pMMR) | −0.032 | −0.041 | +0.015 | −0.034 |
| proximal | +0.092 | +0.008 | −0.022 | +0.061 |
| distal | −0.063 | −0.006 | +0.015 | −0.042 |
| CIMP+ | +0.094 | +0.167 | −0.106 | +0.130 |
| CIMP− | −0.011 | −0.022 | +0.047 | −0.011 |

EMT and E2F_TARGETS are **anti-correlated** (r = −0.271).

---

## 4. Scientific interpretation

The four tiers below are kept strictly separate. Nothing in §4.1 depends on interpretation;
nothing in §4.2 or §4.3 is established.

### 4.1 Direct observations

1. The matrix is complete and well-conditioned: 0 non-finite values, 0 zero-variance genes, no sample pair below r = 0.80.
2. Location is tightly controlled across samples (median spread 0.302 log₂); dispersion varies three times more (IQR spread 0.91).
3. Variance is concentrated: 1% of genes carry 13.5% of total variance, 10% carry 45.5%.
4. The PCA spectrum is shallow — PC1 10.3%, and 127 components are needed for 80% of variance.
5. The leading PCs track MMR / CIMP / CIN / BRAF / tumour location (ε² 0.10–0.21). Stage, sex and the discovery/validation partition show ε² ≤ 0.009.
6. `combat_mod` (`modcan1`/`modcan2`) is **perfectly concordant with deposited MMR status**: 0 discordant among 410 S3 samples and 0 among 519 SDRF tumour arrays.
7. ComBat batches differ sharply in dMMR composition (χ² = 136.0, p = 1.6e-25); batch05 is 63.3% dMMR vs 14.3% cohort-wide.
8. The PC–batch association does not survive stratification by `combat_mod` (PC3 ε²: 0.061 → −0.012, p: 2.8e-05 → 0.894).
9. GSM972311 is the singleton at every clustering cut and has the lowest mean correlation (0.834). There is no two-group structure.
10. All five frozen gene sets meet the frozen coverage rule (97.0–99.0%), with zero zero-variance exclusions.
11. Phase 6C's collinearity figures reproduce on S3; the single |r| > 0.5 pair in the 5-set is P53_PATHWAY ↔ PI3K_AKT_MTOR (r = 0.631).
12. Training-fold parameters are bit-identical when held-out samples are deleted from the input (5/5 folds); full-cohort scoring would shift held-out scores by up to 0.0722 SD.
13. 8/453 samples (1.8%) are discordant with recorded sex on **both** sex-linked axes, spread across five batches; 31 males are LOY-like and 37 females are XIST-low.

### 4.2 Biological hypotheses *(generated here, not tested here)*

- **H1.** The dominant transcriptional axis of this cohort is the MSI/CIMP/BRAF/proximal axis, consistent with the CMS1 immune subtype of Guinney et al. (*Nat Med* 2015;21(11):1350–1356). Observations 5 and 7 are consistent with it; nothing here tests it.
- **H2.** The high EMT score in dMMR/CIMP+ tumours (§3.6) may reflect **stromal and immune cell content** rather than epithelial trans-differentiation — the HALLMARK EMT set is rich in ECM and stromal genes, and dMMR tumours are characteristically immune-infiltrated. Bulk expression cannot separate the two; a deconvolution step would be required.
- **H3.** The EMT/E2F anti-correlation (r = −0.271) is consistent with the proliferative-versus-mesenchymal trade-off reported across carcinomas. Descriptive only.
- **H4.** The LOY-like (12.2% of males) and XIST-low (18.7% of females) rates are within published ranges for solid tumours, so most of those 68 samples are plausibly biological rather than technical.
- **H5.** The high-variance gene list being dominated by secretory/differentiation markers suggests a substantial **normal-epithelium and mucinous-content** gradient across samples.

### 4.3 Technical explanations

- **T1.** The tight location control and the complete absence of non-finite values are the expected signature of deposited quantile-normalised RMA output. They are evidence about *processing*, not about biological validity.
- **T2.** The mean–SD rank correlation of +0.633 is residual intensity-dependent heteroscedasticity, standard for log₂ RMA.
- **T3. The batch05 observation (Phase 7C OBS-1) is explained by composition, not by residual technical artefact** — observation 8. This *strengthens* the Phase 7C D3 decision with a positive mechanism, where D3 previously rested on "no independent evidence found". **It does not prove that no batch effect exists:** ComBat was already applied to this deposit, so batch variance was removed by construction and this analysis cannot recover what was removed. The conditional test also stratifies on one variable only.
- **T4.** If the ComBat `mod` matrix did encode MMR status, then MMR-associated variance was explicitly **protected** during correction. That would make the strong MMR/PC association in observation 5 partly a *consequence of the preprocessing choice*, not purely a biological finding — the two cannot be separated from the deposited matrix alone.
- **T5.** The deposited matrix was ComBat-corrected with a two-level `mod` covariate that does **not** encode sex, so any confounding between sex and ComBat batch was unprotected during correction. Relevant to §1.5.

### 4.4 Unresolved questions

| # | Question | Status |
|---|---|---|
| **U4** | Does `combat_mod` *being* MMR status follow from perfect concordance? **No.** Concordance is consistent with the depositors having used `mmr.status` in the mod matrix, and equally with any variable partitioning these 519 arrays identically. Resolving it needs author correspondence — the same route Phase 7A identified for Q1. **The 43 S3 samples with MMR unrecorded are NOT imputed from `combat_mod`; A-19 complete-case is frozen.** |
| ~~**U5**~~ | ~~σ̂ degrees of freedom are not pinned by the M2 freeze record (§3.5).~~ **RESOLVED 2026-11-02 by amendment M2-A1** — pinned to `ddof = 1` in both specifications and the freeze record (spec revision 2). Membership, roles and coverage rules unchanged; Phase 8A re-run produced byte-identical scores. |
| **U6** | The 8 both-axes sex-discordant samples cannot be adjudicated from array data. No action taken. |
| **U7** | Whether GSM972311's isolation is technical or biological remains open — unchanged from Phase 7C, where 29 investigations found no independent evidence of technical invalidity. |
| **U1, U2** | Unchanged from Phase 1: `rfs.delay` time unit `[PROVISIONAL: months]`; death-without-recurrence handling. Neither is touchable in an outcome-blind phase. |

### 4.5 What this phase does **not** claim

No biomarker is identified. No pathway is associated with recurrence. No predictive value of any
kind is established or suggested. Every association reported here is between expression and
**outcome-blind annotation** (batch, MMR, CIMP, CIN, BRAF, stage, sex, age, location, partition).
The recurrence endpoint was never opened — enforced in code, not by convention.

---

## 5. Deliverables

**Script:** `scripts/32_phase8a_eda.py` (958 lines, deterministic, seed 20261102)

**Figures** (`reports/phase8a/figures/`):

| File | Content |
|---|---|
| `F1_expression_distributions.png` | per-sample quantile bands; median vs IQR by batch; gene-mean histogram |
| `F2_gene_variability.png` | gene-SD histogram; mean–variance hexbin; variance concentration curve |
| `F3_sample_correlation.png` | pairwise-r histogram; per-sample mean r; clustered correlation matrix |
| `F4_pca.png` | scree + cumulative; PC1–PC2 by ComBat batch, by deposit partition, by CIMP |
| `F5_pathway_scores.png` | 4-set score boxplots; inter-score correlation matrix; EMT vs E2F by batch |
| `F6_pathway_vs_pca.png` | PC1–PC2 coloured by each of the four primary pathway scores |

**Machine-readable summary:** `reports/phase8a/phase8a_summary.json` — every statistic in this
report, plus the 14 frozen-artifact hashes and the full environment record.

**Tables** (`reports/phase8a/`): `sample_distribution_summary.csv` · `gene_variability.csv`
(all 21,815) · `sample_correlation_pearson.csv` · `sample_mean_correlation.csv` ·
`distributional_heterogeneity.csv` · `sex_positive_control.csv` · `pca_scores_gene_level.csv` ·
`pca_variance_explained.csv` · `pc_annotation_associations.csv` · `pc_by_batch.csv` ·
`cluster_by_batch.csv` · `combat_mod_vs_mmr.csv` · `batch_by_combat_mod.csv` ·
`pc_batch_conditional_on_combat_mod.csv` · `pathway_scores_{4set_primary,5set_sensitivity}_EXPLORATORY_full_cohort.csv` ·
`pathway_score_correlation_{4set_primary,5set_sensitivity}.csv` · `fold_isolation_{4set_primary,5set_sensitivity}.csv` ·
`pathway_vs_pc_and_batch.csv`

**Tests:** `tests/test_phase8a_eda.py`

---

## 6. Decisions requested before Phase 8B / modelling

1. ~~**U5 — pin σ̂ ddof**~~ — **DONE.** Approved and implemented as amendment **M2-A1** (spec revision 2); see `docs/PHASE6C_PATHWAY_PANEL_REVIEW.md` §0.
2. Whether the `combat_mod` ≡ MMR observation should be pursued via author correspondence, which is the only route that can close Phase 7A Q1.
3. Whether the 8 both-axes sex-discordant samples warrant any recorded action beyond this observation. **Exclusion is not proposed**; D2 is frozen.

No further EDA or modelling will proceed without approval.
