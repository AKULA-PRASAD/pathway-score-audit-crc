# Phase 8D — Consolidated Exploratory Analysis and Decision Review

**Version:** 1.0 · **Status:** executed · **Date:** 2026-11-03
**Script:** `scripts/35_phase8d_consolidation.py` · **Tests:** `tests/test_phase8d_consolidation.py`
**Artifacts:** `reports/phase8d/` — `evidence_table.csv` · `consolidated_claims.csv` ·
`phase9_readiness_checklist.csv` · `phase8d_summary.json` · `decision_register.json`

> **Consolidation only.** No new analysis was performed. Every number below is **resolved from a
> source artifact by path and JSON pointer** at run time, not transcribed: `consolidated_claims.csv`
> records all 48 with their origin. No outcome was accessed. **Phase 8D decides nothing** — the
> eight register entries all read `AWAITING_APPROVAL`.

---

## Executive summary

Three exploratory phases over the deposited GSE39582 S3 cohort (n = 453) are mutually consistent:
**15 of 15 cross-phase corroboration checks agree**, including an independent re-derivation of the
Phase 8A PCA (signed correlation 1.0 to 10⁻¹²) and identical pathway coverage between 8A and 8C.

The data are clean and complete: 0 non-finite values, 0 zero-variance genes, no sample pair below
r = 0.80. The dominant transcriptional axis is the **MMR/CIMP/CIN/BRAF/tumour-location** molecular
axis, not stage and not the discovery/validation partition. Batch **location** effects on both
principal components and pathway scores are fully explained by MMR composition — **0 of 10** score
tests survive conditioning. Two residual signals do survive and are new: **AFFX control-probe
structure** and **PC2 dispersion** across batches. Neither is a location effect.

The frozen M2 scoring machinery audits clean: the production scorer agrees with an independent
re-implementation to 3 × 10⁻¹⁵, fold isolation holds under deletion *and* poisoning, and Phase 6C's
collinearity figures reproduce across three phases.

**Nothing found here blocks Phase 9 entirely.** Four items block a *specific* analysis and four are
documented limitations. The substantive open questions are a **score-standardisation policy**
(the five scores span a 3× range in dispersion, which interacts with penalisation), a **4-set/5-set
reading rule** fixed before unblinding, and the carried-forward endpoint and external-validation
constraints.

Three things this report does **not** say: that ComBat removed all batch effects; that `combat_mod`
*is* the depositors' MMR covariate; that any sample was swapped. **No outcome association has been
tested, so no pathway is prognostic, predictive, or a biomarker on this evidence.**

---

## A. Cohort, data and provenance

| | |
|---|---|
| **Analysis population** | S3, **n = 453**, all retained — Phase 7C decision D2, unchanged by 8A/8B/8C/8D |
| **Expression source** | GSE39582, GPL570 (HG-U133 Plus 2.0), deposited matrix used **as-is** |
| **Representation** | frozen median gene-level aggregation, **21,815 genes** (Phase 2B) |
| **Sensitivity representation** | frozen representative-IQR, `discovery_443` selection, same gene index |
| **Preprocessing (documented)** | deposited as "RMA normalisation + ComBat" |
| **Preprocessing (inferred, NOT documented)** | the ComBat **fitting population** and **model matrix**. `combat_mod` (`modcan1`/`modcan2`) is concordant with MMR but its identity is **unresolved** — Phase 7A Q1 |
| **Applied in Phases 8A–8D** | no normalization, no batch correction, no fRMA, no CEL acquisition |

**Documented fact vs inference, kept apart.** That the deposit says "RMA + ComBat" is documented.
*Which* arrays were in the fitting population and *which* covariates entered the model matrix is
not. The ComBat batch assignments themselves were recovered from the ArrayExpress SDRF (absent
from GEO) — a reconstruction, not a GEO-documented fact.

Sources: `docs/PHASE2_GENE_REPRESENTATION.md`, `docs/PHASE3_DATASET.md`,
`docs/PHASE7A_PREPROCESSING_FEASIBILITY.md` §Q1, `data/processed/frozen/frozen_representation.json`.

---

## B. Expression-level findings (Phase 8A)

| Finding | Value | Type | Source |
|---|---|---|---|
| Matrix completeness | 0 non-finite, 0 zero-variance genes | observation | `reports/phase8a/phase8a_summary.json` |
| Per-sample median range | 4.974 – 5.276 log₂ (spread 0.302) | observation | §1.1 |
| Pairwise sample correlation | min 0.798, median 0.945 | observation | §1.3 |
| Mean–SD rank correlation | **+0.633** | association | §1.2 |
| Variance concentration | top 1 % of genes → 13.5 % of variance; top 10 % → 45.5 % | observation | §1.2 |
| Sex control, concordant on both axes | 376 / 453 | observation | §1.5 |

**Descriptive.** The completeness, the tight location control and the correlation distribution are
descriptive properties of the deposited RMA output. They support using the deposit as-is for
relative-risk modelling; they say nothing about biological validity.

**Potential technical implication.** The mean–SD correlation of +0.633 is residual
intensity-dependent heteroscedasticity — low-expressed genes are compressed toward array
background. This is an argument against naive variance-based gene filtering, and it is one reason
the project uses a pre-specified pathway panel rather than a variance filter. The most variable
genes include **RPS4Y1, XIST and UGT2B17** — sex-linked and a common germline deletion
polymorphism — confirming that the highest-variance axis is partly germline and normal-epithelium
content, not disease process.

Full detail: `docs/PHASE8A_EXPLORATORY_ANALYSIS.md` §1.

---

## C. PCA and batch structure (Phase 8B)

### C.1 Dimensionality

PC1 **10.32 %**; **17 / 127 / 237** PCs for 50 / 80 / 90 % of variance. Horn parallel analysis
gives **48** components above the permuted 95th percentile; broken-stick gives 97. The criteria
disagree on how many components matter — that disagreement is preserved rather than resolved.
The spectrum is shallow: **no single dominant technical gradient exists to adjust for.**

Phase 8B re-derived the Phase 8A spectrum through a covariance eigendecomposition rather than an
SVD: variance explained agrees to 4.9 × 10⁻⁹, signed score correlations to 1.0 − 10⁻¹². This is
an algorithmic reproduction, **not** a second data source.

### C.2 Structure

Leading PCs track **MMR / CIMP / CIN / BRAF / tumour location** (ε² 0.10–0.21). Stage, sex and the
**discovery/validation partition** show ε² ≤ 0.009 — consistent with the Phase 1 finding that the
two deposit subsets are **not independent external cohorts**.

### C.3 `combat_mod` and MMR

**0 discordant** among the 410 S3 samples with MMR recorded, and **0 among 519** SDRF tumour
arrays. **This is concordance, not identity.** Any variable partitioning these arrays identically
fits equally well. The 43 S3 samples with MMR unrecorded are **not** imputed from it; A-19
complete-case remains frozen. → **D8D-5**.

### C.4 MMR-stratified results, with their limits

| Stratum | PC3 ε² | p | n | min cell |
|---|---|---|---|---|
| pMMR | −0.012 | 0.894 | 388 | 13 |
| **dMMR** | **+0.130** | **0.028** | 64 | **2** |
| Conditional permutation, all 453 retained | — | **0.381** | 453 | — |

Phase 8A tested the pMMR stratum only. Phase 8B found a **nominal** PC3 × batch association within
dMMR that is **not corroborated**: it is uncorrected for multiplicity, rests on cells of 2, and the
better-powered conditional test is null. It is reported because it exists, **not** because it is
established.

**Detection floor.** The stratified test reaches 80 % power at a batch05 PC3 shift of **≈ 0.3 SD**
and has zero power at 0.2 SD. Combined with the fact that ComBat was applied before deposit —
removing between-batch variance *by construction* — **no claim is made that ComBat removed all
batch effects**, and the null results are bounded statements, not evidence of absence.

### C.5 Residual structure that survives stratification — new

| Test | p (stratified) |
|---|---|
| RNA extraction method ~ PC3 / PC1 | 0.427 / 0.725 |
| **AFFX control PC1 ~ batch** | **0.0016** |
| **PC2 dispersion ~ batch** | **0.0001** |

Batch-associated **variance** heterogeneity on PC2 (Brown–Forsythe p = 1 × 10⁻⁵) and structure in
the 62 Affymetrix control probes both survive MMR conditioning. Phase 8A compared batch **means**
only and could not have seen either. **Caveat:** the AFFX probes went through the same deposited
RMA + ComBat, so they are **not** an uncorrected technical reference. → **D8D-4**.

### C.6 batch05

batch05 is **63.3 % dMMR** against **14.3 %** cohort-wide. Within each MMR stratum its mean PC3 is
close to the other batches and the bootstrap CIs overlap (dMMR +22.4 [18.2, 26.4] vs +18.2
[13.3, 22.8]). **Observation:** the elevation of batch05's overall mean is carried by its dMMR
share. **Explanation — not established:** composition and batch-specific processing are *perfectly
confounded* in this deposit; these data cannot separate them.

### C.7 The eight sex-axis-discordant samples

Eight of 453 are discordant on **both** the Y-linked and XIST axes, consistently across an
11-marker panel, spread over **five** batches. 3/8 carry a Phase 7C soft flag against 29/453
cohort-wide (Fisher OR 9.67, p = 0.011) — but a soft flag reports array **quality**, not
**identity**, and is derived from the same matrix.

**Genotype-level evidence is unavailable to this project.** Expression markers cannot separate a
mix-up from loss of Y or XIST silencing. **No sample swap is inferred. All 453 are retained.**
→ **D8D-6**.

Full detail: `docs/PHASE8B_PCA_BATCH_INVESTIGATION.md` §§2–8.

---

## D. Pathway-score findings (Phase 8C)

### D.1 Frozen membership and method

4-set **primary**: EMT (198), E2F_TARGETS (195), WNT_BETA_CATENIN (41), P53_PATHWAY (194).
5-set **sensitivity**: the above **+ PI3K_AKT_MTOR_SIGNALING (104)**. Coverage 97.0–99.0 %, **0
zero-variance exclusions**, all rules MET — identical in 8A and 8C. Scoring: mean of
training-fold-standardised members, `sigma_ddof = 1` read from the specification (M2-A1), canonical
implementation at `src/pathway_scoring.py`.

Audit: production scorer vs an independent re-implementation **3.0 × 10⁻¹⁵**; vs the Phase 8A
artifacts **5.0 × 10⁻⁷** (CSV rounding); fold isolation **5/5 bit-identical** under deletion *and*
+1000 poisoning of held-out samples.

### D.2 Distributions and within-set coherence

| Pathway | k | SD | SD if independent | implied r̄ |
|---|---:|---|---|---|
| E2F_TARGETS | 195 | 0.554 | 0.072 | **0.304** |
| EMT | 198 | 0.495 | 0.071 | **0.241** |
| WNT_BETA_CATENIN | 41 | 0.253 | 0.156 | **0.041** |
| P53_PATHWAY | 194 | 0.195 | 0.072 | **0.033** |
| PI3K_AKT_MTOR | 104 | 0.185 | 0.098 | **0.025** |

Two regimes: two internally coherent sets, three nearly incoherent. **Low dispersion is not low
biological relevance.** Score SD is a joint property of the gene set's coherence and the averaging
transform — a mean of *k* standardised genes has SD 1/√k if members are independent. These numbers
say how *coordinated* each set is in this tissue under this preprocessing. They say nothing about
recurrence, which has not been examined. → the dispersion range does interact with penalisation,
**D8D-1**.

No sample moves any score's SD by more than **2.2 %**.

### D.3 Correlation structure

4-set: max |r| 0.316, condition 2.949, **0** pairs > 0.5. 5-set: max |r| **0.631**, condition
7.032, **1** pair — **P53 ↔ PI3K**. Phase 6C's figures (recorded at n = 461) reproduce at n = 453.

### D.4 PI3K

**R² = 0.540** on the four primary scores — 46 % of its variance not linearly captured. This is a
**descriptive redundancy measure among scores only**, with no outcome involved. It is **not**
grounds to promote or remove PI3K; the M2-A1 immutability clause stands and PI3K remains
**sensitivity-only**. → **D8D-2**.

### D.5 Representation stability (sensitivity-only)

E2F 0.998 · EMT 0.992 · P53 0.972 · PI3K 0.910 · **WNT 0.858**. The smallest and least coherent
primary set is the one whose value depends materially on the probe→gene rule. These are differences
in summarisation, **not biological differences**. → **D8D-3**.

### D.6 Fold composition

All correlations with the full-cohort fit > **0.9996**; max shift **0.072 SD** — stable but
**non-zero**, which is exactly why standardisation must stay inside the resampling loop.

### D.7 Associations, with circularity identified

**0 of 10** batch and extraction-method associations with scores survive MMR-stratified
permutation. The survivors all involve `qc_any_soft_flag` as **dispersion** effects, and that
variable is **matrix-derived**: the Phase 7B flags are Tukey rules on the median, IQR and
inter-sample correlation of this same matrix, so a flagged sample has unusual z-scores by
construction. `cit.molecularsubtype` is **circular** — the depositors' own expression-derived
clustering of this matrix — and is excluded from all non-circular counts.

Full detail: `docs/PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md` §§2–6.

---

## E. Cross-phase synthesis

Machine-readable: `reports/phase8d/evidence_table.csv` (24 findings, seven columns).

### E.1 Classification

| Class | n | Examples |
|---|---:|---|
| **Reproduced across phases** | 9 | PCA spectrum (8A↔8B, independent algorithm) · pathway coverage (8A↔8C) · 4-set collinearity (6C↔8A↔8C) · the same eight sex-discordant samples (8A↔8B) · fold-shift 0.072 SD (8A↔8C) |
| **New in 8B or 8C** | 8 | Horn/broken-stick dimensionality · AFFX batch structure · PC2 dispersion · implied within-set coherence · PI3K R² · representation stability · independent scorer agreement · dependent-array check |
| **Unresolved** | 3 | `combat_mod` identity · nominal dMMR PC3 association · the eight sex-discordant samples |
| **Not independent — same matrix or QC rules** | 4 | `cit.molecularsubtype` (circular) · `qc_any_soft_flag` (matrix-derived) · AFFX controls (same ComBat) · soft-flag enrichment among the eight |

### E.2 Selected rows

| Finding | Phase | Type | Corroboration | Limitation | Modelling implication | Status |
|---|---|---|---|---|---|---|
| Batch–PC location association vanishes under MMR stratification | 8B | association | pMMR stratum + 10,000-permutation conditional test | ComBat already removed batch variance; floor ≈ 0.3 SD | no batch location covariate indicated | documented limitation |
| Nominal PC3 × batch within dMMR | 8B | association | **not corroborated** (conditional test p = 0.38) | n = 64, min cell 2, uncorrected | too fragile to act on | unresolved |
| AFFX control PC1 ~ batch survives stratification | 8B | association | new in 8B | same ComBat — not an uncorrected reference | candidate prespecified sensitivity analysis | requires approval |
| PC2 dispersion ~ batch survives stratification | 8B | association | new in 8B; 8A compared means only | a dispersion effect is not fixed by a location covariate | → D8D-4 | requires approval |
| Two coherence regimes among the five sets | 8C | observation | new in 8C | property of sets + transform, **not relevance** | low-dispersion scores shrunk harder → D8D-1 | requires approval |
| WNT representation stability 0.858 | 8C | observation | new in 8C (sensitivity-only) | summarisation difference, not biological | → D8D-3 | requires approval |
| `combat_mod` ↔ MMR, 0/519 discordant | 8B | observation | reproduced in 8D | concordance ≠ identity | limits provenance claims only | unresolved |

---

## F. Decision register — **Phase 8D decides nothing**

Machine-readable: `reports/phase8d/decision_register.json` (`decisions_made_by_phase8d: 0`).
Every entry reads `AWAITING_APPROVAL`. **No item blocks Phase 9 entirely.**

| # | Decision | Blocking | Blocks what |
|---|---|---|---|
| **D8D-1** | Score-level standardisation | blocks a specific analysis | the penalised Cox specification only |
| **D8D-2** | 4-set vs 5-set reporting framework | documented limitation | nothing |
| **D8D-3** | WNT representation sensitivity | blocks a specific analysis | interpretation of a WNT coefficient |
| **D8D-4** | Prespecified batch sensitivity | documented limitation | nothing |
| **D8D-5** | `combat_mod` provenance | documented limitation | provenance claims only |
| **D8D-6** | Eight sex-discordant samples | documented limitation | nothing |
| **D8D-7** | Endpoint U1 / U2 | blocks a specific analysis | **all fixed-time absolute-risk analyses** |
| **D8D-8** | External validation | blocks a specific analysis | any confirmatory external claim |

**D8D-1 — Score-level standardisation.** M2 currently freezes `score_standardisation: "NONE
before the Cox model (penalisation handles scale)"`. The five scores span SD 0.554 → 0.185. A
ridge/elastic-net penalty applies the same penalty per unit coefficient, so on unstandardised
scores a low-dispersion predictor needs a *larger* coefficient for the same effect and is shrunk
harder. Leaving the scales as frozen is defensible — but it is a **choice with consequences for
relative shrinkage**, not a neutral default. Any standardisation adopted must be **training-fold
only**; a full-cohort score standardisation would reintroduce the leakage the frozen rule exists to
prevent. *No policy is chosen here.*

**D8D-2 — 4-set vs 5-set.** Adding PI3K raises the condition number to 7.03 and introduces the only
|r| > 0.5 pair. A primary/sensitivity divergence is therefore **foreseeable**, and without a
reading rule fixed in advance, choosing between them after the fact is post-hoc selection.
*Proposed framework, not a change:* report the 4-set as primary always; report the 5-set alongside,
never instead; state in advance that divergence is read as evidence about the **stability** of the
primary estimate, not as grounds to switch panels; pre-specify the divergence metric before
unblinding. **Roles are unchanged.**

**D8D-3 — WNT representation.** Observed median vs representative-IQR correlation **0.858** for
WNT (vs 0.972–0.998 for the others). A WNT result that did not hold under the sensitivity
representation would be as much a representation artefact as a finding. *No new policy established.*

**D8D-4 — Batch sensitivity.** Location effects are fully explained by composition. The two
survivors (AFFX PC1, PC2 dispersion) are **not** location effects, so a batch covariate would not
address them. Batch, extraction method and MMR are **mutually confounded** in this deposit. *No
batch covariate added; the primary model is unaltered.*

**D8D-5 — `combat_mod`.** Outstanding route: **author correspondence** with the CIT/GSE39582
depositors. *No MMR imputed; no ComBat model inferred.*

**D8D-6 — Sex-axis discordance.** **Genotype-level evidence recorded as unavailable.** All eight
retained; no exclusion, no relabelling.

**D8D-7 — Endpoint U1 / U2 (Phase 6A numbering).** U1 (`rfs.delay` unit, `[PROVISIONAL: months]`)
**blocks every fixed-time absolute-risk analysis** — t-year risk, predicted survival/CIF curves,
cloglog calibration at a horizon, time-dependent AUC/Brier, and the 60-unit endpoint — unless
explicitly labelled conditional. Relative-risk analyses are unaffected, which is why the primary
analysis can proceed. U2 (death without recorded recurrence) constrains any analysis treating that
as non-informative censoring; the deposited encoding is censoring-at-death and that censoring is
**demonstrably informative** (46/453 competing deaths, associated with age, stage and treatment).

> ⚠ **Naming hazard, recorded not resolved.** `docs/PHASE1_COHORT_ENDPOINT_INVESTIGATION.md` uses
> the labels **U1 and U2 for different items** (CIT446 stage vs M1; CIT106/CIT250 zero times). The
> **Phase 6A numbering** is the one used in this report and in 8A–8C. Renumbering a historical
> report is out of scope.

**D8D-8 — External validation.** Gate **G4 FAILED** on pre-specified criteria
(`docs/PHASE4_PREACQUISITION.md` §6.1): GSE14333 has **48** events and GSE33113 **18**, against
**≈ 338** required at C = 0.65 under an **approximation** of the Riley framework — 14 % and 5 % of
requirement. Both remain **exploratory candidates and their counts must not be combined**. GSE33113
is a **separate estimand** (stage II only, metastasis-or-recurrence composite, time in days).
GSE17536 is **excluded** for established patient overlap with GSE14333. **No suitable confirmatory
cohort has been identified**, and **this analysis is not externally validated.**

---

## G. Phase 9 readiness checklist

Machine-readable: `reports/phase8d/phase9_readiness_checklist.csv`.

| Status | Items |
|---|---|
| **READY** (9) | cohort frozen (453) · representation frozen · bundle verified · M2 panels, roles and `ddof=1` frozen · canonical scorer in `src/` · fold isolation verified · A-19 complete-case · QC executed and recorded · EDA complete and consolidated |
| **DECISION REQUIRED** (4) | D8D-1 standardisation · D8D-2 reporting framework · D8D-3 representation sensitivity · D8D-4 batch sensitivity |
| **OPEN — documented limitation** (2) | D8D-5 `combat_mod` · D8D-6 sex discordance |
| **CONSTRAINS** (2) | D8D-7 U1 → fixed-time analyses · D8D-7 U2 → censoring assumption |
| **NOT AVAILABLE** (1) | D8D-8 external confirmatory validation |
| **NOT STARTED** (1) | Phase 6B protocol finalisation (resampling scheme, Monte-Carlo stability) |

---

## H. Limitations carried into Phase 9

1. **ComBat was applied before deposit.** Every null batch result bounds what remains detectable in an already-corrected matrix. §C.5 shows it did not remove everything.
2. **Mutual confounding** of batch, RNA extraction method and MMR. No analysis of this deposit can separate them.
3. **Detection floor ≈ 0.3 SD** of PC3 for the stratified batch test.
4. **Sparse strata** — dMMR has 64 samples across 6 usable batches, min cell 2; batch and extraction method have min cell 14.
5. **Multiple exploratory comparisons** throughout 8A–8C. BH q-values are reported where tests were enumerated; single nominal p-values are weak evidence.
6. **Four findings are not independent evidence** — derived from the same expression matrix or the same QC rules (§E.1).
7. **No outcome association has been tested.** No pathway is prognostic, predictive, a biomarker, or clinically useful on this evidence.
8. **Not externally validated.** No confirmatory cohort exists (D8D-8).
9. **Provenance gaps** — ComBat fitting population and model matrix undocumented (D8D-5).
10. **U1 / U2 unresolved** (D8D-7).

---

## I. Decision log — nothing changed in Phase 8D

| Decision | Phase | Status after 8D |
|---|---|---|
| Cohort S3 = 453, all retained | 7C D2 | **UNCHANGED** |
| M2 primary = 4-set; sensitivity = 5-set | 6C | **UNCHANGED** |
| PI3K sensitivity-only, immutability clause | M2-A1 | **UNCHANGED** |
| σ̂ `ddof = 1` | M2-A1 | **UNCHANGED** |
| Canonical scorer `src/pathway_scoring.py` | 8C rec. 4 | **UNCHANGED** |
| `score_standardisation: NONE` | 6C | **UNCHANGED** — D8D-1 open |
| Median representation primary | 2B | **UNCHANGED** |
| Complete-case missing data | A-19 | **UNCHANGED** |
| QC specification and decisions | 7B / 7C | **UNCHANGED** |
| Endpoint definition | 6A | **NOT ACCESSED** |

**Phase 8C recommendations 1–3 were not implemented.** `reports/phase8d/phase8d_summary.json`
records `phase8c_recommendations_implemented_here: []`.

**Phase 8D is complete. Phase 9 has not been started.**
