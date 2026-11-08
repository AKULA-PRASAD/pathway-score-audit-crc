# Phase 8C — Frozen Pathway-Score Analysis

**Version:** 1.0 · **Status:** executed · **Date:** 2026-11-02
**Script:** `scripts/34_phase8c_pathway_scores.py` · **Tests:** `tests/test_phase8c_pathway_scores.py`
**Outputs:** `reports/phase8c/` — 12 CSV + 1 JSON + 3 figures
**Cohort:** S3, n = 453, **all retained** (Phase 7C D2, unchanged)
**Specifications:** `m2_panel_4set.json` / `m2_panel_5set.json`, **revision 2** (amendment **M2-A1**, `sigma_ddof = 1`)

> **This is not a modelling phase.** No outcome was read. No pathway is claimed to be
> prognostic, predictive, a biomarker, or clinically useful — **no outcome association has been
> tested.** No panel change was made or is proposed.

---

## 1. Governing records inspected

`config/pathways/m2_panel_4set.json` (rev 2) · `m2_panel_5set.json` (rev 2) ·
`M2_FREEZE_RECORD.json` incl. amendment M2-A1 · `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md`
(§3.2.2 fold-specific computation, §4 leakage rules) · `docs/PHASE6C_PATHWAY_PANEL_REVIEW.md` ·
`scripts/32_phase8a_eda.py` (`score_panel`) and its artifacts · `docs/PHASE8A_EXPLORATORY_ANALYSIS.md` ·
`docs/PHASE8B_PCA_BATCH_INVESTIGATION.md` · `config/analysis/QC_SPECIFICATION.json` (rev 3) ·
`config/analysis/QC_DECISIONS_RECORD.json` · `data/processed/frozen/frozen_representation.json`.

**There is no separate "Phase 8 roadmap" document** in `docs/`; the forward plan lives in the
Phase 6B modelling protocol and in the open-questions sections of the 8A and 8B reports. Noted
rather than assumed.

The current frozen specifications are the sole source of pathway membership. Nothing was
reconstructed from memory or from another MSigDB release.

### 1.1 Input artifacts

| Artifact | SHA-256 (16) |
|---|---|
| `config/pathways/m2_panel_4set.json` | `cb4de071e2f64a6b` |
| `config/pathways/m2_panel_5set.json` | `10dae4874a918160` |
| `config/pathways/M2_FREEZE_RECORD.json` | `c2b31f3340b36355` |
| `data/processed/analysis/expression_median_S0.parquet` | `a1ca41807c534584` |
| `data/processed/analysis/expression_repiqr_disc443_S0.parquet` | `e636e259b4feabfa` |
| `data/processed/analysis/cohort_membership.csv` | `bd3671f8df64c594` |
| `data/processed/analysis/partition_and_batch.csv` | `57a45bab8efb81f9` |
| `data/processed/analysis/covariates.csv` | `47b2d9829a4fccec` |
| `metadata/GSE39582_sample_metadata.csv` | `37cef7feb767766f` |

### 1.2 Outcome blindness — enforced at the loading boundary

| Layer | Mechanism |
|---|---|
| **File** | `_guard_outcome_files()` patches `pandas.read_csv/read_parquet/read_table` to raise `PermissionError` on any path containing `endpoint.csv` |
| **Column** | `load_metadata_outcome_blind()` drops `rfs.event`, `rfs.delay`, `os.event`, `os.delay_(months)` at read and asserts their absence before returning |

Tests prove both, including a poisoning test that corrupts every outcome value in the source and
demands bit-identical output (§6).

---

## 2. A — Reproduction and audit of the frozen scores

### 2.1 Membership and coverage — verified against the current specifications

| Gene set | source | resolved | present in matrix | used | zero-var dropped | coverage | rule |
|---|---:|---:|---:|---:|---:|---|---|
| EPITHELIAL_MESENCHYMAL_TRANSITION | 200 | 198 | 198 | **198** | 0 | 99.0 % | MET |
| E2F_TARGETS | 200 | 195 | 195 | **195** | 0 | 97.5 % | MET |
| WNT_BETA_CATENIN_SIGNALING | 42 | 41 | 41 | **41** | 0 | 97.6 % | MET |
| P53_PATHWAY | 200 | 194 | 194 | **194** | 0 | 97.0 % | MET |
| PI3K_AKT_MTOR_SIGNALING *(5-set only)* | 105 | 104 | 104 | **104** | 0 | 99.0 % | MET |

Membership matches `M2_FREEZE_RECORD.json` for both panels. Both carry
`specification_revision: 2` and `sigma_ddof: 1`. **Per-sample missingness is 0** and the matrix
holds **0 non-finite values**, so every sample uses the identical gene set — there is no
per-sample coverage variation to report.

### 2.2 Audit checks

| Check | Result |
|---|---|
| Production scorer vs an **independent naive re-implementation** | max \|diff\| **3.0 × 10⁻¹⁵** |
| vs the Phase 8A score artifacts (CSV rounded to 6 dp) | max \|diff\| **5.0 × 10⁻⁷** |
| `ddof = 1` genuinely in force (ddof 0/1 ratio vs √(n/(n−1))) | 1.0011055835**330** vs 1.0011055835**328** |
| Fold isolation under **deletion** of held-out samples | **5/5 bit-identical** |
| Fold isolation under **+1000 poisoning** of held-out samples | **5/5 bit-identical**, and held-out scores do move |
| Zero-variance path exercised (constant gene injected) | 1 injected → **1 dropped** (observed in real data: 0) |
| 4-set vs 5-set: the four shared scores | **max \|diff\| = 0.0, exactly** |

The independent re-implementation uses explicit Python loops and no pandas broadcasting, so
agreement at 3 × 10⁻¹⁵ is evidence that neither path carries a bug the other would cancel. The
5 × 10⁻⁷ figure against Phase 8A is the 6-decimal rounding of the stored CSV, not a discrepancy —
**no numerical discrepancy was found, and no Phase 8A artifact was overwritten** (Phase 8C writes
only to `reports/phase8c/`).

The last row is a structural invariant worth stating: a score depends only on its own members, so
adding a fifth set **cannot** move the other four. It is exactly zero, confirmed.

---

## 3. B — Score distributions and behaviour

| Pathway | k | SD | IQR | skew | kurtosis | SD if members independent | implied mean inter-gene r̄ |
|---|---:|---|---|---|---|---|---|
| EPITHELIAL_MESENCHYMAL_TRANSITION | 198 | 0.4950 | 0.6404 | +0.59 | +0.56 | 0.0711 | **0.241** |
| E2F_TARGETS | 195 | 0.5543 | 0.7207 | −0.69 | +1.09 | 0.0716 | **0.304** |
| WNT_BETA_CATENIN_SIGNALING | 41 | 0.2531 | 0.3357 | +0.22 | −0.04 | 0.1562 | **0.041** |
| P53_PATHWAY | 194 | 0.1949 | 0.2396 | +0.45 | +0.97 | 0.0718 | **0.033** |
| PI3K_AKT_MTOR_SIGNALING | 104 | 0.1851 | 0.2357 | +0.09 | +0.53 | 0.0981 | **0.025** |

All means are 0 by construction under a full-cohort fit. Full quantiles (p0, p1, p5, p25, p50,
p75, p95, p99, p100) are in `score_distributions.csv`.

### 3.1 Separating the data from the transformation — **new in 8C**

Averaging *k* unit-variance z-scores yields SD = 1/√k **if members were independent**. Inverting
`Var = (1 + (k−1)·r̄)/k` recovers the implied mean inter-gene correlation r̄, which separates
*set coherence* from the shrinkage that averaging necessarily imposes.

**Direct observation.** The five sets fall into two regimes. E2F_TARGETS (r̄ = 0.304) and EMT
(0.241) are internally coherent: their members move together, and their scores carry 7–8× the
dispersion independence would give. WNT (0.041), P53 (0.033) and PI3K (0.025) are **nearly
incoherent**: P53 averages 194 genes whose mean pairwise correlation is ≈ 0.03, producing an SD
only 2.7× the independent-member floor.

**Plausible explanation.** The HALLMARK P53 and PI3K/AKT/mTOR sets are broad, pathway-membership
collections rather than tight co-expression modules, and a mean-of-z-scores summary rewards
coordinated sets. **Unresolved:** whether low coherence reflects the gene sets, this tissue, the
deposited preprocessing, or all three. This is a property of **the score as constructed**, not a
statement about recurrence.

### 3.2 Variation and influence

All five scores have non-zero, non-degenerate dispersion across 453 samples. **Whether that
dispersion is *sufficient* for a model is an outcome-dependent question and is not answered
here.**

Leave-one-out influence on score SD — **reported, never an exclusion rule**:

| Pathway | max ΔSD | % of SD | sample | z |
|---|---|---|---|---|
| EMT | 0.00999 | 2.02 % | GSM972479 | +4.36 |
| E2F_TARGETS | 0.01183 | 2.13 % | GSM972357 | −4.47 |
| WNT_BETA_CATENIN | 0.00294 | 1.16 % | GSM972371 | +3.38 |
| P53_PATHWAY | 0.00311 | 1.60 % | GSM972275 | +3.91 |
| PI3K_AKT_MTOR | 0.00233 | 1.26 % | GSM971993 | +3.50 |

No single sample moves any score's SD by more than **2.2 %**. GSM972479 and GSM972357 also appear
in the Phase 8A/8B observations (sex-axis discordance; low mean correlation). **No exclusion
criterion is created from any of this.**

---

## 4. C — Relationships among pathways

| | 4-set (8C, n=453) | 6C recorded (n=461) | 5-set (8C) | 6C recorded |
|---|---|---|---|---|
| max \|r\| | 0.316 | 0.318 | **0.631** | 0.632 |
| median \|r\| | 0.186 | 0.181 | 0.186 | — |
| n(\|r\| > 0.5) | **0** | 0 | **1** | 1 |
| condition number | 2.949 | 2.950 | 7.032 | 6.950 |
| min eigenvalue | 0.4811 | 0.4798 | 0.2777 | 0.2801 |
| PCs for 80 % | 3 | 3 | 3 | 3 |
| Kaiser | 2 | 2 | 2 | 2 |

**The Phase 6C and Phase 8A redundancy findings reproduce.** Residual differences are the cohort
(453 vs 461), not the method.

Eigenvalues — 4-set: 1.419, 1.244, 0.857, 0.481 · 5-set: 1.953, 1.337, 0.943, 0.489, 0.278.

The single \|r\| > 0.5 pair in the 5-set is **P53_PATHWAY ↔ PI3K_AKT_MTOR (r = +0.631)**.

### 4.1 Incremental descriptive contribution of PI3K

Regressing the PI3K score on the four primary scores: **R² = 0.540**, so **46.0 %** of its
variance is not linearly captured by them (residual SD 0.126). Adding it raises the condition
number 2.949 → 7.032 and lowers the minimum eigenvalue 0.481 → 0.278.

> This is a **descriptive redundancy measure computed among scores only**. It says nothing about
> recurrence relevance. It is **not** grounds to promote, remove or reorder any pathway: the M2-A1
> immutability clause stands, and PI3K remains **sensitivity-only**.

---

## 5. D — Technical and biological structure in the scores

40 tests (8 variables × 5 scores), each assessed for **location** (Kruskal–Wallis, ε²) **and
dispersion** (Brown–Forsythe) — the Phase 8B lesson — with BH q-values and reported cell sizes.

### 5.1 Two distinct circularity problems, kept apart

| Class | Variable | Why |
|---|---|---|
| **Circular** | `cit.molecularsubtype` | the depositors' own **expression-derived** clustering of this same matrix. An M2 score is also a function of that matrix, so any association is near-tautological and is **not** independent biological validation. |
| **Matrix-derived** | `qc_any_soft_flag` | the Phase 7B S-S1/S-S2/S-S3 flags are **Tukey rules on the sample median, IQR and mean inter-sample correlation of this same matrix**. A sample flagged for an unusual IQR has unusual z-scores **by construction**, so a dispersion association with the scores is partly **mechanical**. |

Both are excluded from the "non-circular" counts and reported separately.

### 5.2 Non-circular location associations (q < 0.05): 17 of 30 tests

| Variable | Pathway | ε² | q (BH) | n | min cell |
|---|---|---|---|---:|---:|
| MMR status | P53_PATHWAY | **0.172** | 4.7e-16 | 410 | 60 |
| `combat_mod` | P53_PATHWAY | 0.151 | 9.9e-16 | 453 | 65 |
| Tumour location | P53_PATHWAY | 0.066 | 1.9e-07 | 453 | 185 |
| MMR status | PI3K_AKT_MTOR | 0.042 | 7.9e-05 | 410 | 60 |
| `combat_mod` | PI3K_AKT_MTOR | 0.040 | 5.4e-05 | 453 | 65 |
| MMR status | EMT | 0.039 | 1.3e-04 | 410 | 60 |
| RNA extraction method | P53_PATHWAY | 0.037 | 1.9e-03 | 453 | 14 |
| ComBat batch | P53_PATHWAY | 0.033 | 8.2e-03 | 453 | 14 |

`combat_mod` is reported as itself; its interpretation remains **unresolved** (Phase 7A Q1) and no
claim is made that it *is* the MMR variable the depositors used.

### 5.3 Conditional re-test — **batch and extraction-method effects do not survive**

Stratified permutation, batch/method labels permuted **within** `combat_mod` strata, 10,000
permutations:

| Group of tests | n tests | result |
|---|---:|---|
| ComBat batch × 5 scores, RNA extraction method × 5 scores | 10 | **0 survive** (min p: location 0.171, dispersion 0.125) |
| `qc_any_soft_flag` × 5 scores (matrix-derived) | 5 | **4 survive** |

| Surviving test | location p | dispersion p |
|---|---|---|
| `qc_any_soft_flag` ~ E2F_TARGETS | 0.051 | **0.0001** |
| `qc_any_soft_flag` ~ P53_PATHWAY | 0.523 | **0.0001** |
| `qc_any_soft_flag` ~ PI3K_AKT_MTOR | 0.597 | **0.0013** |
| `qc_any_soft_flag` ~ EMT | **0.027** | 0.088 |

**Direct observation.** Every batch and extraction-method association with every pathway score
disappears once MMR composition is conditioned on — consistent with Phase 8B's finding for PC
*location*. The only survivors involve the QC soft flags, and they are **dispersion** effects.

**Explanation — partly mechanical, not established.** The flags are Tukey rules on this same
matrix (§5.1), so flagged samples have unusual z-scores by construction. This is reported as an
observation; it is **not** an independent discovery of array failure, and it is **explicitly not
grounds for exclusion** — Phase 7C D2 retains all 453.

**Limitation.** ComBat was applied to this deposit before release, so batch variance was removed
by construction; a null batch result bounds what remains detectable and is **not** evidence that
batch effects are absent. Sparse strata (min cell 14 for batch and extraction method) and
multiplicity make single nominal p-values weak evidence.

---

## 6. E — Sensitivity and stability

### 6.1 Primary 4-set vs sensitivity 5-set

The four shared scores are **identical** (max \|diff\| = 0.0). The panels differ only by the
presence of the PI3K score and the consequent correlation-structure change (§4.1).

### 6.2 Frozen median vs frozen representative-IQR representation — **sensitivity-only**

Both frozen artifacts are 21,815-gene indices over the same 453 samples, so the panels resolve
identically and **only the probe→gene rule differs**.

| Pathway | Pearson | Spearman | SD (median) | SD (rep-IQR) | mean \|diff\| |
|---|---|---|---|---|---|
| EPITHELIAL_MESENCHYMAL_TRANSITION | 0.992 | 0.991 | 0.4950 | 0.5314 | 0.060 |
| E2F_TARGETS | **0.998** | 0.997 | 0.5543 | 0.5648 | 0.032 |
| WNT_BETA_CATENIN_SIGNALING | **0.858** | 0.841 | 0.2531 | 0.2460 | 0.105 |
| P53_PATHWAY | 0.972 | 0.966 | 0.1949 | 0.2028 | 0.038 |
| PI3K_AKT_MTOR_SIGNALING | 0.910 | 0.900 | 0.1851 | 0.1892 | 0.060 |

Correlation structure: max \|Δr\| **0.106**; condition number 7.032 → **9.530**.

**Direct observation.** The two coherent sets are near-invariant to the representation choice;
the **smallest set (WNT, 41 genes) is the least stable at r = 0.858**, and PI3K is next at 0.910.
This tracks set size and coherence: with few, weakly correlated members, each gene's probe choice
matters more.

**Scope limit.** This is a **representation** sensitivity check. The differences are differences
in probe→gene summarisation, **not biological differences**, and the small SD shifts are
score-scale effects, not evidence about any pathway.

### 6.3 Training-fold composition — the Phase 8A folds, **reused not created**

The 5-way outcome-blind split is regenerated from seed 20261102. **No new folds were chosen and
no scoring choice was tuned.**

| Pathway | max \|shift\| vs full-cohort fit | min corr with full-cohort fit |
|---|---|---|
| EPITHELIAL_MESENCHYMAL_TRANSITION | 0.0722 SD | 0.999985 |
| E2F_TARGETS | 0.0440 SD | 0.999981 |
| P53_PATHWAY | 0.0309 SD | 0.999805 |
| WNT_BETA_CATENIN_SIGNALING | 0.0201 SD | 0.999856 |
| PI3K_AKT_MTOR_SIGNALING | 0.0199 SD | 0.999610 |

Scores are highly stable to training-fold composition (all correlations > 0.9996), while the
shift is **non-zero** — which is exactly why the frozen rule keeps standardisation inside the
resampling loop.

---

## 7. Labelled findings

### Reproduced from Phase 6C / 8A
- Membership, resolved counts and coverage for all five sets.
- 4-set: max \|r\| 0.316, condition 2.949, 0 pairs > 0.5. 5-set: max \|r\| 0.631, condition 7.032, 1 pair > 0.5. PCs for 80 % = 3, Kaiser = 2.
- Full-cohort scores match the Phase 8A artifacts to CSV rounding.
- Fold isolation holds; the max held-out shift of 0.0722 SD reproduces.

### New in Phase 8C
- Independent re-implementation agrees with the production scorer to 3 × 10⁻¹⁵.
- Fold isolation verified under **poisoning**, not only deletion.
- Zero-variance exclusion path **exercised** with an injected constant gene.
- 4-set and 5-set shared scores proven **exactly** identical.
- **Implied mean inter-gene correlation** separating set coherence from averaging shrinkage: two coherent sets (0.24, 0.30), three near-incoherent (0.025–0.041).
- Leave-one-out influence: no sample moves any SD by more than 2.2 %.
- PI3K descriptive redundancy: R² = 0.540 on the four primary scores.
- Location **and dispersion** testing of all score × annotation pairs with BH control.
- Conditional permutation re-test: **0 of 10** batch/extraction associations survive MMR stratification.
- `qc_any_soft_flag` identified as **matrix-derived**, with its dispersion associations flagged as partly mechanical.

### Sensitivity-only
- Median vs representative-IQR representation (§6.2) — WNT least stable (r = 0.858).
- The 5-set panel throughout; PI3K remains sensitivity-only.

---

## 8. Limitations and unresolved questions

1. **No outcome association has been tested.** Nothing here supports any claim of prognostic or predictive value.
2. **Low coherence of three sets (§3.1)** is a property of the score construction. Whether it limits their usefulness is an outcome-dependent question for the modelling protocol.
3. **ComBat was applied pre-deposit.** Null batch results bound what remains detectable; they are not evidence of absence.
4. **`combat_mod` remains uninterpreted** — Phase 7A Q1, closable only by author correspondence.
5. **The QC soft-flag dispersion effect is partly mechanical** (§5.1) and cannot be cleanly separated from genuine array-quality effects using this matrix alone.
6. **Sparse strata**: minimum cell 14 for batch and extraction method; nominal p-values from such cells are weak evidence.
7. **`cit.molecularsubtype` is circular** and is excluded from all non-circular counts.
8. **Representation differences are not biological** (§6.2).
9. The production scorer currently lives inside `scripts/32_phase8a_eda.py` and must be imported by path. For the modelling phase it would be safer in `src/` — **recommendation only, not done here**, because moving it risks the reproducibility of the Phase 8A artifacts.
10. **U1 / U2 unchanged** — `rfs.delay` unit `[PROVISIONAL: months]`; death-without-recurrence handling. Untouchable in an outcome-blind phase.

---

## 9. Decision log — nothing changed

| Decision | Phase | Status after 8C |
|---|---|---|
| M2 primary = 4-set (EMT, E2F_TARGETS, WNT_BETA_CATENIN, P53_PATHWAY) | 6C | **UNCHANGED** |
| M2 sensitivity = 5-set, adding PI3K_AKT_MTOR_SIGNALING | 6C | **UNCHANGED** — sensitivity-only |
| PI3K immutability clause | M2-A1 record | **UNCHANGED** — not promoted, not removed |
| σ̂ `ddof = 1` | M2-A1 | **UNCHANGED** — read from spec, verified in force |
| Cohort S3 = 453, all retained | 7C D2 | **UNCHANGED** — 0 excluded, 0 relabelled |
| QC specification and decisions | 7B / 7C | **UNCHANGED** |
| Median gene-level representation (primary) | 2B | **UNCHANGED** — rep-IQR used only as a sensitivity comparison |
| Analysis bundle | 3 | **UNCHANGED** |
| Complete-case missing data | A-19 | **UNCHANGED** |
| Endpoint definition | 6A | **NOT ACCESSED** |

No pathway was promoted, removed, replaced or reordered. No exploratory result was converted into
an exclusion criterion or a specification change.

---

## 10. Recommendations for the modelling protocol — **requiring separate review**

These are recorded as candidates for the Phase 6B protocol review. **None is applied here and
none should be applied without your explicit decision.**

1. **Low-coherence scores (§3.1).** WNT, P53 and PI3K have implied inter-gene r̄ ≤ 0.041 and correspondingly small dispersion. Under a penalised Cox model, a low-dispersion predictor is shrunk more heavily. Whether the protocol should standardise scores before penalisation — the M2 record currently says **no** score-level standardisation, and that is frozen — deserves explicit review *before* any outcome is seen.
2. **PI3K collinearity (§4.1).** Adding PI3K raises the condition number to 7.03 and shares r = 0.631 with P53. The 5-set panel remains the pre-specified sensitivity analysis; the protocol should state how a divergence between the 4-set and 5-set results will be interpreted, decided in advance.
3. **WNT representation sensitivity (§6.2).** At r = 0.858 between the two frozen representations, WNT is the one score whose value depends materially on the probe→gene rule. A pre-specified representation sensitivity analysis is worth considering.
4. **Scorer location (§8.9).** Move the production scorer into `src/` before modelling.

---

**Phase 8C is complete. No Phase 8D or Phase 9 activity has been started.**
