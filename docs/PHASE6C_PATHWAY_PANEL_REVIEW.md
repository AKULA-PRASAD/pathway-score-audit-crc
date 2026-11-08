# Phase 6C — M2 Pathway Panel Review (rev. 3, **FROZEN** · spec revision 2)

**Date:** 2026-11-02 · **Revision 2** — final specification audit.
**Amended:** 2026-11-02 — **M2-A1**, σ̂ pinned to `ddof=1` (see §0 below). Membership and roles unchanged.
**Status:** ✅ **APPROVED AND FROZEN** — 2026-11-02T21:13:12.824057+00:00
**Outcome-blind:** no recurrence, survival or outcome column was opened at any point. Panel selection used **only** gene-set metadata, published biology, and the expression correlation structure.

**Artifacts:** `config/pathways/m2_candidate_panel.json` · `reports/phase6c/` (GMT files, mapping, overlap, collinearity, provenance)

---

## ⛔ M2 IS FROZEN

**Frozen:** `2026-11-02T21:13:12.824057+00:00` · **outcome-blind** · **frozen BEFORE any recurrence or survival association analysis.**

| Role | Panel | Sets | Union genes | SHA-256 (rev. 2, current) | SHA-256 (rev. 1, superseded by M2-A1) | File |
|---|---|---:|---:|---|---|---|
| **PRIMARY** | **4-set** | 4 | 607 | `cb4de071e2f64a6b…` | `bab99da9f139db34ccc2d5f5…` | `config/pathways/m2_panel_4set.json` |
| **SENSITIVITY** | 5-set | 5 | 702 | `10dae4874a918160…` | `6058f733b2acda1c7a77746a…` | `config/pathways/m2_panel_5set.json` |

The rev. 1 hashes are **not** retired — they remain the correct attestation for every
artifact produced before the amendment (`A19_FREEZE_RECORD.json`,
`QC_DECISIONS_RECORD.json`), and `src/m2_amendment.py` resolves them against the
amendment chain rather than rewriting them.

---

## 0. Amendment M2-A1 — σ̂ degrees of freedom pinned to `ddof=1`

**Date:** 2026-11-02 · **Approved by:** project owner · **Basis:** Phase 8A §3.5 / §4.4 open item **U5**
**Specification revision:** 1 → **2**

| | |
|---|---|
| **Previous state** | **AMBIGUOUS.** The frozen rule specified "per-gene μ̂/σ̂ from training fold only" but did **not** pin the denominator of σ̂. `numpy.std` defaults to `ddof=0` (population SD); `pandas.DataFrame.std` defaults to `ddof=1` (sample SD). Two conforming implementations could produce different numbers from the same frozen specification and the same data. |
| **Approved convention** | **`ddof = 1`** — sample standard deviation, denominator `n_train − 1`. |
| **Rationale** | (1) consistency with Phase 8A, which used `ddof=1` throughout; (2) explicit reproducibility across implementations — the specification must determine the numbers, not the default of whichever library a future implementer happens to call. |
| **Expected numerical consequence** | **Sample**-standard-deviation scaling rather than **population**-standard-deviation scaling. Because every gene in a fold shares the same `n_train`, this is a single constant rescaling of every score in that fold by `√(n/(n−1))` relative to `ddof=0` — measured at **1.0011055835** for n = 453, matching the analytic constant to 10 decimal places. **Rank order and all inter-score correlations are unchanged; numeric score values are not.** |
| **Membership** | **UNCHANGED** — same 4 and 5 sets, same resolved members (198/195/41/194 (+104)), same union counts (607 / 702). |
| **Roles** | **UNCHANGED** — 4-set PRIMARY, 5-set SENSITIVITY; PI3K/AKT/mTOR remains sensitivity-only. |
| **Coverage rules** | **UNCHANGED** — ≥80% of source members **and** ≥20 genes surviving per-fold zero-variance exclusion. |
| **Outcome data used** | **NO.** The ambiguity was found and decided in Phase 8A, which blocked the recurrence endpoint at the I/O layer. The amendment therefore does **not** reclassify M2 as M3-exploratory under the modification policy, which permits "non-outcome corrections … recorded with reason and timestamp, retain the prior version, and re-freeze". |

**Verification.** Phase 8A was re-run after the amendment. All four pathway-score
artifacts are **byte-identical** to the pre-amendment run, confirming the amendment
pins the convention Phase 8A already used and changes no result.

**PRIMARY (4-set):** EMT (198) · E2F_TARGETS (195) · WNT_BETA_CATENIN (41) · P53_PATHWAY (194)
**SENSITIVITY (5-set):** the above **+ PI3K_AKT_MTOR_SIGNALING (104)**

### ⚠️ PI3K/AKT/mTOR — sensitivity-only, permanently

> **`HALLMARK_PI3K_AKT_MTOR_SIGNALING` is a SENSITIVITY-ONLY addition. It must NOT be promoted to the primary panel, and must NOT be removed from the sensitivity panel, on the basis of observed model performance.**
>
> Either change would convert a pre-specified analysis into a post-hoc one, and **reclassifies M2 as M3-exploratory, permanently and irreversibly** (Phase 6B §3.2.4).

### What was used — and what was not

| Used to decide | NOT used |
|---|---|
| effective dimensionality (identical: 3 PCs → 80%, Kaiser 2) | any recurrence or survival association |
| score-matrix conditioning (**2.95** vs 6.95) | any model fit or performance metric |
| smallest eigenvalue (**0.480** vs 0.280) | any outcome column whatsoever |
| count of \|r\| > 0.5 pairs (**0** vs 1) | |
| parameter budget against 136 events | |
| published recurrence-prediction evidence strength | |

### Preserved provenance

| Item | Value |
|---|---|
| Database | MSigDB Human Hallmark **v2026.2.Hs**, released **2026-10-15** |
| GMT symbols SHA-256 | `56bba1300555fb6302844fc8cdca6340…` |
| GMT entrez SHA-256 | `a26578a72aed3171fffce4d2cda8672c…` |
| Mapping | Entrez-primary (MSigDB Entrez → GPL570 `ENTREZ_GENE_ID` → our symbol), symbol fallback |
| Target index | `frozen_gene_to_probeset_median.csv`, 21,815 HGNC-style symbols |
| Scoring | mean of **training-fold**-standardised member genes |
| Freeze record | `config/pathways/M2_FREEZE_RECORD.json` |

**Modification policy:** immutable. Only **non-outcome** corrections (e.g. a database version fix) are permitted — recorded with reason and timestamp, prior version retained, re-frozen. **Any revision informed by outcome data reclassifies M2 as M3-exploratory.**

**Verification:** `tests/test_phase6c_m2.py` — 22 tests across freeze integrity, specification validity, membership reproducibility, outcome-blindness and immutability policy.

---

## 0. FINAL SPECIFICATION AUDIT (rev. 2)

### 0.1 Correction — the "unmatched" enumeration over-reported

Rev. 1 listed e.g. `CCN1`, `CCN2`, `H2AX`, `RACK1` as unmatched. **The audit found this enumeration was wrong, though the counts were always right.**

| | |
|---|---|
| **The bug** | `unmatched` was computed as *MSigDB symbol strings* minus *our resolved symbols*. A gene recovered under a **legacy alias** was therefore reported as unmatched |
| **Example** | MSigDB `CCN1` (Entrez 3491) matched our `CYR61` via the Entrez path. The gene **was** included — only its modern name is absent from our vocabulary |
| Verified | `CCN1→CYR61`, `CCN2→CTGF`, `P3H1→LEPRE1` all present in resolved members |
| **What was correct** | every count. 198/200 for EMT etc. were right throughout |
| **What was wrong** | only the enumeration, which over-reported by up to 2× |

I also tested whether the Entrez bridge itself was dropping genes (rebuilding it from all probe rows rather than one per symbol): **0 genes gained, 20,471 keys either way.** The bridge was sound.

### 0.2 TRUE unmatched — reconciles exactly

| Gene set | source | resolved | **true unmatched** | previously reported | genes |
|---|---:|---:|---:|---:|---|
| EMT | 200 | 198 | **2** | 5 | FZD8, MATN2 |
| E2F_TARGETS | 200 | 195 | **5** | 11 | CIT, MMS22L, NOP56, PMS2, SRSF2 |
| WNT_BETA_CATENIN | 42 | 41 | **1** | 1 | FZD8 |
| P53_PATHWAY | 200 | 194 | **6** | 12 | IFI30, KRT17, MXD4, PVT1, RACK1, TAX1BP3 |
| PI3K_AKT_MTOR | 105 | 104 | **1** | 2 | MAP2K3 |

**source − resolved = true unmatched for every set.** The definition is now explicit in both spec files: *"a source gene with NO match by Entrez AND no match by symbol. Genes recovered under a legacy alias are MATCHED, not unmatched."*

### 0.3 Mapping audit (requested items 3–4)

| Check | Finding |
|---|---|
| **Multiple probes → one gene** | 10,751 genes have >1 probe (max 15, HFE). **Resolved upstream by frozen median aggregation** — one value per gene |
| **One probe → multiple genes** | 2,619 ambiguous probes exist on the array; **0 are in our gene index** (excluded by construction, Phase 2). No probe contributes to two genes |
| **Duplicate membership within a set** | **0** in all five sets |
| **Unmatched identifiers** | §0.2 — enumerated with Entrez IDs |
| **Final unique genes per score** | EMT 198 · E2F 195 · WNT 41 · P53 194 · PI3K 104 |
| **Each gene contributes exactly once** | ✅ **verified.** Expression index is unique (21,815 rows); for every set, `members = unique_members = present_in_matrix = unique_present`. Scores are computed on the **gene-level** matrix, so probe multiplicity cannot re-enter |

### 0.4 Five-set vs four-set

| Metric | **5-set** | **4-set** |
|---|---:|---:|
| sets | 5 | **4** |
| union unique genes | 702 | 607 |
| **max \|r\|** | 0.632 | **0.318** |
| median \|r\| | 0.181 | 0.181 |
| pairs \|r\| > 0.5 | **1** | **0** |
| pairs \|r\| > 0.7 | 0 | 0 |
| **condition number** | 6.95 | **2.95** |
| **smallest eigenvalue** | 0.280 | **0.480** |
| PCs for 80% | 3 | 3 |
| PCs for 90% | 4 | 4 |
| Kaiser (eigenvalues > 1) | 2 | 2 |

5-set correlation matrix (all 461, **descriptive only**):

| | EMT | E2F | WNT | P53 | PI3K |
|---|---:|---:|---:|---:|---:|
| **EMT** | 1.000 | −0.274 | 0.103 | 0.318 | 0.074 |
| **E2F** | −0.274 | 1.000 | 0.094 | 0.182 | 0.475 |
| **WNT** | 0.103 | 0.094 | 1.000 | 0.180 | 0.077 |
| **P53** | 0.318 | 0.182 | 0.180 | 1.000 | **0.632** |
| **PI3K** | 0.074 | 0.475 | 0.077 | **0.632** | 1.000 |

**Biological coverage.** Both panels retain four distinct axes: stromal/mesenchymal (EMT, the CMS4 representative), proliferation (E2F), canonical Wnt driver (WNT), DNA-damage/checkpoint (P53). The 5-set adds growth signalling (PI3K/AKT/mTOR).

**Redundancy.** PI3K is the least independent member: r = 0.632 with P53 and 0.475 with E2F. It is the **only** source of a >0.5 pair in the 5-set panel.

**Effective dimensionality is IDENTICAL** — 3 PCs for 80%, 4 for 90%, Kaiser 2 — in both panels. **Adding PI3K does not add a dimension; it adds a parameter and halves the conditioning.**

### 0.5 RECOMMENDATION — the **4-set** panel as primary

> **Recommended primary M2 panel: EMT · E2F_TARGETS · WNT_BETA_CATENIN · P53_PATHWAY** (`config/pathways/m2_panel_4set.json`)
> **Pre-specified sensitivity: the 5-set panel** (`config/pathways/m2_panel_5set.json`)

**Reasoning — entirely outcome-blind:**

1. **No loss of effective dimensionality.** Both carry 3–4 real dimensions (Kaiser 2). PI3K contributes a parameter, not information.
2. **Materially better conditioning.** Condition number 2.95 vs 6.95; smallest eigenvalue 0.480 vs 0.280. Coefficients are correspondingly more stable and more interpretable.
3. **Eliminates the only collinear pair.** 0 pairs above \|r\| = 0.5, vs 1 (P53↔PI3K at 0.632).
4. **One fewer parameter** against 136 events — 10–13 total rather than 11–14.
5. **Weakest literature support of the five.** PI3K_AKT_MTOR was graded *moderate–strong* for CRC progression but **weak** for recurrence prediction (§4) — the weakest on the recurrence axis of any retained set.

**Reasons 1–4 are structural and were computed without touching any outcome column. Reason 5 is literature-based.** No recurrence or survival variable was consulted at any point in this decision.

**Why keep the 5-set as sensitivity rather than discarding it:** growth signalling is a genuinely distinct biological axis, and pre-specifying both panels now means the choice cannot later be made on the basis of results.

---

## 1. Database verification

**Claim checked, not assumed.** A bogus version string returns HTTP 404, so the 200s below are meaningful.

| Field | Verified value |
|---|---|
| Collection | MSigDB Human **Hallmark (H)** |
| Version | **v2026.2.Hs** ✅ exists |
| **Release date** | **2026-10-15** (server timestamp) — matches "October 2026" |
| Sets in collection | **50** |
| Symbols GMT | `h.all.v2026.2.Hs.symbols.gmt`, 48,694 B, SHA-256 `56bba1300555fb63…` |
| Entrez GMT | `h.all.v2026.2.Hs.entrez.gmt`, 43,213 B, SHA-256 `a26578a72aed3171…` |
| **All 10 candidates present** | ✅ **yes** |

---

## 2. Target index — identifier type determined BEFORE matching

| Property | Finding |
|---|---|
| Source | `frozen_gene_to_probeset_median.csv` (Phase 2B) |
| n | **21,815** |
| **Identifier type** | **HGNC-style gene symbols**, taken from the GPL570 platform record's `Gene Symbol` column |
| Entrez IDs present? | **0** · Ensembl? **0** · RefSeq? **0** |
| **Legacy symbols** | **LOC 1,547 · FLJ 65 · KIAA 89** |

> **These are NOT interchangeable with a 2026 MSigDB vocabulary.** The GPL570 annotation is a snapshot; 1,701 legacy placeholder symbols confirm it predates current HGNC nomenclature. Symbol-only matching therefore under-matches **systematically, not randomly** — it drops exactly the genes whose names were revised.

### 2.1 Two mapping paths compared

A second, symbol-independent path exists: the GPL570 record carries an `ENTREZ_GENE_ID` column.

| Path | Route | Bridge coverage |
|---|---|---|
| **Primary — Entrez** | MSigDB Entrez → GPL570 `ENTREZ_GENE_ID` → our symbol | **20,471 / 21,815 (93.8%)** |
| Fallback — symbol | direct string match | — |

| Gene set | source n | **symbol-path** | **Entrez-path** |
|---|---:|---:|---:|
| WNT_BETA_CATENIN | 42 | 41 (97.6%) | 41 (97.6%) |
| TGF_BETA | 54 | 52 (96.3%) | 52 (96.3%) |
| EMT | 200 | 195 (97.5%) | **198 (99.0%)** |
| KRAS_SIGNALING_UP | 200 | 191 (95.5%) | **198 (99.0%)** |
| PI3K_AKT_MTOR | 105 | 103 (98.1%) | **104 (99.0%)** |
| P53_PATHWAY | 200 | 188 (94.0%) | **194 (97.0%)** |
| E2F_TARGETS | 200 | 189 (94.5%) | **195 (97.5%)** |
| ANGIOGENESIS | 36 | 35 (97.2%) | 35 (97.2%) |
| HYPOXIA | 200 | 189 (94.5%) | **198 (99.0%)** |
| INFLAMMATORY_RESPONSE | 200 | 198 (99.0%) | **200 (100%)** |
| **mean** | | **96.4%** | **98.2%** |

**Entrez-primary with symbol fallback is adopted.** It recovers 1.8 percentage points on average and up to 3.5 (KRAS_UP) — precisely the revised-symbol genes.

### 2.2 Final resolved membership and ambiguities

| Gene set | source | resolved | unmatched |
|---|---:|---:|---:|
| WNT_BETA_CATENIN_SIGNALING | 42 | **41** | 1 |
| TGF_BETA_SIGNALING | 54 | **52** | 2 |
| EPITHELIAL_MESENCHYMAL_TRANSITION | 200 | **198** | 2 |
| KRAS_SIGNALING_UP | 200 | **198** | 2 |
| PI3K_AKT_MTOR_SIGNALING | 105 | **104** | 1 |
| P53_PATHWAY | 200 | **194** | 6 |
| E2F_TARGETS | 200 | **195** | 5 |
| ANGIOGENESIS | 36 | **35** | 1 |
| HYPOXIA | 200 | **198** | 2 |
| INFLAMMATORY_RESPONSE | 200 | **200** | 0 |

**Mapping ambiguities (one Entrez → several of our symbols): 0.** Every unmatched gene is **enumerated** in `m2_candidate_panel.json` under `unmatched`, never silently dropped.

---

## 3. Overlap and collinearity

### 3.1 Gene-level overlap — low, as Hallmark intends

Max pairwise **Jaccard = 0.054** (EMT ↔ Angiogenesis, 12 shared genes). Union of all ten = **1,238 genes**; only **143 (11.6%)** appear in more than one set. Hallmark sets are curated for non-redundancy and the data confirm it.

### 3.2 Score-level correlation — **severe, and this is what matters**

Gene-level disjointness does **not** imply score independence: co-regulated biology produces correlated scores from disjoint genes. Measured directly (outcome-blind, descriptive mean-of-z scores on the existing bundle — **not** the Phase 6B fold-wise scores):

| Pair | r |
|---|---:|
| **EMT ↔ ANGIOGENESIS** | **+0.914** |
| **KRAS_UP ↔ INFLAMMATORY** | **+0.892** |
| KRAS_UP ↔ EMT | +0.800 |
| KRAS_UP ↔ ANGIOGENESIS | +0.785 |
| HYPOXIA ↔ INFLAMMATORY | +0.770 |
| KRAS_UP ↔ HYPOXIA | +0.765 |
| HYPOXIA ↔ ANGIOGENESIS | +0.738 |
| EMT ↔ HYPOXIA | +0.736 |

**9 pairs exceed \|r\| = 0.7; 18 exceed 0.5.** Condition number of the score correlation matrix: **77.0**; smallest eigenvalue 0.067.

### 3.3 Effective dimensionality

| PCs | cumulative variance |
|---:|---:|
| 1 | **51.4%** |
| 2 | 69.5% |
| 3 | **80.4%** |
| 5 | 91.9% |

**Eigenvalues > 1 (Kaiser): 3.** Ten scores carry roughly **3–5 independent dimensions.**

### 3.4 Correlation clustering — and what it corresponds to

Average linkage on 1 − \|r\|, cut at \|r\| ≥ 0.6 (identical at ≥ 0.5):

| Cluster | Members |
|---|---|
| **C-I** | **TGF_BETA · EMT · KRAS_UP · ANGIOGENESIS · HYPOXIA · INFLAMMATORY** |
| C-II | PI3K_AKT_MTOR · P53_PATHWAY |
| C-III | E2F_TARGETS |
| C-IV | WNT_BETA_CATENIN |

> **C-I is not an artifact — it is the CMS4 mesenchymal subtype.** Guinney et al. describe CMS4 as showing *"epithelial mesenchymal transition (EMT), activation of TGFβ, angiogenesis, matrix remodeling pathways and complement inflammatory system"* with *"a gene expression profile compatible with stromal infiltration."* Six of the ten proposed sets are measuring **one** biological axis.

**This convergence of an unsupervised correlation structure with an independently published subtype is the strongest argument for reducing the panel — and it required no outcome data.**

---

## 4. Literature support — progression vs recurrence-prediction evidence

**The distinction requested is maintained:** evidence that a pathway is implicated in CRC *progression biology* is weaker than evidence that it *predicts recurrence* in stage II–III.

| Gene set | CRC progression biology | Direct recurrence-prediction evidence | Key source |
|---|---|---|---|
| **EMT** | **Strong** — defining feature of CMS4 | **Strong** — CMS4 shows *"clearly worse progression free and overall survival in stage II-III"* | Guinney 2015 |
| **TGF_BETA** | **Strong** — CMS4 defining | **Indirect** — via CMS4, not independently validated as a standalone predictor | Guinney 2015 |
| **ANGIOGENESIS** | **Strong** — CMS4 defining | **Indirect** — via CMS4 | Guinney 2015 |
| **INFLAMMATORY** | **Strong** — complement/inflammatory in CMS4; immune infiltration prognostic in CRC | **Indirect** — via CMS4 | Guinney 2015 |
| **HYPOXIA** | **Strong** | **Moderate–direct** — hypoxia-related signatures reported as independent prognostic factors in stage II/III CRC | hypoxia-signature literature (§8) |
| **KRAS_SIGNALING_UP** | **Strong** — KRAS mutation is a canonical CRC driver | **Weak for the transcriptional signature.** KRAS *mutation* prognosis ≠ KRAS-UP *signature* prognosis | — |
| **WNT_BETA_CATENIN** | **Very strong** — APC/Wnt is the initiating event in most CRC | **Weak.** An initiating event is near-universal and therefore a poor discriminator *between* CRC tumours | — |
| **E2F_TARGETS** | **Strong** — proliferation | **Moderate** — proliferation features recur in published CRC recurrence signatures | recurrence-signature literature (§8) |
| **P53_PATHWAY** | **Strong** — TP53 mutation frequent | **Weak–moderate** as a transcriptional signature | — |
| **PI3K_AKT_MTOR** | **Moderate–strong** | **Weak** | — |

**Honest summary:** only **EMT** has strong evidence on *both* axes. Several candidates are well-supported as progression biology but have **no independent evidence as recurrence predictors** — and WNT, despite being the most canonical CRC pathway, is arguably a *poor* discriminator precisely because it is near-universal.

**None of this used GSE39582 outcomes.**

---

## 5. Is a 10-set panel justified at 136 events?

**No.** Three independent reasons, none outcome-derived:

1. **Parameter budget.** M1 contributes ~6–9 parameters. Adding 10 pathway scores gives 16–19 against 136 events. Phase 5 established that sample size here is a binding constraint.
2. **The panel carries 3–5 real dimensions, not 10** (§3.3). Six of the ten sets measure one axis. Ten correlated scores spend parameters without adding information.
3. **Condition number 77** with nine \|r\| > 0.7 pairs makes individual coefficients unstable and uninterpretable, even under penalisation. Ridge will stabilise prediction but will distribute one effect arbitrarily across six collinear scores — any per-pathway interpretation becomes meaningless.

---

## 6. RECOMMENDED PANEL — 5 sets *(recommendation only; NOT frozen)*

One representative per correlation cluster, chosen on **biology and literature**, never on outcome association:

| # | Gene set | Cluster | n resolved | Rationale |
|---|---|---|---:|---|
| **1** | **HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION** | C-I representative | **198** | Canonical CMS4 marker; **the only candidate with strong progression AND recurrence evidence** |
| **2** | **HALLMARK_E2F_TARGETS** | C-III | **195** | Proliferation axis, **independent** of C-I (r ≤ 0.47, several negative) |
| **3** | **HALLMARK_WNT_BETA_CATENIN_SIGNALING** | C-IV | **41** | Lowest correlation with everything (max \|r\| 0.18); canonical CRC driver. Retained for **biological completeness and independence**, with its weak discrimination rationale stated |
| **4** | **HALLMARK_P53_PATHWAY** | C-II | **194** | DNA-damage/checkpoint axis |
| **5** | **HALLMARK_PI3K_AKT_MTOR_SIGNALING** | C-II | **104** | Growth-signalling axis; r = 0.63 with P53 — the only retained pair above 0.5 |

**Dropped (5):** TGF_BETA, ANGIOGENESIS, HYPOXIA, INFLAMMATORY, KRAS_SIGNALING_UP — **all C-I members represented by EMT** (r = 0.60–0.91 with it).

**Resulting budget:** ~6–9 clinical + 5 pathway = **11–14 parameters** against 136 events.

**Alternative 4-set option:** drop PI3K_AKT_MTOR (the only retained r > 0.5 pair). More conservative; loses a distinct growth-signalling axis.

> **Explicit statement: this panel was selected using gene-set metadata, published CRC biology, and the outcome-blind expression correlation structure. No recurrence, survival or outcome variable was consulted. Should this panel later be revised using outcome information, M2 is permanently reclassified as M3-exploratory** (Phase 6B §3.2.4).

---

## 7. Scoring method

### 7.1 Options assessed

| Method | Mechanism | **Fold-safe?** |
|---|---|---|
| **ssGSEA** (Barbie 2009) | rank-based ECDF difference; scores conventionally normalised **across the sample set** | ❌ **No** — score depends on which samples are present |
| **GSVA** (Hänzelmann 2013) | kernel density estimation of the expression distribution **across samples**, then rank statistic | ❌ **No** — explicitly sample-set-dependent by design |
| **z-score / PLAGE** (GSVA pkg) | standardise across samples, then combine / SVD | ❌ **No** — same dependence; PLAGE adds an SVD rotation |
| **singscore** | rank-based, single-sample | ⚠️ ranks are within-sample, but normalisation references a sample set |
| **Mean of standardised member genes** | per-gene centre/scale, then average over members | ✅ **Yes** — the only estimated quantities are per-gene mean and SD, which project cleanly from train to test |

### 7.2 FINALIZED scoring specification

**Method: mean of training-fold-standardised member genes.**

```
score_P(i) = mean over { g ∈ P ∩ index : σ̂_g^train > 0 } of  (x_gi − μ̂_g^train) / σ̂_g^train
```

| Element | Specification |
|---|---|
| **μ̂_g, σ̂_g** | estimated on **TRAINING-FOLD SAMPLES ONLY** |
| **Application to held-out samples** | the **training-derived** μ̂ and σ̂ are applied unchanged. Held-out data never contributes to either |
| **Zero-variance genes** | σ̂_g^train = 0 → gene **excluded for that fold**; the exclusion is recorded per fold |
| **Missing gene members** | genes absent from the frozen index are **excluded and enumerated**, never imputed |
| **Minimum set coverage** | **≥80% of source members must resolve** AND **≥20 genes must remain** after per-fold zero-variance exclusion. Otherwise the score is **dropped for that fold** and the drop recorded. *(All four primary sets pass comfortably: 97.0–99.0% resolved, 41–198 genes)* |
| **Score standardisation** | **none** before the Cox model — penalisation handles scale. If added later, training-fold-only |
| **Interpretation** | the mean standardised expression of the set's members, **in training-fold SD units**. Positive = coordinated above-average expression relative to the training fold |
| **Each gene once** | guaranteed by construction — see §0.3 |

### 7.3 Why training-derived parameters do NOT leak

> **Leakage is the flow of information from held-out samples into anything fitted.** It is not about whether held-out samples are *transformed*; it is about whether they *inform* the transformation.

| | |
|---|---|
| **μ̂_g, σ̂_g are estimated from training samples only.** | Held-out outcomes, and held-out *expression values*, contribute **nothing** to them |
| **Applying them to held-out data is prediction, not fitting.** | Exactly as applying a fitted Cox coefficient to a held-out patient is prediction. The parameter already exists; the held-out sample is passed through it |
| **The held-out sample does not change the parameter.** | If a held-out sample were added, μ̂ and σ̂ would be unchanged — they are functions of the training fold alone |
| **This is what GSVA/ssGSEA cannot do** (§7.1) | Their scores are defined *relative to the sample set in the matrix*, so a held-out sample's score depends on which other samples accompany it. There is no fixed parameter to carry across — which is why they are rejected here |

**The test to apply:** *if I added one more held-out sample, would any fitted quantity change?* For the z-mean: **no**. For GSVA/ssGSEA: **yes** — and that is leakage.

### 7.4 Fold-internal collinearity assessment

Score correlation and conditioning are **re-assessed inside each training fold**, because the correlation structure is itself an estimate:

| Diagnostic | Computed where | Reported as |
|---|---|---|
| pairwise score correlation matrix | **training fold only** | distribution across folds (median, min, max per pair) |
| condition number of the score matrix | **training fold only** | distribution across folds |
| smallest eigenvalue | **training fold only** | distribution across folds |
| variance inflation factors | **training fold only** | distribution across folds |

**Pre-specified thresholds for flagging** (diagnostic, not automatic action): median across folds of max \|r\| > 0.8, or condition number > 30. **Exceeding a threshold triggers a reported warning, not an automatic panel change** — changing the panel in response to in-fold diagnostics computed alongside outcome modelling would break pre-specification.

### 7.5 Descriptive vs fold-specific — kept strictly separate

| | **§0.4 / §3.2 correlations** | **§7.4 diagnostics** |
|---|---|---|
| Samples | **all 461** (S0 superset) | **training fold only** |
| Purpose | **panel review** — choosing between candidate panels | **modelling diagnostic** |
| Outcome used | **none** | none (diagnostics only) |
| Enters a model? | **NO** | no — reported alongside |
| Status | **descriptive, one-off, pre-modelling** | **per-fold, repeated** |

> **No quantity computed on all 461 samples enters any model, any fold, or any score.** The all-461 matrix exists solely to compare the 5-set and 4-set panels before modelling begins. The scores actually used are computed fold-wise from training data only.

## 8. Reference table

| # | Source | Identifier | Used for |
|---|---|---|---|
| R1 | Liberzon A et al. *The Molecular Signatures Database Hallmark gene set collection.* Cell Syst 2015;1(6):417–425 | [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004) | Hallmark collection design and non-redundancy rationale |
| R2 | MSigDB Human Hallmark **v2026.2.Hs**, released 2026-10-15 | [h.all.v2026.2.Hs.symbols.gmt](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.2.Hs/h.all.v2026.2.Hs.symbols.gmt) · SHA-256 `56bba130…` | **Authoritative membership source** |
| R3 | Guinney J et al. *The consensus molecular subtypes of colorectal cancer.* Nat Med 2015;21(11):1350–1356 | [10.1038/nm.3967](https://doi.org/10.1038/nm.3967) · PMC4636487 | CMS4 = EMT + TGFβ + angiogenesis + inflammatory + stroma; worse RFS in stage II–III |
| R4 | Hänzelmann S, Castelo R, Guinney J. *GSVA: gene set variation analysis for microarray and RNA-seq data.* BMC Bioinformatics 2013;14:7 | [10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7) | GSVA method; sample-set dependence |
| R5 | Barbie DA et al. *Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1.* Nature 2009;462:108–112 | [10.1038/nature08460](https://doi.org/10.1038/nature08460) | ssGSEA origin |
| R6 | GSVA Bioconductor package, v2.4.9 (Bioc 3.22) | [bioconductor.org/packages/GSVA](https://bioconductor.org/packages/GSVA) | Version check for rejected option |
| R7 | Salazar R et al. *Gene expression signature to improve prognosis prediction of stage II and III colorectal cancer.* J Clin Oncol 2011;29(1):17–24 | [10.1200/JCO.2010.30.1077](https://doi.org/10.1200/JCO.2010.30.1077) | ColoPrint; recurrence-prediction precedent in the target population |

**Honest limitation on §4:** the per-pathway evidence grades are my reading of the cited literature plus targeted searches. **R3 and R7 were verified to title/journal/DOI level; the hypoxia- and recurrence-signature claims in §4 rest on search-result summaries whose full texts I did not retrieve.** They are marked "moderate" rather than "strong" for that reason.

---

## 9. Status and next step

| Item | Status |
|---|---|
| Database + version verified | ✅ |
| Identifier type determined before matching | ✅ |
| Mapping completed, ambiguities enumerated | ✅ 0 ambiguities |
| Overlap + collinearity quantified | ✅ |
| Literature reviewed, progression vs recurrence separated | ✅ |
| Scoring method recommended | ✅ z-mean, fold-safe |
| Software identified | ✅ none required |
| **Panel frozen** | ❌ **NOT FROZEN — awaiting approval** |

**On approval**, the five Phase 6B fields A-13/A-14 are filled by writing the approved panel and scoring specification to a timestamped, version-controlled file, after which it is immutable except under the non-outcome revision rule (Phase 6B §3.2.4).
