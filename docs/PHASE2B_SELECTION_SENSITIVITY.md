# Phase 2B — Gene-Level Representation: Evidence, Interpretation and Approved Decisions

**Date:** 2026-11-02 · **Revision 2** — reviewed and approved with corrections
**Status:** Representation approved and frozen. See §0.
**Reproduce:**
```
.venv/bin/python scripts/20_phase2b_selection_sensitivity.py
.venv/bin/python scripts/21_phase2b_probe_quality.py
```
Evidence in `reports/phase2b/`. Six candidate matrices now exist; **the original four are untouched.**

**Boundaries honoured:** no normalization, no batch correction, no outcome column read in any Phase 2/2B module, no survival modelling, no feature selection, no final selection. Scripts 01–14 and the Phase 1 manifest unmodified. Probe cache unchanged.

---

## 0. DECISION SUMMARY — read this first

Five categories, deliberately never merged.

### ✅ APPROVED DECISIONS

| # | Decision | Artifact |
|---|---|---|
| D1 | **Probe-selection population = `discovery_443`** — all 443 discovery primary tumours. The 123 validation samples are **not** used for probe selection. S3 exclusions are **not** used to define it. | `data/processed/frozen/selection_population_discovery443.txt` |
| D2 | **Primary representation = `median` aggregation** | `frozen_gene_to_probeset_median.csv` |
| D3 | **Principal sensitivity = `representative_iqr` selected on `discovery_443`** | `frozen_gene_to_probe_discovery443.csv` |
| D4 | **`mean` retained as historical evidence, not a primary strategy** | original candidate preserved |
| D5 | **α-stratified hybrid NOT approved; no α threshold implemented** | — |
| D6 | **Mapping count = 8,893 unmapped** (8,894 counted a SOFT terminator) | §5 |

**Why `discovery_443` is the right selection population.** It is the *development* population. The representation must be defined without consulting the held-out validation population, or the holdout is no longer held out — the gene→probe map would encode information from the samples later used to assess generalisation. `discovery_443` is also defined purely by the GEO `dataset` field (study design), which keeps the selection strictly outcome-blind; `s3_discovery_354` does not (§1.2).

**The frozen map was created before any application to validation samples**, and runtime assertions confirm no validation sample entered the selection (§0.1).

### ⚠️ PROVISIONAL DECISIONS

| # | Decision | Status |
|---|---|---|
| P1 | `median` as the uniform aggregation rule | **Provisional.** Chosen as a single, uniform, reportable rule. **Not claimed to be biologically optimal** — no evidence in this project establishes that, and none could without an independent reference measurement. |

### 🔬 OBSERVED EVIDENCE (measurements, §1–§6)

Selection-set sensitivity · residual batch association · probe intensity and variability · Affymetrix suffix distribution · within-gene probe correlation · α distribution · probe accounting · transportability requirements.

### 📐 STATISTICAL INTERPRETATION (what the measurements support)

Each section states what it does and does not support. The α section (§4) was **materially corrected** in this revision — see §4.0.

### ❓ SCIENTIFIC LIMITATIONS (§9)

No reference measurement exists; residual-batch findings cannot speak to uncorrected cross-cohort effects; annotation-version and probe-availability dependencies are unverified for any specific external cohort.

---

## 0.1 Frozen-representation verification

| Check | Result |
|---|---|
| selection population size == 443 | **PASS** |
| no validation sample in selection population | **PASS** |
| no non-tumoural sample in selection population | **PASS** |
| selection population defined by GEO `dataset` field only | **PASS** |
| no endpoint column read during freezing | **PASS** |
| frozen map covers every uniquely-mapped gene (21,815) | **PASS** |
| median groupings cover every uniquely-mapped probe (43,163) | **PASS** |
| re-selection from discovery_443 is deterministic | **PASS** |
| ties in selection | **0** |

Selection-population SHA-256: `26943e88503f9d32506d07f9570f08d3…`
The frozen discovery map differs from the all-585 map in **542 genes**; the all-585 and s3-354 mappings are **preserved** for audit at `reports/phase2b/representative_probe_by_selection_set.csv`.

---

## 1. Representative-probe choice vs selection set

IQR-based selection run on three sample sets. All three emit all 585 samples; only the *scoring* is scoped.

| Comparison | Genes differing | % of all 21,815 | **% of 10,751 multi-probe** |
|---|---:|---:|---:|
| all_585 vs discovery_443 | 542 | 2.48% | **5.04%** |
| all_585 vs s3_discovery_354 | 735 | 3.37% | **6.84%** |
| discovery_443 vs s3_discovery_354 | 523 | 2.40% | **4.86%** |

Single-probe genes can never differ, so the multi-probe percentage is the meaningful figure. **Selection is ~93–95% stable** across sample sets.

### 1.1 Which genes are affected — a counter-intuitive result

| | median within-gene ρ |
|---|---:|
| genes that **changed** probe | **0.406** |
| genes that **did not** change | **0.328** |

Affected genes have *higher* probe agreement, not lower. This is self-consistent: when a gene's probes agree, their IQRs are close, so the ranking is a near-tie and a modest change in the sample set flips it. **The genes whose selection is unstable are precisely those where the choice matters least** — the competing probes measure the same thing.

Most-affected genes (largest IQR gap between the two winners), from `genes_affected_by_selection_set.csv`:

| gene | k probes | within-gene ρ | all_585 pick | s3_354 pick | IQR(585) | IQR(354) |
|---|---:|---:|---|---|---:|---:|
| CAPN9 | 2 | 0.917 | 208063_s_at | 210641_at | 1.592 | 1.413 |
| AFF4 | 8 | 0.102 | 239439_at | 232865_at | 1.080 | 0.940 |
| ETS2 | 3 | 0.559 | 201328_at | 241193_at | 0.970 | 0.848 |
| RBMS1 | 6 | 0.917 | 209868_s_at | 225265_at | 1.448 | 1.331 |
| PAPSS2 | 3 | 0.892 | 203058_s_at | 203060_s_at | 1.247 | 1.157 |

### 1.2 A finding that bears on which set to freeze on

**`s3_discovery_354` is not strictly outcome-blind.** S3 membership is *defined using the endpoint*: of the 8 discovery samples it drops, 3 are dropped for missing `rfs.event`/`rfs.delay`, **2 because `rfs.event == 1`, and 3 because `rfs.event == 0`**.

Selecting probes on S3 therefore conditions the expression representation on outcome values. The effect is tiny (8 of 362 samples) but the principle is not: it is exactly the kind of leak that is invisible later. **`discovery_443` is defined purely by the GEO `dataset` field — study design, not outcome — and is strictly outcome-blind.**

---

## 2. Batch association of selected probes (descriptive only)

ComBat batch labels recovered from the ArrayExpress SDRF (**566 of 585** samples; the 19 non-tumoural mucosa samples are absent from ArrayExpress). Nine batches. **No correction performed.**

Variance fraction explained by batch (one-way η², per probe, multi-probe genes):

| | n | median η² | mean η² | p90 |
|---|---:|---:|---:|---:|
| **selected** probes | 10,751 | **0.0050** | 0.0080 | 0.0184 |
| **not-selected** probes | 21,348 | **0.0040** | 0.0065 | 0.0146 |

η² by IQR decile — monotone increasing: 0.0021 (decile 1) → 0.0063 (decile 9).

| | Spearman |
|---|---:|
| IQR vs batch η² | **+0.342** |
| mean expression vs batch η² | +0.257 |

| | η² > 0.25 |
|---|---:|
| selected probes | **0.0%** |
| not-selected probes | **0.0%** |

**Answer to the question asked: yes, directionally — high-IQR probes *are* disproportionately batch-associated (ρ = +0.342). But the absolute magnitude is negligible: median 0.5% of variance, and not a single probe exceeds 25%.** On these data, IQR selection is not meaningfully picking up batch structure.

**The critical caveat.** These values are measured on data the depositors **already ComBat-corrected**. Residual batch variance is small *by construction* — ComBat was fitted to remove exactly this. This analysis therefore says:

- ✅ IQR selection is not chasing *residual* batch structure in this deposit.
- ❌ It says **nothing** about how batch-sensitive these probes would be in uncorrected data.
- ❌ It says **nothing** about cross-cohort batch behaviour, which is the Phase 4 question.

A +0.342 rank correlation on a near-zero quantity should not be read as reassurance about transport.

**Stated explicitly, as required: residual batch association in already-ComBat-processed data is NOT evidence that uncorrected cross-cohort batch effects are negligible.** ComBat was fitted on these very arrays to remove between-batch variation; finding little left is the expected consequence of that fit, not an independent finding about batch sensitivity. An external cohort will have been processed by a different laboratory, on a different schedule, with its own normalization — none of which this measurement observes. Cross-cohort batch behaviour remains entirely unmeasured and is a Phase 4 obligation.

---

## 3. Probe quality: selected vs not-selected (multi-probe genes)

| Measure | Selected (median) | Not selected (median) |
|---|---:|---:|
| mean expression | **6.062** | 4.894 |
| IQR | 0.688 | 0.433 |
| variance | 0.290 | 0.120 |
| mean sibling ρ | **0.395** | 0.348 |

IQR and variance are higher **by construction** — that is the selection criterion, not evidence of quality.

**Low-expression frequency:**

| Threshold | Selected | Not selected |
|---|---:|---:|
| mean log2 < 4.0 | **13.7%** | 33.3% |
| mean log2 < 5.0 | **29.4%** | 51.8% |
| mean log2 < 6.0 | **48.8%** | 66.1% |

**Annotation specificity** (Affymetrix suffix: `_at` one transcript; `_s_at` documented transcript family; `_x_at` unpredictable cross-hybridisation):

| suffix | not-selected % | selected % | n total |
|---|---:|---:|---:|
| `_at` | 60.18 | 59.36 | 19,229 |
| `_s_at` | 25.96 | **31.38** | 8,915 |
| `_x_at` | 9.94 | **5.65** | 2,729 |
| `_a_at` | 3.81 | 3.55 | 1,196 |

**This is the most reassuring result in Phase 2B.** The obvious failure mode of an IQR criterion is preferentially selecting cross-hybridising probes, which are variable for technical reasons. **The opposite happens: `_x_at` probes are selected at 5.65% versus 9.94% among non-selected — roughly half the rate.** IQR selection is enriched for `_s_at` rather than `_x_at`.

**Stated explicitly, as required: variability is not biological quality.** A cross-hybridising or outlier-prone probe can be highly variable for purely technical reasons, and an IQR criterion *cannot distinguish that from real biology*. The suffix evidence shows this failure mode is not dominant here; it does not show that IQR measures biological quality. Sibling-ρ is only marginally better for selected probes (0.395 vs 0.348), so selected probes are **not** markedly more representative of their gene's consensus.

---

## 4. Probe agreement — Cronbach's α as descriptive evidence only

### 4.0 Correction to revision 1

**Revision 1 of this report overstated what α establishes, and that interpretation is withdrawn.** Specifically withdrawn:

| Withdrawn claim (rev. 1) | Why it was wrong |
|---|---|
| "Averaging is **statistically justified** for 4,689 genes (43.6%) where α ≥ 0.7" | α ≥ 0.7 is a psychometric convention for questionnaire scales. It does not establish that averaging probes yields a more accurate measure of transcript abundance. |
| "The averaged value is **more reliable than any single probe**, median reliability gain **+0.159**" | **The 'gain' was circular.** It was computed as α − r̄, where r̄ was back-solved from α and k via the Spearman–Brown identity. Recomputing it from α and k alone reproduces the reported values to within **8.9 × 10⁻¹⁶**. It contained no information beyond α and k and measured no external accuracy. |
| "Averaging is **not justified** for 3,541 genes (32.9%) with α < 0.5" | Low co-variation among probes is not proof that an average is invalid; it may reflect isoform-specific probes, differing dynamic range, or probe-level noise. |
| "α rises with k **exactly as theory predicts**" | Presented as confirmation. It is largely an artefact — see §4.3. |

### 4.1 What Cronbach's α is, and what it assumes

α is an **internal-consistency** statistic: it summarises how much a set of items co-vary across observations. Applying it to microarray probes assigned to one gene requires assumptions that are **not met here**:

| Assumption | Status for these data |
|---|---|
| (a) items measure **one latent construct** | **Unverifiable.** Probes may target different transcripts, isoforms or 3′UTR regions of a locus. |
| (b) measurement errors are **uncorrelated** | **Violated.** Probes share the array, sample preparation, hybridisation, scanning — and the depositors' joint RMA + ComBat. All induce correlated error, which **inflates α**. |
| (c) items are at least **tau-equivalent** | Not established. |
| (d) a **reference measurement** exists against which accuracy can be assessed | **Absent.** No qPCR or RNA-seq on these samples. No accuracy claim is testable. |

Because of (b) in particular, α here is an **upper-biased** description of probe co-variation, not a reliability coefficient in the psychometric sense.

### 4.2 Observed distribution — descriptive

α across the 10,751 multi-probe genes: median **0.653**, IQR 0.410–0.820, range −1.130 to 0.998.

**Descriptive bands. These are not decision thresholds** — the conventional 0.7 cut comes from psychometrics and has **no established biological meaning** for microarray probe sets.

| Band | Genes | % |
|---|---:|---:|
| α ≤ 0 — probes unrelated or anti-correlated | 360 | 3.3% |
| 0.0–0.5 — low agreement | 3,181 | 29.6% |
| 0.5–0.7 — moderate agreement | 2,521 | 23.4% |
| 0.7–0.8 — high agreement | 1,655 | 15.4% |
| α ≥ 0.8 — very high agreement | 3,034 | 28.2% |

### 4.3 α is mechanically driven by probe count — new evidence

| k probes | genes | median α | median implied mean inter-probe r |
|---:|---:|---:|---:|
| 2 | 5,256 | 0.560 | **0.389** |
| 3 | 2,814 | 0.654 | 0.387 |
| 4 | 1,376 | 0.701 | 0.369 |
| 6 | 327 | 0.770 | 0.358 |
| 8 | 49 | 0.796 | 0.327 |
| 9 | 30 | 0.793 | **0.299** |

**α rises with k while the underlying mean inter-probe correlation falls.** Genes with more probes do not have better-agreeing probes — they have *worse*-agreeing ones. The rise in α is the Spearman–Brown formula, not a measurement.

This is the clearest demonstration that **α cannot serve as an aggregation decision rule**: a gene can cross 0.7 purely by having more probes assigned to it, irrespective of whether those probes agree.

### 4.4 What this section does and does not support

**Does support:** probe sets vary widely in how much they co-vary, from anti-correlated to near-identical. That variation is real and is consistent with the within-gene ρ findings in §1 and Phase 2.

**Does not support:** that averaging is more accurate for any gene; that a single probe is more accurate for any gene; any α threshold as a selection rule.

**A real answer would require an independent reference measurement on these samples — qPCR, RNA-seq, or a spike-in standard. None exists.** The aggregation choice is therefore a design decision made on grounds of uniformity, reportability and transportability, **not** an empirically optimised one.

## 5. Mapping-count reconciliation (8,893 vs 8,894)

**Complete, mutually exclusive accounting of all 54,675 expression probes:**

| Category | n | % |
|---|---:|---:|
| unmapped (0 gene symbols) | **8,893** | 16.27% |
| unique (exactly 1 symbol) | **43,163** | 78.94% |
| ambiguous (2+ symbols) | **2,619** | 4.79% |
| **TOTAL** | **54,675** | 100% |

Probes in no class: **0**. Probes in more than one class: **0**. Annotation probe set equals expression probe set: **True**.

Of the 43,163 unique probes: 11,064 in single-probe genes (11,064 genes) + 32,099 in multi-probe genes (10,751 genes) = 43,163. ✓

**Cause of the discrepancy.** The raw annotation table parses to **54,676 rows**. One is the SOFT control record **`!platform_table_end`**, whose `Gene Symbol` is NaN. Script 11 filters with `dropna(subset=["ID"])`, which **keeps** that row because the string `"!platform_table_end"` is not NaN — so it counts 54,676 "probes" and 8,894 "unmapped". `src/annotation.py` filters IDs beginning with `!`, giving 54,675 and 8,893.

**8,893 is correct. 8,894 counts a file terminator as a probe.** Script 11 is unmodified — this is a reconciliation, not a fix.

---

## 6. Transportability

| | representative_iqr / _mean | median / mean |
|---|---|---|
| **frozen object** | list of 21,815 `(gene → probe_id)` pairs | `gene → {probe_id,…}` groupings over 43,163 probes |
| **needs target-cohort data to apply?** | **No** | No, *if* every probe is present |
| **depends on annotation version?** | only to name genes | **Yes** — groupings change with it |
| **re-derivation in target cohort?** | none | recompute per cohort |
| **probes required** | **21,815** | **43,163** |

**Both can be frozen on discovery.** They differ in what breaks:

- **representative**: needs its 21,815 specific probes present.
- **median/mean**: needs *all* probes of each gene present. If any are missing, the aggregate is computed over a different probe set than in discovery — **the same gene name silently denotes a different quantity.**

### 6.1 Limitations, affecting both

1. **Annotation version.** GPL570 symbols have been revised since deposit. Freezing the gene→probe list freezes the symbol vocabulary with it; a target cohort annotated against a newer build may disagree about which gene a probe belongs to.
2. **Probe availability.** A same-platform GPL570 cohort carries all 54,675 probes, so neither strategy loses probes. A *different* Affymetrix platform carries only a subset. **This must be verified per cohort in Phase 4, not assumed** — I have not verified the overlap for any specific candidate cohort.
3. **Neither strategy addresses the deposited RMA + ComBat scale.** Freezing a probe list does **not** make values comparable across cohorts. That is the Phase 4 problem and is untouched here.

---

## 7. Approved representation

### 7.1 The decision as approved

| Role | Strategy | Selection population | Frozen object |
|---|---|---|---|
| **PRIMARY** | `median` aggregation | **none** — median has no selection step | `gene → {probe_id,…}` groupings, 21,815 genes over 43,163 probes |
| **SENSITIVITY** | `representative_iqr` | **`discovery_443`** | `gene → probe_id`, 21,815 genes, 0 ties |
| historical evidence | `mean` | — | original candidate matrix preserved, **not** a primary strategy |
| **not approved** | α-stratified hybrid | — | **not implemented**; no α threshold exists in any code path |

### 7.2 What is and is not claimed for `median`

**`median` is a provisional, uniform aggregation rule.** It is chosen because it is a single rule applied identically to every gene, it is simple to state in a methods section, it uses all uniquely-mapped probes rather than discarding siblings, and it is robust to one discordant probe in a way `mean` is not.

**It is explicitly NOT claimed to be biologically optimal.** No evidence in this project establishes that any aggregation is closer to true transcript abundance, and — per §4 — no such evidence is obtainable from these data, because no independent reference measurement exists for these samples.

The revision-1 argument that a constructed average is inherently inferior to a representative probe was **withdrawn** (§4.0); so too is the α-based argument that averaging is superior. Neither direction is supported.

### 7.3 Why `representative_iqr` remains the principal sensitivity analysis

It differs from `median` in a structurally important way: it preserves a single measured probe per gene rather than constructing a value, and its frozen object is smaller and annotation-version-robust (§6). If a downstream finding holds under both, the aggregation choice is not driving it. If it does not, that is a reportable result about the finding's fragility — which is exactly what a sensitivity analysis is for.

Its selection is frozen on `discovery_443` and verified free of validation samples (§0.1).

### 7.4 `mean` — retained, not primary

Retained as historical evidence for auditability. It is identical to `median` for the 5,256 two-probe genes (where the median *is* the mean), correlates at r = 0.99 overall, and is strictly less robust to a discordant probe. It adds no information as a distinct strategy.

### 7.5 α-stratified hybrid — not approved

Not implemented. No α threshold appears in any module, script or configuration file. The reasoning is in §4: α is descriptive, is inflated by correlated probe error, rises mechanically with probe count, and carries no validated biological threshold. Using it to switch aggregation rules per gene would encode a psychometric convention as a biological decision.

## 8. Scientific limitations

1. **No reference measurement exists.** No qPCR, RNA-seq or spike-in standard accompanies these samples. **No aggregation strategy can be shown to be more accurate than another on these data.** The choice of `median` is a design decision, not an empirical optimum.
2. **α is inflated and mechanically k-dependent** (§4.1–4.3). It describes probe co-variation under violated assumptions and must not be read as reliability.
3. **Residual batch association says nothing about cross-cohort batch effects** (§2). The data were already ComBat-corrected by the depositors.
4. **Variability is not biological quality.** The IQR criterion cannot distinguish technical variability from biological signal. The `_x_at` evidence (§3) shows the worst failure mode is not dominant; it does not show IQR tracks biology.
5. **Annotation-version dependence.** GPL570 symbols have been revised since deposit. `median` groupings are more exposed to this than a frozen probe list, because changing a symbol changes which probes are averaged.
6. **Probe availability in external cohorts is unverified.** A same-platform GPL570 cohort carries all 54,675 probes; a different Affymetrix platform carries a subset. **This has been verified for no candidate cohort** and is a Phase 4 obligation. `median` requires all 43,163 probes; `representative_iqr` requires 21,815.
7. **The deposited RMA + ComBat scale is untouched.** Freezing a representation does not make values comparable across cohorts. That remains the Phase 4 problem.
8. **11,512 probes (21.1%) are unrepresentable at gene level** — 8,893 unmapped, 2,619 ambiguous. Signal carried only on those probes is invisible to every gene-level matrix. The probe-level cache remains canonical and available.
9. **Selection-set instability is concentrated where it matters least** (§1.1) — but 5.04% of multi-probe genes do change probe between `all_585` and `discovery_443`, and the frozen map is the `discovery_443` one.

---

## 8. Status

| Check | Result |
|---|---|
| Phase 0 suite | 21/21 passed |
| Raw data SHA-256 | both OK |
| Probe-level cache | **unchanged** |
| Original 4 candidates | **preserved**; 2 new ones written separately |
| Scripts 01–14, Phase 1 manifest | unmodified |
| Outcome columns read in Phase 2/2B code | **none** |
| Normalization / batch correction | **none** |
| Samples in every matrix | **all 585** |

**Not started:** survival modelling, feature selection, expression filtering, normalization strategy, external validation, Phase 3.
