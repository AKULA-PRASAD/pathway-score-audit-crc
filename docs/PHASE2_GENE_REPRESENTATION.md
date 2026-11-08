# Phase 2 — Gene-Level Expression Representation

**Date:** 2026-11-02 · **Status:** Candidates built and compared.
**Superseded in part** — see the notice at §6 and `PHASE2B_SELECTION_SENSITIVITY.md` for the approved representation.
**Reproduce:**
```
.venv/bin/python scripts/17_audit_probe_gene_mapping.py
.venv/bin/python scripts/18_build_candidate_matrices.py
.venv/bin/python scripts/19_compare_candidate_matrices.py
```

## Provenance — carried forward, not re-derived

GSE39582 was deposited as **`RMA normalisation (Bioconductor R affy package) + ComBat (R sva package)`**, applied jointly across all 585 arrays by the depositors.

**Phase 2 performs no normalization, no batch correction and no expression filtering.** Nothing in this phase has been independently normalized or batch-corrected. The nine recovered ComBat batch assignments and the discovery/validation imbalance remain Phase 4 material and are untouched here.

Verified at runtime: **no Phase 2 module or script references any endpoint column** (`rfs.*`, `os.*`). Aggregation is outcome-blind by construction.

---

## 1. Probe-to-gene mapping

Every probe falls into exactly one of three disjoint classes. **Nothing is silently discarded** — the excluded probes are enumerated by ID in `reports/phase2/probes_unmapped.csv` and `probes_ambiguous.csv`.

| Class | n | % | Handling |
|---|---:|---:|---|
| **unique** (1 symbol) | 43,163 | 78.9% | Used for gene-level aggregation |
| **unmapped** (0 symbols) | 8,893 | 16.3% | **Excluded by construction** — no gene to assign to. Retained in the probe cache. Mean log2 3.91 vs 5.49 for mapped probes, consistent with control/EST probes |
| **ambiguous** (≥2 symbols) | 2,619 | 4.8% | **Excluded by construction** — assigning one gene would be arbitrary. Up to 22 symbols on one probe (e.g. `DDR1 /// MIR4640`). Retained in the probe cache |
| **Total** | **54,675** | 100% | exhaustive, matches array size |

Probe sets reconcile exactly: annotation and expression matrix contain **identical** probe ID sets.

**Genes reachable from uniquely-mapped probes: 21,815**

| | n |
|---|---:|
| single-probe genes | 11,064 |
| multi-probe genes | **10,751** |
| probes in multi-probe genes | 32,099 |
| max probes for one gene | 15 (HFE) |

### 1.1 Reproduction of the earlier audits

| Script 11/12 figure | Reproduced |
|---|---|
| unmapped / unique / ambiguous probes | 8,893 / 43,163 / 2,619 ✓ |
| unique genes | 21,815 ✓ |
| genes with multiple probes | **10,751** ✓ |
| probe pairs evaluated | **41,616** ✓ |
| median within-gene Spearman ρ | **0.331803** ✓ *exact* |

**On the ρ figure.** Script 12 computed it over the **566 primary tumours**; Phase 2's matrix holds all **585** samples. Both are reported: recomputing on the same 566 samples reproduces **0.331803** exactly; the 585-sample value is **0.333532**. The +0.001729 difference is entirely the 19 non-tumoural samples that Phase 2 retains by instruction. **A sample-set difference, not a discrepancy.**

---

## 2. Expression distributions and low-expression probes

Probe-level, all 54,675 probes × 585 samples, values exactly as deposited: **0 non-finite, 0 missing**, range −0.399 to 16.112 (log2).

**Low-expression counts — reported, NOT filtered:**

| Threshold | Probes below | % |
|---|---:|---:|
| mean log2 < 4.0 | 19,021 | 34.8% |
| mean log2 < 5.0 | 29,102 | 53.2% |
| mean log2 < 6.0 | 37,031 | 67.7% |
| sd < 0.2 | 7,927 | 14.5% |
| sd < 0.3 | 21,911 | 40.1% |
| sd < 0.5 | 38,447 | 70.3% |

**No expression filter is applied.** These counts describe what a filter *would* remove. Applying one is a separate decision requiring approval.

### 2.1 Within-gene probe agreement — the Phase 0 finding, confirmed

Median within-gene ρ stratified by the **lower** expression (left) and **lower** SD (right) of each probe pair:

| decile | by mean expr | by SD |
|---|---:|---:|
| 1 (lowest) | 0.109 | 0.086 |
| 5 | 0.318 | 0.395 |
| 10 (highest) | **0.797** | **0.821** |

Restricted to expressed **and** variable pairs (mean ≥ 6, sd ≥ 0.5): 3,619 pairs, **median ρ = 0.820**.

**Low within-gene agreement tracks low expression and low variance — it is not evidence of a defective probe-to-gene mapping.** The GPL570 annotation is sound.

---

## 3. Candidate matrices

Four candidates, each **21,815 genes × 585 samples**, float64, written separately to `data/processed/candidates/` (126.6 MB each). **The probe-level cache is unchanged and remains canonical.**

| Strategy | Criterion | Probes consumed | Ties |
|---|---|---:|---:|
| `representative_iqr` | per gene, probe with largest IQR across samples | 21,815 | **0** |
| `representative_mean` | per gene, probe with largest mean expression | 21,815 | **0** |
| `median` | median across the gene's probes | 43,163 | — |
| `mean` | arithmetic mean across the gene's probes | 43,163 | — |

Ties break deterministically on the lexicographically smallest probe ID; **zero ties occurred**.

### 3.1 QC — all candidates

| strategy | genes | samples | non-finite | missing | min | max | mean | sd |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| representative_iqr | 21,815 | 585 | 0 | 0 | 0.843 | 16.112 | 5.608 | 2.073 |
| representative_mean | 21,815 | 585 | 0 | 0 | 0.843 | 16.112 | 5.946 | 2.272 |
| median | 21,815 | 585 | 0 | 0 | 1.346 | 16.112 | 5.382 | 2.006 |
| mean | 21,815 | 585 | 0 | 0 | 1.346 | 16.112 | 5.381 | 1.954 |

Identical sample IDs and column order across all four; identical gene index; **all 585 samples retained**; sample set verified against the Phase 1 cohort manifest. **No clinical cohort exclusion applied to the expression matrices** — the manifest flags 453 samples as the Phase 1 primary cohort, and all 585 are kept here.

### 3.2 Reproducibility

Every strategy rebuilds **bit-identically** from the probe cache. All four deterministic.

---

## 4. Comparative results

### 4.1 Half the genes cannot differ at all

| | n | % |
|---|---:|---:|
| single-probe genes | 11,064 | 50.7% |
| multi-probe genes | 10,751 | 49.3% |

Single-probe genes are **bit-identical across all four strategies** (verified). **All disagreement is confined to the 10,751 multi-probe genes.**

### 4.2 Pairwise agreement on multi-probe genes

| pair | median r | frac r<0.9 | frac r<0.7 | median \|diff\| | p95 \|diff\| | identical genes |
|---|---:|---:|---:|---:|---:|---:|
| rep_iqr vs rep_mean | **1.0000** | 0.411 | 0.281 | 0.000 | 3.305 | **5,801** |
| rep_iqr vs median | 0.8883 | 0.536 | 0.203 | 0.595 | 2.972 | 118 |
| rep_iqr vs mean | 0.9004 | 0.499 | 0.046 | 0.707 | 2.633 | 0 |
| rep_mean vs median | 0.8320 | 0.673 | 0.310 | 0.833 | 3.385 | 0 |
| rep_mean vs mean | 0.8837 | 0.565 | 0.101 | 0.934 | 2.996 | 0 |
| **median vs mean** | **0.9908** | 0.253 | 0.051 | 0.009 | 0.947 | **5,256** |

Two structural facts explain the extremes:

- **median vs mean are near-identical** because **5,256 genes have exactly two probes**, and for a two-probe gene the median *is* the mean by definition. That is exactly the 5,256 identical genes observed. Mean adds essentially nothing over median as a distinct option.
- **rep_iqr vs rep_mean has median r = 1.0000** because the two criteria pick the **same probe for 5,801 of 10,751** multi-probe genes (54.0%); where they agree the columns are identical. Where they differ (46%), they can diverge sharply — p95 absolute difference 3.3 log2 units.

### 4.3 Disagreement is largest where probes disagree

Strategy agreement (`rep_iqr` vs `mean`) by within-gene probe-agreement decile:

| ρ decile | median within-gene ρ | median strategy agreement |
|---|---:|---:|
| 1 | 0.001 | 0.870 |
| 2 | 0.095 | 0.836 |
| 3 | 0.158 | 0.831 |
| 5 | 0.290 | 0.860 |
| 7 | 0.480 | 0.902 |
| 10 | 0.879 | **0.976** |

| | n genes | median strategy agreement |
|---|---:|---:|
| probe ρ < 0.3 | 4,979 | **0.848** |
| probe ρ ≥ 0.7 | 1,874 | **0.968** |

**The trend is upward but not strictly monotone** — decile 1 (ρ ≈ 0.001) sits at 0.870, above deciles 2–4. I do not have an established explanation for that inversion and am not going to invent one; the overall relationship is nonetheless clear.

**The operative conclusion:** where a gene's probes agree, the strategy choice barely matters (r ≈ 0.97). Where they disagree — **4,979 genes, 46% of multi-probe genes** — the strategy choice decides what the gene-level value means.

### 4.4 Distributional consequence

Gene-level mean expression, multi-probe genes only:

| | rep_iqr | rep_mean | median | mean |
|---|---:|---:|---:|---:|
| median | 6.062 | 6.812 | 5.491 | 5.519 |
| IQR | 4.74–7.26 | 5.26–8.11 | 4.27–6.71 | 4.41–6.61 |

Representative strategies sit higher because selecting one probe avoids averaging a bright probe with dim siblings. `rep_mean` sits highest by construction — it selects *for* brightness.

---

## 5. Limitations

1. **Aggregation discards real measurements.** Representative strategies use 21,815 of 43,163 uniquely-mapped probes — **21,348 probes of real data are not represented** in the gene-level value. Median/mean use all of them but return a value that may correspond to no probe's biology.
2. **11,512 probes (21.1%) cannot be represented at any gene level** — 8,893 unmapped and 2,619 ambiguous. If a real signal lives on an ambiguous probe, no gene-level matrix can see it. The probe-level cache remains available precisely for this reason.
3. **GPL570 annotation is a snapshot.** Probe-to-gene assignments are as deposited in the platform record; gene symbols have been revised since.
4. **No expression filter has been applied**, so all four candidates contain ~35% of probes with mean log2 < 4. A filter would change every downstream result and is an unmade decision.
5. **The deposited values carry RMA + ComBat applied by the depositors.** Nothing here is independently normalized. Cross-cohort portability is unaddressed and remains the Phase 4 problem.
6. **The non-monotone agreement pattern in §4.3 decile 1 is unexplained.**

---

> ## ⚠️ SUPERSEDED — see `PHASE2B_SELECTION_SENSITIVITY.md`
>
> **The recommendation in §6 below was NOT adopted.** Two of its arguments
> were subsequently withdrawn:
>
> 1. The claim that a constructed gene-level average is inherently inferior
>    because it "corresponds to no measured transcript" — refuted in
>    Phase 2B §4.
> 2. The counter-argument built on Cronbach's α in Phase 2B revision 1,
>    which was **also** withdrawn as circular and assumption-violating
>    (Phase 2B §4.0). Neither direction is supported by these data.
>
> **The approved representation is:** `median` aggregation (primary,
> provisional) with `representative_iqr` selected on `discovery_443`
> (principal sensitivity). See Phase 2B §0 and §7.
>
> §6 is retained unaltered as a record of the reasoning at the time.

## 6. Recommendation — `representative_iqr`, with evidence

**I recommend `representative_iqr` as the primary gene-level representation**, with `median` retained as the sensitivity comparator.

**1. Where the choice matters, averaging is the weaker option.** 46% of multi-probe genes have probe ρ < 0.3. Averaging discordant probes produces a value that corresponds to no measured transcript. A representative probe always corresponds to a real measurement.

**2. Median's robustness advantage applies where it is least needed.** Median is more robust than mean to a discordant probe — but the comparison in §4.3 shows that where probes agree, *all* strategies agree (r ≈ 0.97). Median's benefit is concentrated exactly where the choice does not matter.

**3. IQR is the better representative criterion than mean.** `rep_mean` selects the brightest probe, which says nothing about whether it varies between tumours — a saturated, invariant probe carries no information. IQR selects the probe that actually varies, and being the middle 50% it is not driven by outliers. Both are equally outcome-blind.

**4. It transfers.** A representative-probe matrix is fully specified by a **list of 21,815 probe IDs**. That list applies unchanged to any other GPL570 cohort, making the Phase 4/6 transportability work a direct lookup. A median/mean matrix depends on which probes are present and measured in each cohort, so the aggregation itself shifts between cohorts — an avoidable source of cross-cohort drift in a project whose entire question is transportability.

**5. It is deterministic.** Zero ties; bit-identical rebuild verified.

**Honest counter-argument:** representative_iqr ignores 21,348 probes of real data, and for a gene whose probes genuinely measure the same transcript, averaging would reduce measurement noise. That is a real cost. The `median` sensitivity analysis exists to quantify it: if a finding holds under both, the aggregation choice is not driving it; if it does not, that is a reportable result about the finding's fragility.

**`mean` should be dropped as a distinct option** — it is identical to `median` for the 5,256 two-probe genes and r = 0.99 overall, while being strictly less robust. It adds a comparator without adding information.

---

## 7. Status

Four candidates built, compared and retained **separately**. **No strategy selected.**

| Check | Result |
|---|---|
| Phase 0 test suite | 21/21 passed |
| Raw data SHA-256 | both OK |
| Probe-level cache | **unchanged**, 54,675 × 585 |
| Scripts 01–14, Phase 1 manifest | unmodified |
| `config/analysis_config.yaml` | no scientific parameter |
| Outcome columns referenced in Phase 2 code | **none** |
| Normalization / batch correction applied | **none** |
| Samples in expression matrices | **all 585** |

**Not started:** survival modelling, feature selection, expression filtering, normalization strategy, external validation, Phase 3.
