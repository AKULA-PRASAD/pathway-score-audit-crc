# Single-cell extension plan (S1–S5): which cells express the pathway genes?

**Status:** registered plan, written **before** any single-cell data were downloaded or inspected. Awaiting owner approval for the downloads listed in §3. Nothing in this plan changes the primary study, its frozen outputs or its estimand.

**Author:** Prasad Akula · **Date:** 8 November 2026

## 1. Why

The bulk analyses established that the Hallmark EMT score tracks estimated stromal content in four cohorts: gene-disjoint ρ 0.94–0.98, each above its own random-gene-set baseline (`docs/FINAL_RESEARCH_CONTRIBUTION.md`). Bulk data cannot say *which cells* express the EMT genes. Open question 1 of the project is therefore:

> Is the EMT-associated expression in tumour cells or in stroma?

Single-cell RNA-seq measures expression per cell, and cell types are annotated by the original authors. That lets us measure what bulk data could only estimate:
- which compartment contributes the EMT-set transcripts;
- whether malignant epithelial cells carry a coherent EMT-like programme;
- whether a pseudo-bulk EMT score tracks the **counted** fibroblast fraction rather than a computational estimate.

## 2. Questions and decision rules

All analyses are outcome-free: no survival or clinical-outcome data are used. Gene sets are the four frozen Hallmark sets (`config/pathways/m2_panel_4set.json`). The six fibroblast-marker genes shared with MCP-counter are removed for the gene-disjoint variant (EMT-ex6), as in the bulk analysis.

**Compartments.** Each cell takes the original authors' annotation, mapped before analysis to four compartments:
- epithelial (including malignant epithelial);
- stromal (fibroblasts, myofibroblasts, pericytes, smooth muscle);
- immune;
- endothelial.

The mapping table will be committed before any expression values are read. Tumour samples only, unless stated otherwise.

| ID | Question | Statistic | Supported if |
|---|---|---|---|
| **S1** | Do stromal cells, rather than epithelial cells, supply most EMT-set transcripts in tumours? | Per patient: the share of total EMT-set UMIs from each compartment, normalised by the same share for 1,000 random gene sets matched on mean expression (enrichment ratio) | Stromal enrichment > epithelial enrichment in the discovery cohort (paired Wilcoxon, two-sided p < 0.05), and in the same direction in at least one replication cohort |
| **S2** | Does S1 hold for the gene-disjoint EMT-ex6 set? | As S1 | As S1 |
| **S3** | Do malignant epithelial cells carry a coherent EMT-like programme? | Per patient, within epithelial cells only: mean pairwise Spearman correlation of EMT-set genes across cells, versus the 95th percentile of matched random sets. Ambient-RNA control: repeat after excluding genes with fibroblast-specific expression (the top decile of the stromal-to-epithelial ratio in normal mucosa) | "Coherent tumour-cell EMT-like programme detected" if coherence exceeds the random 95th percentile in more than half of patients in the discovery cohort and in at least one replication cohort, with and without the ambient-RNA control. Otherwise "not detected at this resolution" (not "absent") |
| **S4** | In pseudo-bulk, does the EMT score track the counted fibroblast fraction? | Per patient: sum all tumour cells into one pseudo-bulk profile; score it exactly as in the bulk study; Spearman ρ with the fraction of stromal cells counted from annotations; bootstrap 95 % CI | ρ > 0 with a CI excluding 0 in the discovery cohort. Reported alongside the cohort's random-set baseline |
| **S5** | Where do the E2F, Wnt and p53 set transcripts come from? | As S1, for each set | Descriptive; no decision rule |

**Multiplicity.** S1–S4 are four pre-specified tests, reported with unadjusted p-values and Holm-adjusted p-values. S5 is descriptive.

**Interpretation limits, fixed in advance:**
- scRNA-seq measures mRNA, not protein or pathway activity.
- Tissue dissociation under-represents fibroblasts, so compartment shares are biased against stroma. This makes S1 conservative and S4 noisier.
- Author annotations are used as given; re-clustering is a sensitivity analysis only.
- The single-cell patients are not the GSE39582 patients, and stage mixes differ.
- A negative S3 result does not show that tumour-cell EMT is absent.

## 3. Data (public, processed files only)

| Role | Accession | Content | Download |
|---|---|---|---|
| Discovery | GSE132465 (Lee et al., *Nat Genet* 2020; SMC cohort) | 63,689 cells, 23 patients, tumour + 10 matched normal; raw UMI counts + cell annotation | ≈ 130 MB |
| Replication 1 | GSE144735 (Lee et al. 2020; KUL3 cohort) | 27,414 cells, 6 patients, tumour core/border + normal; raw UMI counts + annotation | ≈ 57 MB |
| Replication 2 (optional) | GSE178341 (Pelka et al., *Cell* 2021) | 62 patients, MMRd and MMRp; HDF5 counts + metadata | ≈ 1.2 GB |
| Not used | GSE200997 | Counts and annotations; kept in reserve if a replication cohort fails QC | — |

Raw reads are controlled-access and are **not** requested. Every downloaded file gets a `PROVENANCE.json` with its URL, date and SHA-256, under `data/external/single_cell/` (gitignored, as for the other external data).

**Data QC gate, applied before any test:**
- annotations map to all four compartments;
- each patient retained has at least 50 epithelial and 20 stromal cells;
- positive controls behave as expected: EPCAM, KRT8 and KRT18 epithelial; COL1A1, DCN and LUM stromal; PTPRC immune; PECAM1 and VWF endothelial. Each must be most expressed in its own compartment.

A cohort that fails the gate is excluded and reported, as GSE103479 was in the bulk study.

## 4. Methods and software

- **Software:** Python with Scanpy and AnnData in a separate environment, with pinned versions recorded in `results/single_cell/ENVIRONMENT.txt`. One fixed seed (20261108) for random gene sets and bootstrap.
- **Expression:** counts normalised per cell to 10,000 and log1p-transformed for coherence (S3); raw UMIs for the transcript shares (S1, S2, S5) and pseudo-bulk (S4).
- **Random sets:** 1,000 per test, drawn from genes expressed in at least 1 % of cells, matched to the target set's size and mean-expression decile.
- **Outputs:**
  - `scripts/single_cell/` (one script per step, plus `run_all.sh`);
  - `results/single_cell/` (tables, JSON);
  - `figures/single_cell/`;
  - an independent validation script with fixed checks;
  - a results report that applies the decision rules above verbatim.

## 5. What this adds

- **Science:** a measured answer to the project's open question 1 at single-cell resolution, without using outcome data.
- **Tool (conditional):** only if S1 or S3 is informative, a `pathway-audit sc-sources` command that reports which annotated cell types supply each gene set's transcripts. It would have its own tests, and would be added only after the registered analysis is complete.
- **Skills:** Scanpy and AnnData, single-cell QC, pseudo-bulk analysis and ambient-RNA controls.

## 6. Deviations

Any change to this plan after data access will be recorded here, dated and justified, before the affected results are reported.

- **2026-11-08, after reading annotation files only (no expression values):** the compartment mapping is fixed in `config/single_cell/compartment_map.json`. Two clarifications follow from the annotations. (1) GSE144735 has tumour core, tumour border and normal samples; tumour core ("Tumor") is used, as the plan specifies tumour samples. (2) GSE178341's released labels do not separate malignant from normal epithelium, so for S3 its "malignant epithelial" set is all epithelial cells from tumour specimens, which may include some normal epithelium. Enteric glia and Schwann cells are excluded from the four compartments.
- **2026-11-08, before any expression-based result:** for S3, malignant epithelial cells are capped at 1,000 per patient by a seeded random subsample (seed 20261108), so the 370,000-cell atlas can be processed in 17 GB of memory. This bounds computation only; patients with fewer cells use all of them. S3 ranks use average ranks for ties.
