# Single-cell extension S1–S5: results

**Date:** 8 November 2026. **Plan:** `docs/SINGLE_CELL_EXTENSION_PLAN.md`, registered and committed before any data access (`86cdc4f`); compartment mapping committed before any expression value was read (`7b78178`). **Outcome-free:** no survival or clinical-outcome data were used.

## Answer in one paragraph

In three public colorectal-cancer single-cell cohorts, the Hallmark EMT gene set's transcripts come predominantly from **stromal cells**, not epithelial cells (S1, supported), and this holds with the six genes shared with fibroblast markers removed (S2, supported). Within malignant epithelial cells, the EMT genes were **not** more coherently co-expressed than expression-matched random genes in any patient (S3: not detected at this resolution). A pseudo-bulk EMT score built from all tumour cells tracks the **counted** stromal-cell fraction (S4, supported), well above its random-gene-set baseline. This turns the bulk finding, that the EMT score tracks *estimated* stromal content, into a measured, cell-resolved result. It does not show that tumour-cell EMT is absent.

## Data and quality control

| Cohort | Role | Tumour patients | Passed the gate | Tumour cells |
|---|---|---|---|---|
| GSE132465 (Lee et al. 2020) | discovery | 23 | 18 | 47,285 |
| GSE144735 (Lee et al. 2020) | replication 1 | 6 | 4 | 8,254 |
| GSE178341 (Pelka et al. 2021) | replication 2 | 62 | 39 | 257,251 |

All three cohorts passed the data-quality gate: all four compartments present, and every positive-control marker most expressed in its own compartment (9/9, 9/9 and 8/8 markers present). Patients were retained if they had at least 50 epithelial and 20 stromal cells. GSE144735 retained only four patients, so it can confirm direction but not significance.

## Results against the registered decision rules

| Question | Result | Decision |
|---|---|---|
| **S1** Do stromal cells, rather than epithelial cells, supply most EMT-set transcripts? | Median enrichment versus matched random sets, stromal vs epithelial: **5.08 vs 0.45** in discovery (paired Wilcoxon p = 7.6 × 10⁻⁶; Holm 2.3 × 10⁻⁵); 2.84 vs 0.19 in GSE144735; 9.80 vs 0.66 in GSE178341 (p = 3.6 × 10⁻¹²) | **Supported** |
| **S2** Same, for the gene-disjoint EMT-ex6 set | 4.15 vs 0.48 in discovery (p = 7.6 × 10⁻⁶); same direction in both replication cohorts | **Supported** |
| **S3** Do malignant epithelial cells carry a coherent EMT-like programme? | EMT-set coherence exceeded the random 95th percentile in **0 of 18**, 0 of 4 and 0 of 39 patients; with the ambient-RNA control in 2 of 18, 0 of 4 and 2 of 39 | **Not detected at this resolution** |
| **S4** Does the pseudo-bulk EMT score track the counted stromal fraction? | Discovery ρ = **0.97** (bootstrap 95 % CI 0.86–1.00; p = 6.4 × 10⁻¹¹; Holm 2.6 × 10⁻¹⁰), against a random-gene-set upper bound of 0.69. GSE178341 ρ = 0.67 (0.46–0.81; random upper bound 0.29). GSE144735 (n = 4) uninformative | **Supported** |
| **S5** Where do the E2F, Wnt and p53 transcripts come from? | No strong compartment concentration (median enrichments 0.74–1.23), apart from a mild endothelial excess for Wnt (1.47–1.73) | Descriptive |

Endothelial cells were also mildly enriched for EMT-set transcripts (median 1.59 in discovery, 3.85 in GSE178341), consistent with EMT-set genes such as collagens and matrix genes being expressed by vessels as well as fibroblasts.

## Interpretation

- **What is now measured, not estimated:** in bulk colorectal tumours, the Hallmark EMT score is driven mainly by stromal-cell expression. The bulk analyses inferred this from computational composition estimates; here it is shown at single-cell resolution with author-annotated cell types, in three cohorts, and with counted cell fractions.
- **What S3 does and does not say:** within malignant epithelial cells the EMT genes did not behave as a coordinated programme; without the ambient control, their coherence was below the median of expression-matched random gene sets in every patient (61 of 61). With the control it moved closer to random (below the random median in 35 of 61) but still rarely exceeded the 95th percentile (4 of 61). That is evidence against a strong, uniform tumour-cell EMT programme in these samples. It is not evidence that tumour-cell EMT never occurs: partial or rare EMT states, transcripts lost to 3′ capture, or changes visible only in protein or spatial data would not be detected.
- **Bias that works against the finding:** tissue dissociation under-represents fibroblasts, so the true stromal contribution is, if anything, larger than measured.

## Limits (fixed in the plan)

- mRNA only: no protein, pathway activity or spatial information.
- Author cell annotations used as given; GSE178341 cannot separate malignant from normal epithelium, so its S3 set may include normal cells.
- The single-cell patients are not the GSE39582 patients, and stage mixes differ.
- S3 uses at most 1,000 malignant epithelial cells per patient (seeded subsample), a computational bound declared before analysis.
- One implementation correction before reporting: S3 was first run on all patients with at least 50 malignant cells, then restricted to the gate-passing patients that the plan specifies. Both runs gave 0 patients above the random 95th percentile.

## Files

- Scripts: `scripts/single_cell/` (`run_all.sh`), implemented with pandas, SciPy and h5py (the plan named Scanpy; see its deviation log); environment in `results/single_cell/ENVIRONMENT.txt`.
- Results: `results/single_cell/summary.json`, `qc_gate.json`, `s1_s2_s5_enrichment.csv`, `s3_coherence.csv`, `s4_pseudobulk.csv`.
- Validation: `results/single_cell/validation_checks.json` (16/16), including a raw re-read of two discovery patients that reproduces the streamed sums exactly.
- Figure: `figures/single_cell/FS1_single_cell_sources.png`.
- Data provenance: `data/external/single_cell/<accession>/PROVENANCE.json` (data not redistributed).
