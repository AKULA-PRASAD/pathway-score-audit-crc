# pathway-audit: methods and decision rules

All thresholds are fixed in code (`src/pathway_audit/`) and listed here. They are conventions chosen in advance, not values tuned to any dataset.

| Component | Method | Rule |
|---|---|---|
| Score | Mean of per-gene z-scores across samples (ddof = 1), then standardised | Scored only if ≥ 80 % of set genes are present and ≥ 10 genes |
| QC gate | Pearson r of 8 canonical co-expressed pairs: COL1A1–COL1A2, COL1A1–COL3A1, COL1A2–COL3A1, SPARC–COL5A1, MKI67–TOP2A, BUB1–CCNB1, CD2–CD3E, CDX2–VIL1 | PASS if median r ≥ 0.5 with ≥ 4 pairs present |
| Coherence | Mean pairwise inter-gene r | "incoherent" if < 0.05 |
| Gene drivers | Item–rest r (gene vs score without it); leave-one-gene-out r with the full score | "driven by few genes" if min leave-one-out r < 0.95 |
| Modules | Average-linkage clustering on 1 − r, cut at 0.7 (average r ≥ 0.3), modules ≥ 10 genes (≥ 5 for sets < 60 genes) | Clusterwise Jaccard over patient bootstraps (Hennig 2007): ≥ 0.75 stable, ≥ 0.6 pattern, ≥ 0.5 weak, otherwise dissolved |
| Composition | ESTIMATE stromal/immune (re-implementation validated to 1e-14 against the package example); MCP-counter marker means (gene mode) | max \|ρ\| ≥ 0.7 composition-dominated; ≥ 0.4 associated; otherwise largely independent |
| Gene-disjoint test | ρ between the score without fibroblast markers and the fibroblast-marker mean | Required (≥ 0.7) for "holds without shared genes" |
| Genotype effects | OLS difference in score SD, unadjusted and adjusted for stromal + immune (+ user covariates), percentile bootstrap CI | robust: adjusted CI excludes 0 with the same sign. Explained: unadjusted CI excludes 0, adjusted includes 0, and \|adjusted\| < 0.5 \|unadjusted\|. Groups need ≥ 10 samples |
| Score relationships | Partial Spearman given stromal + immune ranks; null from random gene sets of median set size | "exceeds null" outside the null's central 95 % |
| Survival screen | Cox models with and without scores; single-level paired Harrell bootstrap (event-stratified, scores re-standardised in each resample) | Corrected ΔC only; no interval |

## Validation

- **Unit tests on synthetic data with planted structure** (`tests/test_core.py`).
- **ESTIMATE re-implementation** checked against the package's own example (`validate_estimate`).
- **Real-data case study** on five cohorts (`examples/case_study_crc.py`). Its results agree with the independent analyses in the parent project:
  - the gene-disjoint EMT–fibroblast ρ (0.94 in GSE39582) matches the independent challenge;
  - the BRAF → p53 "explained" verdict matches E2;
  - the survival screen (+0.0048) is close to the frozen nested estimate (+0.0043).

## Single-cell source check (`sc-sources`)

**Input:** single-cell UMI counts (cells × genes), a cell-type label per cell, and optionally a group label (for example patient).

**Method:**
1. Keep genes detected in at least 1 % of cells; assign each gene to a decile of total expression.
2. For each group and cell type with at least 20 cells, sum the counts per gene.
3. For a gene set, compute its share of transcripts in each cell type: the set's summed counts in that cell type divided by its summed counts over all cell types.
4. Compute the same share for 1,000 random gene sets matched to the target on size and expression decile (target genes excluded), and average.
5. Enrichment = observed share / random share. With groups, report the median across groups.

**Reading it:** enrichment ≥ 2 means a cell type supplies at least twice the share that expression-matched random genes would; ≤ 0.5 means it is depleted. Values near 1 mean no preference.

**Validated against:** the registered single-cell analysis of GSE132465 (EMT set: stromal 5.10 vs 5.08, epithelial 0.45 vs 0.45; `examples/case_study_single_cell.py`).
