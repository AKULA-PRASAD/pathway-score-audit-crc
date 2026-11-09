# Resume and LinkedIn text (evidence-backed)

Every number below is traceable to the evidence table at the end. Pick the 3–4 bullets that best fit the role, and keep the first bullet for the tool.

## Project title line
**Pathway-score audit of colorectal cancer transcriptomes (`pathway-audit`)** · Python, R, survival analysis, GEO and TCGA

## Resume bullets (one or two lines each)

- **Built `pathway-audit`, a tested Python package and CLI** that checks whether a gene-set score in bulk or single-cell tumour data reflects tissue composition, a few driver genes, genotype confounding or a broken dataset; HTML reports, built-in demo, 13 tests, CI on GitHub Actions.
- **Replicated across 4 cohorts and 2 platforms (GEO microarray, TCGA RNA-seq)** that the Hallmark EMT score tracks stromal fibroblast content in bulk colorectal tumours (gene-disjoint ρ 0.94–0.98, each above a cohort-specific random baseline).
- **Confirmed the finding at single-cell resolution in 3 public cohorts (up to 370,000 cells, streamed in blocks):** EMT-set transcripts came mainly from stromal cells (5.1× the expression-matched share vs 0.45× in epithelium), under a plan registered before data access.
- **Caught a defective public dataset (GSE103479)** with a positive-control co-expression gate (median r 0.07 vs 0.78–0.94), preventing a false non-replication.
- **Separated real genotype effects from confounders:** the TP53 → p53-score link replicated in 2 external cohorts and survived hypermutation adjustment, while the BRAF–p53 link was explained by MMR status and composition, and TCGA BRAF effects by hypermutation.
- **Designed a pre-specified nested-bootstrap evaluation (750 × 100 resamples) in 452 colon cancers** showing four pathway signatures add no reproducible recurrence discrimination over clinical factors (corrected ΔC +0.004; all intervals include 0).
- **Ran the work as an audited, reproducible pipeline:** analysis plans registered before execution, hash-frozen code, a 60/60 verifier gate, independent code re-checks and 1,097 automated tests.

**Skills line:** Python (pandas, NumPy, SciPy, scikit-learn, lifelines, h5py), R (survival), single-cell RNA-seq analysis, Cox models and penalised regression, bootstrap validation, GEO/TCGA data, ESTIMATE/MCP-counter, packaging and CLI design, pytest, GitHub Actions CI, Git.

## LinkedIn "Projects" entry (about 600 characters)
I built pathway-audit, a tested Python tool that checks what a cancer pathway score really measures before it is interpreted. Applying it to colorectal cancer, I replicated across four cohorts, including TCGA, that the widely used Hallmark EMT score tracks tumour stroma in bulk tissue, confirmed in three single-cell cohorts that its transcripts come mainly from stromal cells, caught a corrupted public dataset, and separated real genotype effects from confounders. A pre-specified, optimism-corrected test in 452 patients showed that four pathway signatures add no reproducible recurrence prediction beyond clinical factors.

## 30-second spoken version
"I tested whether four well-known cancer pathway scores improve recurrence prediction in colon cancer, using a pre-specified, optimism-corrected bootstrap. They didn't. So I asked what these scores actually measure in a bulk tumour sample. Across four cohorts, the EMT score turned out to track stromal content, and the p53 score tracked TP53 mutation. I packaged those checks into a tested Python tool, pathway-audit, so other analysts can run them before interpreting a pathway score."

## What not to claim
- Clinical utility.
- Pathway activity or protein function.
- Tumour-cell EMT, or that the EMT score "is not" tumour-cell EMT: bulk data cannot separate the two.
- External validation of the prediction model.
- A PyPI package or published container image: neither exists. The repository is public (https://github.com/AKULA-PRASAD/pathway-score-audit-crc), MIT-licensed and CI passes; the Docker image builds and runs locally (tests and demo pass inside it) but is not on a registry.
- Scanpy-based analysis: the single-cell work used pandas, SciPy and h5py directly.
- "Novel algorithm": the tool combines established methods in a new, tested workflow.

## Evidence table (for interviews)

| Claim | Source |
|---|---|
| 13/13 tests; demo recovers 4/4 planted truths; CI passing | `tools/pathway-audit/tests/test_core.py` (with `PATHWAY_AUDIT_REFS` set; 12 pass + 1 skipped without references); GitHub Actions run history |
| Single-cell 5.1× vs 0.45×; 3 cohorts; 370,000 cells | `results/single_cell/summary.json`; `docs/SINGLE_CELL_RESULTS.md` |
| ρ 0.94–0.98 in 4 cohorts; random baselines | `reports/research_extensions/challenge_tests.json`, `e3_tcga_summary.json`, `tcga_challenge_tests.json` |
| GSE103479 median r 0.07 | `results/research_extensions/e3_positive_control_qc.csv` |
| TP53 → p53 score replication; hypermutation | `results/research_extensions/e3_replication_tests.csv`, `e3_tcga_replication_tests.csv`; `tcga_challenge_tests.json` (C5) |
| BRAF–p53 +0.89 → +0.03 SD after adjustment | `results/research_extensions/e2_genotype_adjusted.csv` |
| ΔC +0.004316; 750 × 100; five intervals | `reports/phase9/corrected_rerun_report/corrected_interval_results.json`; `docs/PHASE9_CORRECTED_FINAL_REPORT.md` |
| Verifier 60/60; 1,097 tests | `docs/PHASE9_CORRECTED_FINAL_REPORT.md` |
