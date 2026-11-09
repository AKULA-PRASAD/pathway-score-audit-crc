# Skills evidence and role fit

Audience: me, and anyone assessing what this single project does and does not demonstrate for MS-level US roles. Every claimed competency below is linked to a file, figure, result or recorded decision.

## 1. Target roles

| Role family | Typical expectations at entry level | How this project relates |
|---|---|---|
| Bioinformatics Analyst | Processing public omics data, QC, annotation, reproducible pipelines, Python/R, version control, clear reporting | Strong match: probe mapping, QC, provenance, tests, documentation |
| Computational Biology Analyst | Turning biological questions into analyses, gene-set methods, interpreting results in biological context | Strong match: pathway selection, bulk-tissue reasoning, subtype and batch interpretation |
| Bioinformatics Data Scientist | Statistical modelling and machine learning on biological data, validation, leakage control | Good match: penalised survival models, nested resampling, optimism correction; no deep learning or large-scale ML |
| Computational Biology Data Scientist | Modelling biological variation, uncertainty and experimental design | Good match: uncertainty quantification, feasibility analysis, batch sensitivity |
| Biomedical Data Scientist (only where the work involves biological data and ML) | Clinical and molecular data integration, survival outcomes, careful claims | Partial match: clinical-plus-molecular survival modelling; no EHR or imaging data |

These expectations summarise common requirements for these titles. They are not drawn from a systematic survey of current postings.

## 2. Skills-to-evidence matrix

| Competency | Evidence in this project | Where to look |
|---|---|---|
| Python scientific computing | Analysis code, pathway scoring, bootstrap engine, figure generation | `src/`, `scripts/`, `portfolio/src/` |
| Biological data processing | GEO SOFT parsing, verified cache, cohort and endpoint construction without altering source data | `docs/PHASE0_CHANGELOG.md`, `docs/PHASE1_COHORT_ENDPOINT_INVESTIGATION.md`, `docs/PHASE3_DATASET.md` |
| Gene-expression analysis | Probe-to-gene mapping (54,675 → 43,163 → 21,815), two frozen representations, QC with soft flags, PCA and batch investigation | `docs/PHASE2_GENE_REPRESENTATION.md`, `docs/PHASE7C_QC_REPORT.md`, `docs/PHASE8B_PCA_BATCH_INVESTIGATION.md`; Figures 2–3 |
| Pathway-level feature engineering | Outcome-blind Hallmark panel from correlation clusters; fold-internal mean-of-z scoring; coherence analysis | `docs/PHASE6C_PATHWAY_PANEL_REVIEW.md`, `docs/PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md`; Figure 4; diagram D2 |
| Survival analysis | Cause-specific Cox models, competing-risk reasoning, informative-censoring assessment, rank-based measures under an undocumented time unit | `docs/PHASE6A_ENDPOINT_RESOLUTION.md`; Figures 1, 5, S2 |
| Machine learning and regularisation | Ridge-penalised Cox with nested cross-validated λ; scale-equivariant penalty design | `docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md` §3.1.1 |
| Leakage prevention | Training-fold-only standardisation and tuning; fold-isolation tests; scorer re-implementation agreement within 3×10⁻¹⁵ | `docs/PHASE8D_CONSOLIDATED_EDA_REVIEW.md`; tests in `tests/` |
| Model validation | Optimism-corrected C; nested bootstrap 5 × 750 × 100; calibration slopes; feasibility-based decision against weak external validation | `docs/PHASE9_CORRECTED_FINAL_REPORT.md` §6–§10; `docs/PHASE5_VALIDATION_PRECISION_ASSESSMENT.md`; Figures 7–9 |
| Statistical uncertainty | Percentile intervals for corrected ΔC; Monte Carlo precision rule; distinguishing estimands; diagnostic intervals | `reports/phase9/corrected_rerun_report/corrected_interval_results.json`; Figures 8, S1 |
| Reproducible research | Hash-frozen code and configuration, checkpointing, gated verifier, portable independent audit, deterministic figures | `docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md`; `scripts/48_phase9_independent_audit.py`; `portfolio/REPRODUCIBILITY_GUIDE.md` |
| Software testing | 1,097 automated tests; portfolio QC with 151 source-linked checks | `tests/`; `portfolio/PORTFOLIO_QC_REPORT.md` |
| Git and release management | Versioned history, private GitHub Release with a hash-verified checkpoint archive | `docs/LARGE_ARTIFACT_ARCHIVE.md`; `reports/phase9/ARCHIVE_MANIFEST_SHA256.txt` |
| Scientific communication | Final report, research report, citation audit, figures, presentation, interactive portfolio | `portfolio/RESEARCH_REPORT.md`; `docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md`; `portfolio/presentation/` |
| Scientific judgement | Recorded decisions DC1–DC7 (pre-specification, penalty scope, CIT internal, seed adjudication, precision top-up, rejecting misleading results, commissioning audits) | `docs/PORTFOLIO_PRESENTATION_BLUEPRINT.md` §F; `docs/OWNER_DECISION_REGISTER.md` |

## 3. What this project does not establish

- **No R/Bioconductor work.** Expression was used as deposited; the fRMA pilot was prepared but not executed.
- **No RNA-seq.** There is no alignment, quantification or count-based modelling (DESeq2, edgeR, limma-voom).
- **No single-cell, spatial or multi-omic analysis**, and no variant calling or genomics pipelines.
- **No workflow managers** (Nextflow, Snakemake), containers, HPC schedulers or cloud deployment.
- **No deep learning** and no large-scale machine learning.
- **No SQL, databases or production data engineering.**
- **No wet-lab work** and no prospective or clinical study design.
- **No external validation** and no demonstrated predictive or clinical utility.
- **One dataset only:** it shows depth on one cohort, not breadth across data types.

These gaps are best addressed by additional projects, for example a bulk RNA-seq differential-expression pipeline in Nextflow with Bioconductor, or a single-cell analysis that would also address this project's cell-composition question.

## 4. Honest positioning statement

This project demonstrates rigorous, reproducible analysis of public gene-expression and clinical data. It covers biologically grounded feature design, penalised survival modelling, leakage-controlled validation, uncertainty quantification and careful scientific interpretation. It reports a null result, and its strength is the reliability of that result. It is evidence of analytical judgement and depth on one dataset, not of expertise across every bioinformatics technology.
