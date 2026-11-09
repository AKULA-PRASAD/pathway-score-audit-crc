# Changelog

How the project developed, using the dates recorded in the commits and in each document. Work was developed in the working tree and committed in logical groups after the owner reviewed it, which is why several commits share a timestamp.

## 2026-11-09: public release
- Repository made public: https://github.com/AKULA-PRASAD/pathway-score-audit-crc.
- Website published on GitHub Pages (https://akula-prasad.github.io/pathway-score-audit-crc/), built from `portfolio/web/index.html` by `portfolio/src/build_site.py` and `.github/workflows/pages.yml`; footer links to the code and the slides.

## 2026-11-08 (evening): single-cell extension and pathway-audit 0.2.0
- Branch pushed to the private GitHub repository, renamed `pathway-score-audit-crc`; GitHub Actions CI passes on Python 3.10 and 3.12.
- Single-cell extension S1–S5 run in GSE132465, GSE144735 and GSE178341 after the plan and the compartment mapping were committed (`docs/SINGLE_CELL_RESULTS.md`).
- `pathway-audit` 0.2.0: single-cell source check (`sc-sources`), 13 tests; reproduces the registered single-cell result.
- Slides 16, 18 and 19, README, interview guide and resume text updated.
- Web page updated for the single-cell results and `pathway-audit` 0.2.0 (interaction checks 76/76) and republished; branch merged into the default branch.
- `pathway-audit` Docker image built and verified locally (Colima, Docker Engine 29): tests (12 passed, 1 skipped) and demo pass inside the container; documents and slide 16 updated. Not published to a registry.

## 2026-11-08: follow-up research, the tool and the portfolio
- Registered post-hoc analyses X1–X9 investigating every limitation (`docs/EXPLORATORY_ANALYSIS_PLAN.md`).
- Registered research extensions E1–E3, with replication in GSE41258, GSE143985 and TCGA-COAD, a positive-control QC gate, and post-review challenge tests (`docs/RESEARCH_EXTENSION_REGISTER.md`).
- `pathway-audit` 0.1.0: a Python package and CLI with tests, a demo and a five-cohort case study (`tools/pathway-audit/`).
- Portfolio: presentation, web page, resume text and portfolio QC (`portfolio/`).
- MIT licence chosen; work committed in groups (`fe7f88f`…`28a9271`).
- Repository usability: data-dependent tests skip cleanly in a fresh clone; documentation and script maps; interview guide updated.
- Single-cell extension plan S1–S5 registered before any data access (`docs/SINGLE_CELL_EXTENSION_PLAN.md`).

## 2026-11-08 (early): primary study closed
- Phase 9 freeze approved; Phases 1–9 source, outputs and evidence committed in groups (`7e0df40`…`97ed2a0`).
- Final scientific report, independent checkpoint-level audit, results figures and retrospective.

## 2026-11-06 to 11-07: corrected nested bootstrap
- Corrected nested-bootstrap pipeline with independent random streams per draw, and the Decision 13 freeze package (`2f5f26c`).
- Path C sizing-decision writer for the pre-specified 500 → 750 top-up (`c432eee`).

## 2026-11-03: Phase 9 protocol and first execution
- Modelling protocol v3.1 and its audit; primary execution of M1 and M2, and the S-7 sensitivity analysis.
- Methodological audit and plan (v2.1) for a confidence interval on the optimism-corrected ΔC.

## 2026-11-02: infrastructure and Phases 1–8
- Phase 0 reproducible infrastructure (`1cd2bba`).
- Cohort and endpoint, gene representation, analysis dataset, validation feasibility, frozen pathway panel, missing data, preprocessing and QC, and outcome-blind exploration (Phases 1–8; documents dated 2–3 November, committed on 8 November).
