# Documentation map

This folder holds the project's full evidence trail. Most readers need only the first five documents. The rest record how each decision was made, and are here for auditors and for anyone checking a specific number.

## Start here

| Read | For |
|---|---|
| [`FINAL_RESEARCH_CONTRIBUTION.md`](FINAL_RESEARCH_CONTRIBUTION.md) | The whole project on one page: what was established, what was not, and why |
| [`manuscript/MANUSCRIPT_DRAFT.md`](manuscript/MANUSCRIPT_DRAFT.md) | The work written as a short paper (draft, not peer reviewed) |
| [`PHASE9_CORRECTED_FINAL_REPORT.md`](PHASE9_CORRECTED_FINAL_REPORT.md) | The source of record for the primary prediction result |
| [`SCIENTIFIC_RESULTS_AND_FIGURES.md`](SCIENTIFIC_RESULTS_AND_FIGURES.md) | The primary results with the nine main and two supplementary figures |
| [`SCIENTIFIC_LIMITATION_RESOLUTION_REGISTER.md`](SCIENTIFIC_LIMITATION_RESOLUTION_REGISTER.md) | The status of every limitation: addressed, partly addressed or open |

The reusable tool is documented separately in [`../tools/pathway-audit/README.md`](../tools/pathway-audit/README.md).

## 1. Primary study, by phase (2–8 November 2026)

| Phase | Documents | What was decided or done |
|---|---|---|
| 0 | `PHASE0_CHANGELOG.md` | Reproducible infrastructure |
| 1 | `PHASE1_COHORT_ENDPOINT_INVESTIGATION.md` | Cohort and recurrence endpoint |
| 2 | `PHASE2_GENE_REPRESENTATION.md`, `PHASE2B_SELECTION_SENSITIVITY.md` | Probe-to-gene representation (median of uniquely mapped probes) |
| 3 | `PHASE3_PLAN.md`, `PHASE3_DATASET.md` | Analysis dataset, built and verified |
| 4 | `PHASE4_PLAN.md`, `PHASE4_PREACQUISITION.md`, `PHASE4_FRMA_PILOT.md` | External-validation feasibility; fRMA pilot prepared, not run |
| 5 | `PHASE5_VALIDATION_PRECISION_ASSESSMENT.md` | Precision-based feasibility of validation |
| 6 | `PHASE6A_ENDPOINT_RESOLUTION.md`, `PHASE6B_MODELING_PROTOCOL_DRAFT.md`, `PHASE6C_PATHWAY_PANEL_REVIEW.md`, `PHASE6D_MISSINGNESS_AUDIT.md` | Endpoint details, protocol structure, the frozen four-set panel, missing data |
| 7 | `PHASE7A_PREPROCESSING_FEASIBILITY.md`, `PHASE7B_QC_PROTOCOL.md`, `PHASE7C_QC_REPORT.md` | Preprocessing and quality control |
| 8 | `PHASE8A_EXPLORATORY_ANALYSIS.md`, `PHASE8B_PCA_BATCH_INVESTIGATION.md`, `PHASE8C_FROZEN_PATHWAY_SCORE_ANALYSIS.md`, `PHASE8D_CONSOLIDATED_EDA_REVIEW.md` | Outcome-blind exploration and batch structure |
| 9 | `PHASE9_MODELING_PROTOCOL_DRAFT.md`, `PHASE9_FINAL_AUDIT.md`, `PHASE9_EXECUTION_REPORT.md`, `PHASE9_S7_SENSITIVITY_REPORT.md` | Modelling protocol, its audit, and the first execution |
| 9 (interval) | `PHASE9_DELTAC_INTERVAL_AUDIT.md`, `PHASE9_DELTAC_INTERVAL_PLAN.md`, `PHASE9_DECISION13_FREEZE_PROCEDURE.md`, `PHASE9_THETA_ALIGNMENT_AUDIT.md` | The nested-bootstrap interval: method, plan, freeze and alignment |
| 9 (report) | `PHASE9_CORRECTED_FINAL_REPORT.md`, `PHASE9_SCIENTIFIC_INTERPRETATION_REVIEW.md`, `SCIENTIFIC_RESULTS_AND_FIGURES.md`, `SCIENTIFIC_RESULTS_FIGURE_AUDIT.md`, `BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md` | Final report, interpretation review, figures and citation checks |

## 2. Follow-up work (registered before running, outcome-free)

| Strand | Documents |
|---|---|
| Limitation investigation X1–X9 | `EXPLORATORY_ANALYSIS_PLAN.md` (registered plan) → `ADDITIONAL_ANALYSES_METHODS_AND_RESULTS.md` (results) → `EXPLORATORY_INDEPENDENT_AUDIT_REPORT.md` (review); `EXTERNAL_COHORT_ELIGIBILITY_AUDIT.md` (X8) |
| Interpretation | `SCIENTIFIC_INTERPRETATION_AND_CONCLUSION.md`, `SCIENTIFIC_LIMITATION_RESOLUTION_REGISTER.md`, `RESEARCH_QUESTION_REASSESSMENT.md` |
| Research extensions E1–E3 | `BIOLOGICAL_CONTRIBUTION_STRENGTHENING_PLAN.md` → `RESEARCH_EXTENSION_REGISTER.md` (registered) and `RESEARCH_EXTENSION_ADDENDUM_TCGA.md` (registered) → `INDEPENDENT_SCIENTIFIC_CHALLENGE.md` (challenge) → `FINAL_RESEARCH_CONTRIBUTION.md` (synthesis) |
| Single-cell extension S1–S5 | `SINGLE_CELL_EXTENSION_PLAN.md` (registered before data access) → `SINGLE_CELL_RESULTS.md` (results) |

## 3. Project records and release

| Documents | Purpose |
|---|---|
| `OWNER_DECISION_REGISTER.md`, `release/G1_G13_IMPLEMENTATION.md` | Decisions made by the project owner, and the wording fixes they approved |
| `PROJECT_COMPLETION_INVENTORY.md`, `FINAL_PROJECT_HANDOVER.md`, `PROJECT_RETROSPECTIVE_AND_ACHIEVEMENTS.md` | Inventory, handover and retrospective of the primary study |
| `release/RELEASE_READINESS.md`, `PUBLIC_RELEASE_REPRODUCIBILITY_AUDIT.md`, `LARGE_ARTIFACT_ARCHIVE.md` | Release checklist, public-release audit, and how the >100 MB checkpoints are archived outside Git |
| `PORTFOLIO_PRESENTATION_BLUEPRINT.md` | Plan behind the slide deck in `../portfolio/` |

## Glossary of internal labels

| Label | Meaning |
|---|---|
| M1 / M2 | Clinical Cox model / clinical model plus the four pathway scores (ridge-penalised) |
| ΔC, θ̂ | Gain in Harrell's C from M1 to M2; θ̂ is the pre-specified optimism-corrected point estimate (+0.004316) |
| B_out × B_in | Outer × inner bootstrap resamples of the nested bootstrap (750 × 100) |
| A-19 | The 452-patient complete-case analysis population (136 recurrences) |
| S-7 | Sensitivity analysis using a representative-probe gene representation on the same arrays |
| CIT partition | The original study's internal validation subset; shares processing with the discovery set, so it is internal only |
| LOBO | Leave-one-processing-batch-out refits |
| `combat_mod` | The covariate protected during the original ComBat batch correction; it matches recorded MMR status |
| dMMR / pMMR | Mismatch-repair deficient / proficient |
| DC1–DC7 | Portfolio labels for seven key methodological decisions of the primary study (listed in `PORTFOLIO_PRESENTATION_BLUEPRINT.md`) |
| D1–D5, G1–G13 | Owner closure decisions and the approved wording fixes |
| X1–X9 | Registered post-hoc analyses of each limitation |
| E1–E3, R1–R6 | Research extensions, and the six pre-registered replication hypotheses within E3 |
| C1–C6 | Post-review challenge tests (for example C2: gene-disjoint EMT–fibroblast test; C5: TCGA hypermutation; C6: cohort-specific random-gene-set baselines) |
| S1–S5 | Registered single-cell questions (see `SINGLE_CELL_RESULTS.md`) |
