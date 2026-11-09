# Script map

Scripts are numbered in the order they were run. Each later phase reads only the frozen outputs of earlier phases. Scripts 01–12 are kept byte-identical to the original baseline (their hashes are checked by `tests/baseline_checksums.json`), so they are described here rather than edited.

## Primary study

| Scripts | Stage | What they do |
|---|---|---|
| `01`, `09` | Data download | Download the GSE39582 series from GEO (`01`) and the GPL570 platform annotation (`09`) |
| `02`–`08` | Metadata and expression audits | Sample metadata; recurrence-time records by subset, including zero times (`03`, `04`); the stage II–III cohort (`05`); expression table structure and probe consistency (`06`–`08`) |
| `10`–`12` | Probe annotation audits | Platform annotation table; probe-to-gene mapping; within-gene probe correlation |
| `13`, `14` | Phase 0 | Build and independently verify the probe-level expression cache |
| `15`, `16` | Phase 1 | Investigate the cohort and endpoint; build the cohort manifest |
| `17`–`22` | Phase 2 | Probe-to-gene mapping audit; candidate gene-level matrices; selection sensitivity and probe quality; freeze the approved representation |
| `23`, `24` | Phase 3 | Build and verify the analysis bundle |
| `25`, `26` | Phase 5 | Precision-based feasibility of external validation, and its methodological audit |
| `27`–`29` | Phase 6 | Freeze the M2 pathway specification; missingness audit; freeze the complete-case population (A-19) |
| `30`, `31` | Phase 7 | Preprocessing feasibility; quality control |
| `32`–`35` | Phase 8 | Outcome-blind exploration, PCA and batch structure, frozen pathway scores, consolidated review |
| `36`–`39` | Phase 9 | Protocol planning and audit; primary execution of M1 and M2 (`38`); representation sensitivity S-7 (`39`) |
| `40` | Phase 9, historical | First nested-bootstrap interval driver (superseded) |
| `41`–`48` | Phase 9, corrected campaign | Corrected nested-bootstrap driver (`41`); launcher (`42`); verifier (`43`); freeze tooling (`44`); sizing decision (`45`); final interval report (`46`); θ̂ alignment check (`47`); independent checkpoint-level audit (`48`) |
| `49` | Figures | Results figures generated from frozen outputs |

Core library code lives in `src/` (pathway scoring, standardisation, Phase 9 models) and `src/phase9_corrected/` (the nested-bootstrap engine).

## Registered follow-up work

| Folder | Contents | Run |
|---|---|---|
| `exploratory/` | X1–X8 (X9, the recurrence-time-unit search, was a document search recorded in `reports/exploratory/x9_time_unit/`): molecular anchoring, competing risks (Python and R), batch frailty (R), MMR provenance, ESTIMATE re-implementation and MCP-counter, scores versus composition, composition-adjusted Cox, external-cohort audit; shared helpers, figures and validation | `exploratory/run_all.sh` (validation 22/22) |
| `research_extensions/` | E1 gene architecture and clusterwise stability; E2 composition- and subtype-conditioned structure; E3 replication in GSE41258, GSE143985, GSE103479 and TCGA-COAD, with the positive-control QC gate and challenge tests; figures and validation | `research_extensions/run_all.sh` (validation 28/28) |

The reusable tool is not here: see `../tools/pathway-audit/`.
