# Independent audit of the exploratory analyses (X1–X9)

**Date:** 8 November 2026.
**Reviewer:** a separate review process given only the registered plan, the input files and the output files. It was told not to use the analysis scripts: it wrote its own implementations and shared only the frozen project loaders (`load_context`, `score_panel`, `build_clinical`, `rcs_knots`). Its scripts are kept outside the repository, in the session scratch area. It made no repository changes.
**Automated validation:** `scripts/exploratory/x_validate.py` → `reports/exploratory/validation_checks.json`: **22 / 22 PASS** after the corrections below.

## 1. Independent recomputation

| Analysis | Quantity | Reported | Independent | Result |
|---|---|---|---|---|
| X1 | p53 score, TP53 mutant vs WT (n, Δmedian, rank-biserial, p) | 156/130, −0.448, −0.296, 1.7 × 10⁻⁵ | identical | MATCH |
| X1 | Wnt and p53 scores, BRAF mutant vs WT | −0.455 / +0.482 (rank-biserial) | identical | MATCH |
| X2 | Recurrences / deaths without recurrence | 136 / 46 | 136 / 46 | MATCH |
| X2 | Aalen–Johansen recurrence at 60 (all) | first run 0.3075 | 0.3069 | **MISMATCH → corrected** (now 0.3069) |
| X2 | Death without recurrence at 60, stage III | first run 0.0441 | 0.0517 | **MISMATCH → corrected** (now 0.0517) |
| X2 | Naive 1 − KM at 60 | 0.3134 | 0.3134 | MATCH |
| X5 | ESTIMATE stromal and immune, all 452 | — | max \|diff\| 5 × 10⁻¹⁰ | MATCH |
| X5 | MCP-counter fibroblasts, all 452 | — | max \|diff\| 1.8 × 10⁻¹⁵ | MATCH |
| X6 | ρ(EMT, MCP fibroblasts); EMT ∩ ESTIMATE stromal | 0.94981; 26 | 0.94981; 26 | MATCH |
| X7 | R0 EMT / E2F HR; R1 EMT / E2F / stromal HR; VIFs | 1.122 / 0.821; 1.810 / 0.778 / 0.609; 12.8 / 19.0 | identical | MATCH |
| X8 | GSE14333 and GSE103479 stage II–III n / events | 185/48; 156/42 | identical | MATCH |

**Cause of the X2 mismatch.** lifelines' `AalenJohansenFitter` breaks ties between event types by random jitter, and 109 of 134 distinct event times are tied. The estimator was replaced with R `survival::survfit` multistate, which is tie-aware and deterministic. All X2 tables and figure XF2 were regenerated.

## 2. Critique and resolution

| # | Issue raised (severity order) | Resolution |
|---|---|---|
| 1 | The unpenalised E2F HR (0.82, p 0.028) could be misread as overturning the ridge primary estimate | Every unpenalised HR is labelled exploratory and shown with the ridge HR. Added an MMR-adjusted sensitivity analysis: E2F 0.87 (0.72–1.05). Stated explicitly that this is not evidence of prognostic value |
| 2 | EMT vs stromal coefficients in R1/R2 are collinear sign-splitting | Coefficients declared non-interpretable individually. Block likelihood-ratio tests added (stromal + immune block p 0.081). No model promoted |
| 3 | CIT-subtype tests are circular, and the BRAF/CIMP/CIN findings are one axis | BH also reported across the 20 non-CIT tests (only KRAS–E2F moves to q 0.050). The cluster is described as one finding. The possibility of immune genes in the p53 set is stated |
| 4 | ComBat protected a covariate aligned with MMR, so MMR differences may be inflated | Caveat added to X4 and X7 and to the interpretation, citing Nygaard et al. 2016. Correction: in this cohort the MMR-missing patients were split between modcan2 (38) and modcan1 (5), not all assigned to pMMR |
| 5 | Shared measurement in X6; EMT-ex is not independent | Caveats added. The "consistent with, not established" wording is kept |
| 6 | X2 tie handling; death-without-recurrence models have about 4.6 events per term | Estimator corrected. Death models marked unstable and descriptive |
| 7 | Frailty is weakly identified with 9 batches; coarse grid | Grid refined (bounds 0.0275–0.50). Classified PARTIALLY ADDRESSED by the pre-registered rule, never "resolved" |
| 8 | X7 bootstrap leakage | **None found**: knots, gene-level and score standardisation, and composition standardisation are refitted in each resample |
| 9 | ESTIMATE purity calibration | Purity labelled an uncalibrated, relative index |
| 10 | X8: GSE14333 includes rectal tumours and overlaps Moffitt; "exact score reproduction" is about measurement equivalence | Both stated in the audit document |
| 11 | Confirmatory-looking language | Every document and figure caption carries "post-hoc exploratory, same data as primary" |

## 3. Requirement checks

| Requirement | Evidence | Result |
|---|---|---|
| New analyses checked against their input files | Validation checks 2–21; independent recomputation (§1) | PASS |
| Sample and event counts | 452 / 136; 46 deaths without recurrence; batch events sum to 136; TP53 156/130 | PASS |
| Pathway definitions and gene mappings | Frozen scores reproduce Phase 8C to 2.2 × 10⁻¹⁶; overlaps recomputed independently | PASS |
| No patient or outcome leakage | X1, X4, X5, X6 are outcome-free. X7 bootstrap refits every step per resample (independently checked). No tuning on outcome | PASS |
| No discovery/validation leakage | CIT partition used only descriptively (batch table); never as validation | PASS |
| External cohorts independent | Independence assessed per cohort; none used for validation; no pooling | PASS |
| Confidence intervals and tests | Recomputed independently (§1); one estimator error found and corrected | PASS after correction |
| Plots match underlying data | Figures read only the CSV/JSON outputs; values in captions come from those files | PASS |
| Limitations not presented as resolved | Register has no RESOLVED rows; deck slide 16 shows before → after with remaining gaps | PASS |
| Scripts reproducible | `scripts/exploratory/run_all.sh` re-runs X1–X8, figures and validation; fixed seeds | PASS |
| Frozen outputs unchanged | Protected-file and figure hashes re-verified after all work (see final report) | PASS |
| References and data sources | Downloads recorded in `data/external/exploratory/PROVENANCE.json` (URL, release, date, SHA-256, licence) | PASS |

## 4. What this audit cannot certify

- Whether the deposited GEO characteristics are themselves correct. For example, GSE24551 appears to have swapped event and time fields.
- The provenance of `combat_mod`.
- Any claim about cell composition or pathway activity beyond computational estimates.
- Rendering of the updated presentation in PowerPoint or Keynote (P6 remains blocked).
