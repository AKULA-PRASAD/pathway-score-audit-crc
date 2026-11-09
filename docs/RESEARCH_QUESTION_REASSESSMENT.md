# Research question reassessment

**Date:** 8 November 2026.

**Status:**
- Assessment document. The original registered study is unchanged and remains primary.
- The extension question below was registered separately, before analysis, in `docs/RESEARCH_EXTENSION_REGISTER.md` (SHA-256 recorded in `reports/research_extensions/registration.json`).

## 1. The original question and its answer

**Original (pre-specified, primary):**
> Do four biologically selected tumour pathway-expression scores (EMT, E2F targets, Wnt/β-catenin, p53) add reproducible recurrence-discrimination information beyond clinical variables in GSE39582?

**Answer:** no reproducible improvement was demonstrated.
- 452 patients, 136 recurrences.
- Corrected ΔC **+0.004316**; all five Monte Carlo replicate 95 % intervals include zero (about −0.009 to +0.05).
- Calibration slopes 0.833 (M1) and 0.862 (M2).
- No primary pathway hazard ratio excluded 1.

This answer is final for the original question. Nothing in the extensions revisits it.

## 2. Why another prediction question is not the strongest next step

| Candidate | Requirement | Available | Verdict |
|---|---|---|---|
| Subtype × pathway interactions for recurrence | Adequate events per subtype | 14–35 recurrences per CIT subtype | Underpowered |
| MMR-stratified pathway models | Events in dMMR | 10 recurrences | Underpowered |
| Regularised gene-level Cox model | Many events per effective parameter; nested validation | 136 events | High overfitting risk; low expected value after a null aggregate result |
| Clinical + composition model | — | Done in X7: composition added no discrimination (corrected C 0.601 vs 0.604) | Already answered |
| External predictive validation | ≈ 340 events, compatible endpoint, all predictors | No eligible public cohort (X8) | Not feasible |

## 3. Candidate new questions

| Candidate | Power | Novelty for this project | Uses outcome? | Selected |
|---|---|---|---|---|
| A. "How do subtype and composition influence the recurrence associations of pathway programmes?" (owner's suggestion) | Low (events per stratum) | Moderate | Yes | No, in this form |
| **B. "What do predefined Hallmark pathway scores measure in bulk CRC tissue, how much of their variation reflects genotype, molecular subtype and estimated composition, and is that structure reproducible in independent tumours?"** | **High: n = 452, no events needed; independent cohorts with genotype data** | **High: gene-level and conditional structure not previously examined here; external replication not previously attempted** | **No** | **Yes** |
| C. "Which preprocessing choice changes pathway scores most?" | High | Low (Phase 2B and 8C already showed representation r ≥ 0.99) | No | No (largely answered) |

**Why B was chosen:**
- It is adequately powered.
- It is outcome-free, so it cannot touch the frozen estimands.
- It answers the question every bulk pathway-score study implicitly assumes: that a score measures its named biology.
- It builds on, rather than repeats, X1 (molecular anchoring), X5 (composition estimates) and X6 (EMT–stroma).

Candidate A survives in a narrowed form: subtype and composition are used to *interpret* the scores, not to model recurrence within strata.

## 4. Registered extension question

> What do predefined Hallmark pathway scores measure in bulk colorectal cancer tissue? How much of their variation reflects tumour genotype, molecular subtype and estimated tissue composition, and is that structure reproducible in independent tumours?

Answered by three registered extensions:
- **E1:** within-pathway gene architecture.
- **E2:** composition- and subtype-conditioned structure.
- **E3:** independent outcome-free replication in GSE41258, GSE143985 and GSE103479.

The question was registered **after** the primary results and the X1–X9 results were known. It is therefore a post-hoc extension, not a pre-specified study aim.
