# Independent scientific challenge of research extensions E1–E3

**Date:** 8 November 2026.

**Reviewer:** a separate review process.
- It read the registered plan and the outputs, then wrote independent code. It shared only the frozen project loaders and not the extension scripts.
- Its scripts are kept outside the repository, in the session scratch area. It made no repository changes.

**Challenge tests reproduced in the repository:** `scripts/research_extensions/e_challenge.py` → `reports/research_extensions/challenge_tests.json`. These are declared deviations, run after review.

**Validation:** `reports/research_extensions/validation_checks.json`.

## 1. Independent recomputation: every quantity matched

| Quantity | Reported | Independent |
|---|---|---|
| Within-subtype ρ(EMT, MCP fibroblasts), C1–C6 | 0.911 / 0.915 / 0.913 / 0.889 / 0.943 / 0.847 | identical to 4 decimals |
| p53 score on BRAF, unadjusted → adjusted (n = 355, 35 mutant) | +0.887 → +0.030 SD | identical |
| Wnt score on BRAF, unadjusted → adjusted | −0.834 → −0.648 SD | identical |
| Wnt module, BRAF and dMMR rank-biserial | −0.800 / −0.802 | identical |
| R1 in GSE41258 (91 mutant / 53 wild-type) | −0.221 (−0.428 to −0.006) | −0.221 (−0.422 to −0.029) |
| R6 ρ in GSE41258 / GSE143985 / GSE103479 | 0.954 / 0.943 / 0.088 | identical |
| GSE103479 COL1A1–COL1A2 r | 0.004 | 0.004 |

## 2. Challenges, results and resulting classification

| # | Challenge | Test and result | Effect on interpretation |
|---|---|---|---|
| 1 | **Is GSE103479 usable?** | Positive-control co-expression median r 0.07 (others 0.78–0.84). Same-gene probes are uncorrelated (9 COL1A1 probes, r −0.2 to 0.2), and XIST probes do not separate annotated sex. Some rows labelled ZNF70, ATRX and so on separate sex by 5–8 log2 units, and adjacent probe IDs correlate at r ≈ 0.998 despite naming unrelated genes. **The deposited series-matrix values appear misaligned to their probe IDs** (row-label defect; sample columns align with annotation) | GSE103479 excluded from all gene-level inference. R3 and R4 (BRAF hypotheses), testable only there, become **not testable**. The X8 audit used only GSE103479's sample characteristics, never its expression, so it is unaffected |
| 2 | **Did E1 meet its registered stability rule?** | Median ARI: Wnt 1.00; p53 0.64; E2F 0.64; EMT 0.00 (ARI is degenerate with a single module). Post-hoc clusterwise Jaccard: EMT_M1 0.93, E2F_M1 0.86, Wnt_M1 0.79, p53_M1 0.75 stable; p53_M2 0.59 weak; p53_M3 0.63; E2F_M2 0.48 dissolved | **H-E1 not established (not evaluable as registered).** All module findings are exploratory. The Jaccard grading is a declared deviation and does not override the registered rule |
| 3 | **Is the "rise" in partial correlations real?** | Null: partial ρ(EMT, random 194-gene set given stromal + immune) has median 0.27 (95 % 0.13–0.38, 1,000 draws). Only EMT–p53 (0.55) exceeds the null; EMT–E2F 0.17, EMT–Wnt 0.23, E2F–p53 0.25 and Wnt–p53 0.26 lie inside it. Removing p53_M1 genes leaves EMT–p53 partial ≈ 0.50 | The "relationships strengthen after conditioning" framing is **withdrawn**: it is a generic consequence of conditioning a small residual. A residual EMT–p53 association remains: **exploratory association** |
| 4 | **Is R6 circular?** | 6 of the 8 MCP fibroblast markers are Hallmark EMT genes. Gene-disjoint EMT (192 genes) vs the fibroblast markers: ρ 0.943 (GSE39582), 0.950 (GSE41258), 0.939 (GSE143985). Random 198-gene sets give ρ ≈ 0.13 | R6 survives without shared genes. It shows that the EMT score is dominated by an ECM/collagen co-expression programme that also defines the fibroblast estimate, **in three cohorts**. It does not show an independently measured cell fraction |
| 5 | **Is the Wnt-module BRAF effect just dMMR?** | 27 of 35 BRAF-mutant tumours are dMMR. BRAF coefficient: −1.38 SD alone; −0.65 (−0.90 to −0.40) with dMMR; −0.55 (−0.78 to −0.32) with dMMR + composition. Within pMMR: rank-biserial −0.62 (p = 0.002), but only 8 BRAF-mutant tumours | About half the effect is shared with dMMR; a residual component remains on small numbers. **Exploratory association**; external test not possible (challenge 1) |
| 6 | **Is the p53_M2–MMR association just TP53 status?** | Within TP53-wild-type tumours, dMMR vs pMMR rank-biserial 0.62 (22 vs 106). Adjusted model: dMMR +1.24 SD, TP53 −0.56 SD. TP53 is missing for 30 of 59 dMMR tumours | Not explained by TP53 status in complete cases. **Exploratory association** (module unstable; complete-case analysis) |
| 7 | **Is "EMT tracks stroma within subtypes" informative?** | All six subtypes ρ ≥ 0.85, but CIT C4 is the stromal subtype derived from these arrays, and within-subtype score SDs stay large. Gene-disjoint scores give within-subtype ρ 0.58–0.66 (reviewer) | It rules out the relationship being an artefact *between* subtypes; it is not independent evidence. **Supported by computational estimates** |
| 8 | **Duplicate patients in GSE41258?** | 186 samples from 183 patients. De-duplicated: R1 −0.222 (−0.419 to −0.012); R5 +0.205 (−0.028 to 0.427); R6 0.953 | No change in any verdict |
| 9 | **Leakage, multiplicity, overreach** | No recurrence data used anywhere. E2 adjusted contrasts are estimation only (no multiplicity control) and must not be read as tests. CIMP → EMT adjusted (+0.14) and TP53 → E2F adjusted (+0.31) are borderline single-cohort results, and R2 did not replicate. Composition adjustment may over-adjust | Labelled accordingly |

## 3. Final classification of findings

| Finding | Classification |
|---|---|
| The bulk Hallmark EMT score is dominated by a stromal/ECM co-expression programme (gene-disjoint ρ 0.94–0.98 with fibroblast markers in GSE39582, GSE41258, GSE143985 and TCGA-COAD) | **Supported by independent external evidence** (three cohorts, two platforms) |
| The EMT score reflects fibroblast *content* | **Supported by computational estimates** (no histology) |
| EMT–stroma relationship holds within each CIT subtype | **Supported by computational estimates** |
| p53 target-gene score lower in TP53-mutant tumours | **Directly demonstrated in GSE39582**; **supported by independent external evidence** (GSE41258 −0.22; TCGA-COAD −0.42, robust to composition; GSE143985 null) |
| The p53-score elevation in BRAF-mutant and CIMP+ tumours is explained by dMMR and estimated immune/stromal composition (+0.89 → +0.03 SD) | **Exploratory association** (GSE39582; composition estimated) |
| BRAF/CIMP → lower Wnt score persists after dMMR and composition adjustment | **Exploratory association** (single cohort; external test not possible) |
| TP53 → higher E2F score (R2); MSI → higher p53 score (R5) | **Not established** (not replicated) |
| BRAF → Wnt lower (R3) and BRAF → p53 higher (R4) externally | **Not established**: TCGA effects (−0.42, +0.23) vanish after accounting for hypermutation (−0.11 and +0.12, CIs span 0) |
| p53 set contains separable canonical-target and stromal/immune modules (H-E1) | **Not established** (registered stability not met) |
| Gene modules (Wnt_M1 AXIN2/NKD1/TCF7/NOTCH1/PTCH1; p53_M1 stromal-like; p53_M2 canonical targets) and their associations | **Exploratory association / hypothesis** |
| Score correlations strengthen after composition conditioning | **Not established** (artefact; withdrawn) |
| Residual EMT–p53 association after conditioning | **Exploratory association** |
| GSE103479 deposited matrix has a probe-label defect | **Directly demonstrated** (positive controls; independent diagnosis) |

## 4. Addendum (TCGA-COAD, registered separately)

After this challenge, R1–R6 were tested in TCGA-COAD under `docs/RESEARCH_EXTENSION_ADDENDUM_TCGA.md`:
- 286 primary tumours, all stages.
- Positive-control QC median r 0.94.
- Results: R1 −0.42 (−0.54 to −0.29); R3 −0.42 (−0.58 to −0.24); R4 +0.23 (+0.05 to +0.40); R6 ρ 0.98 (gene-disjoint 0.98); R2 +0.09 (CI spans 0); R5 not testable (5 tumours with MSI annotation).

The classifications in §3 were updated accordingly. The TCGA analysis itself has not been independently re-implemented. The `pathway-audit` tool reproduces its R1 direction through a separate code path (TP53 effect −0.71 SD unadjusted, −0.56 adjusted).

## 5. Independent re-check of the TCGA analysis

A separate code re-check recomputed the TCGA analysis with independently written code, without reading the analysis script. It was run within this project, not by an external reviewer. **Every quantity matched** (n 286; TP53 162/105; BRAF 44/223; R1 −0.415; R3 −0.419; R4 +0.233; R6 0.981; gene-disjoint 0.979; control median 0.941). The reviewer's challenges are reproduced in `reports/research_extensions/tcga_challenge_tests.json`:

| Challenge | Result | Classification |
|---|---|---|
| Is R3 / R4 driven by hypermutation (MSI proxy; MSI annotation present for only 6 of 286 tumours)? | 31 of 44 BRAF-mutant tumours have > 300 mutated genes. Excluding them: Wnt −0.11 (−0.44 to 0.21), p53 +0.12 (−0.27 to 0.49). Adjusting for log mutation burden: −0.17 and +0.10, CIs span 0 | **Confounded**; R3, R4 not established |
| Is R1 robust? | Excluding hypermutated tumours −0.37 (−0.53 to −0.20); adjusted for burden −0.64 SD (−0.88 to −0.39) | **Replicated, robust** |
| R2 | Null under the registered rule; +0.20 (0.04 to 0.37) after excluding hypermutated tumours | Not replicated (registered); post-hoc sensitivity noted |
| How large is R6 relative to chance? | Random 189-gene sets correlate with the fibroblast estimate at median ρ 0.64 (95 % 0.46–0.76) in TCGA; −0.14 (−0.35 to 0.05) in GSE39582; −0.50 (−0.72 to −0.15) in GSE41258; 0.00 (−0.26 to 0.26) in GSE143985 | **Replicated**: EMT (0.94–0.98) exceeds every cohort's baseline, but in TCGA much of a correlation of this kind reflects a sample-wide factor |
| Stage restriction | R6 within stage II–III TCGA tumours (n 190): 0.98 | Unchanged |

The earlier statement "random gene sets ≈ 0.13" applied to one cohort and one null construction. It has been replaced everywhere by the cohort-specific baselines above.
