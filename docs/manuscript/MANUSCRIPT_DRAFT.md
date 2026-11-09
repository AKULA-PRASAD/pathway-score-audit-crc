# Pathway-expression scores in bulk colorectal cancer: no added recurrence discrimination, and a stromal reading of the Hallmark EMT signature across four cohorts

**Prasad Akula** · MS in Bioinformatics, Northeastern University

*Draft manuscript, 8 November 2026. Not peer reviewed. Every number is traceable to the files cited in square brackets.*

## Abstract

**Background.** Gene-set scores such as the Hallmark EMT, E2F, Wnt/β-catenin and p53 signatures are widely used as biological readouts and prognostic features in bulk tumour transcriptomes. Whether they add recurrence information beyond clinical variables, and what they measure in bulk tissue, are rarely tested together.

**Methods.**
- In 452 stage II–III colon cancers (GSE39582; 136 recurrences), four outcome-blind Hallmark scores were added to a clinical Cox model.
- The optimism-corrected gain in Harrell's C was estimated with a paired nested bootstrap (750 × 100; five Monte Carlo replicates; pre-specified precision rule).
- Registered, outcome-free follow-up analyses examined what the scores measure: gene-level structure, composition-conditioned genotype effects using ESTIMATE and MCP-counter, and pre-registered replication in GSE41258, GSE143985 and TCGA-COAD.
- A positive-control co-expression gate screened datasets.

**Results.**
- The corrected gain was +0.0043; all five 95 % intervals included zero (about −0.009 to +0.05).
- The EMT score correlated 0.94–0.98 with fibroblast markers after removing shared genes, in all four usable cohorts and within every molecular subtype.
- The p53 target score was lower in TP53-mutant tumours (replicated in two of three external cohorts, robust to hypermutation).
- Its elevation in BRAF-mutant tumours was explained by mismatch-repair status and estimated composition. Lower Wnt-target expression in those tumours was not, but its external test in TCGA was confounded by hypermutation.
- One public dataset (GSE103479) failed the quality gate. Its expression values appear misaligned to probe identifiers.

**Conclusions.** The four pathway scores did not demonstrate added recurrence discrimination. In bulk colorectal tissue, the Hallmark EMT score is predominantly a stromal programme. Pathway scores should be audited for composition, genotype confounding and data quality before biological interpretation. We provide a tool for doing so.

## 1. Introduction

Colorectal cancers are driven by recurrent alterations of the Wnt, p53 and RAS/RAF pathways (TCGA 2012) and stratified into consensus molecular subtypes (Guinney et al. 2015). Bulk transcriptome signatures of these processes are attractive biomarkers, but two problems are underappreciated:
- **Optimism.** Small cohorts make any added feature look predictive.
- **Bulk composition.** The mesenchymal CMS4 signal in bulk tumours arises largely from stroma (Calon et al. 2015; Isella et al. 2015).

We asked a pre-specified prediction question and then, separately, what the scores measure.

## 2. Methods

**Cohort.** GSE39582 (Marisa et al. 2013): 585 arrays, of which 452 complete-case stage II–III patients with 136 recurrences.
- Endpoint: cause-specific hazard of recurrence. The time unit is undocumented (provisionally months); 46 deaths without recurrence were censored [docs/PHASE9_CORRECTED_FINAL_REPORT.md].

**Features.**
- Probe mapping: 54,675 → 21,815 genes (median per gene).
- Four Hallmark sets (EMT 198, E2F 195, Wnt 41, p53 194 genes), chosen from correlation clusters without outcome data.
- Scores: mean of gene z-scores, standardised inside each training fold.

**Models.**
- M1: unpenalised Cox model with age (spline), sex, stage, location and adjuvant chemotherapy.
- M2: M1 plus ridge-penalised scores (λ by inner cross-validation).

**Estimand.** Corrected ΔC, estimated with a paired nested bootstrap [docs/PHASE9_DELTAC_INTERVAL_PLAN.md]. Verification: a frozen verifier (60/60) and an independent recomputation of 138 quantities.

**Follow-up analyses** (registered before execution; outcome-free) [docs/RESEARCH_EXTENSION_REGISTER.md; docs/RESEARCH_EXTENSION_ADDENDUM_TCGA.md]:
- gene-level item–rest correlations and modules with bootstrap stability;
- composition estimates (ESTIMATE re-implemented and validated to 1.4 × 10⁻¹⁴; MCP-counter);
- genotype effects adjusted for composition and mismatch repair;
- partial correlations against random-gene-set nulls;
- gene-disjoint marker tests;
- six pre-registered replication hypotheses with a fixed decision rule, in GSE41258 (MSKCC, U133A), GSE143985 (U133 Plus 2.0) and TCGA-COAD (RNA-seq).

Cohorts were never pooled. A post-hoc positive-control co-expression gate (median r of eight canonical pairs ≥ 0.5) was applied after an implausible result.

**Tool.** The audit is implemented as `pathway-audit` (Python; command-line interface; unit tests on synthetic data with planted structure) [tools/pathway-audit/].

## 3. Results

**3.1 No demonstrated added discrimination** [reports/phase9/corrected_rerun_report/corrected_interval_results.json]

| Quantity | Value |
|---|---|
| Apparent C, M1 vs M2 | 0.6242 vs 0.6382 |
| Corrected C, M1 vs M2 | 0.6016 vs 0.6059 |
| Corrected ΔC | +0.004316; all five Monte Carlo replicate intervals include zero |
| Calibration slopes | 0.833 (M1), 0.862 (M2) |

No primary pathway hazard ratio excluded 1. Post-hoc, composition estimates added no discrimination, and competing-risk (Fine–Gray) estimates matched the cause-specific ones.

**3.2 The EMT score is a stromal programme in bulk tissue** [results/research_extensions/; reports/research_extensions/challenge_tests.json; e3_tcga_summary.json]

| Cohort | n | EMT vs fibroblast markers, gene-disjoint ρ (95 % CI) |
|---|---|---|
| GSE39582 | 452 | 0.94 (0.93–0.96) |
| GSE41258 | 186 | 0.95 (0.92–0.97) |
| GSE143985 | 91 | 0.94 (0.90–0.96) |
| TCGA-COAD | 286 | 0.98 (0.97–0.98) |

Each value exceeds its cohort's random-gene-set baseline: the upper 97.5 % bound of 200 random sets of the same size was 0.05, −0.15, 0.26 and 0.76 respectively. The TCGA baseline (median 0.64) shows that RNA-seq gene-set scores share a strong sample-wide factor, so the TCGA correlation is less informative than its magnitude suggests. Within each of six CIT subtypes in GSE39582, ρ was 0.85–0.94.

**3.3 Genotype and the p53 and Wnt scores**

| Effect | GSE39582 | External |
|---|---|---|
| p53 score, TP53 mutant vs wild-type | −0.45 SD (−0.39 within pMMR) | GSE41258 −0.22 (−0.43 to −0.01); TCGA-COAD −0.42 (−0.54 to −0.29); GSE143985 −0.06 (n.s.) |
| p53 score, BRAF mutant | +0.89 SD unadjusted → +0.03 after MMR and composition adjustment | TCGA-COAD +0.23 raw, +0.12 (−0.27 to 0.49) excluding hypermutated tumours |
| Wnt score, BRAF mutant | −0.83 SD → −0.65 after adjustment | TCGA-COAD −0.42 raw, but −0.11 (−0.44 to 0.21) after excluding hypermutated tumours: confounded |

Lower Wnt expression was concentrated in canonical β-catenin targets (AXIN2, NKD1, TCF7). Gene modules within the p53 set did not meet the pre-registered stability criterion and are treated as hypotheses.

**3.4 A defective public dataset.** GSE103479 failed the quality gate: median control-pair r 0.07, against 0.78–0.94 elsewhere. Same-gene probes were uncorrelated, and adjacent probe identifiers carried near-identical values, consistent with values misaligned to probe IDs. Without the gate, this dataset would have suggested that the EMT–stroma relationship does not replicate.

## 4. Discussion

- **Prediction.** Adding four mechanistically motivated pathway scores to clinical variables did not measurably improve recurrence discrimination. Most of the apparent gain was optimism.
- **What the scores measure.** In bulk tissue, the scores encode composition and genotype-linked programmes:
  - the EMT signature behaves as a fibroblast/ECM index, consistent with Calon and Isella and now replicated across platforms;
  - p53 target expression tracks TP53 status;
  - some apparent genotype effects are confounded by mismatch repair and composition.
- **Implications.** A rising "EMT score" in bulk colorectal data should not be read as tumour-cell EMT. Associations between signatures and mutations should be adjusted for composition and MMR before interpretation.

**Limitations:**
- Composition was estimated computationally from the same arrays, not measured.
- Follow-up analyses were post hoc relative to the primary study.
- TCGA-COAD includes all stages.
- Pathway activity was not measured.
- No compatible external cohort allowed predictive validation.

## 5. Data and code availability

**Data.** All data are public:
- GEO: GSE39582, GSE41258, GSE143985, GSE103479;
- UCSC Xena: TCGA-COAD.

**Code.** Analysis code, registered plans and provenance records are in the project repository. `tools/pathway-audit/` provides the reusable tool.

## References

1. Marisa L, et al. *PLoS Med* 2013;10:e1001453.
2. Guinney J, et al. *Nat Med* 2015;21:1350–1356.
3. Calon A, et al. *Nat Genet* 2015;47:320–329.
4. Isella C, et al. *Nat Genet* 2015;47:312–319.
5. TCGA Network. *Nature* 2012;487:330–337.
6. Liberzon A, et al. *Cell Syst* 2015;1:417–425.
7. Yoshihara K, et al. *Nat Commun* 2013;4:2612.
8. Becht E, et al. *Genome Biol* 2016;17:218.
9. Hennig C. *Comput Stat Data Anal* 2007;52:258–271.
10. Riley RD, et al. *Stat Med* 2022;41:1280–1295.
