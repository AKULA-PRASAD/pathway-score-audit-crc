# Final research contribution

**Date:** 8 November 2026.
**Scope:** the original registered study, the post-hoc limitation investigation (X1–X9), the registered research extensions (E1–E3 plus the TCGA-COAD addendum), and the `pathway-audit` tool (`tools/pathway-audit/`).
**Supporting documents:**
- `docs/RESEARCH_QUESTION_REASSESSMENT.md`
- `docs/BIOLOGICAL_CONTRIBUTION_STRENGTHENING_PLAN.md`
- `docs/RESEARCH_EXTENSION_REGISTER.md`
- `docs/INDEPENDENT_SCIENTIFIC_CHALLENGE.md`

## 1. Original research question

> Do four predefined tumour pathway-expression scores (EMT, E2F targets, Wnt/β-catenin, p53) add reproducible recurrence-discrimination information beyond clinical variables in stage II–III colon cancer (GSE39582)?

## 2. Biological motivation

These are four axes of colorectal-cancer biology, one gene set each:
- Wnt/β-catenin initiation, with APC/CTNNB1 altered in about 80 % of tumours;
- E2F-driven proliferation;
- the p53 checkpoint, with TP53 mutated in about 60 % of non-hypermutated tumours;
- the mesenchymal programme linked to the poor-prognosis CMS4 subtype.

The panel was constrained and fixed before any outcome was examined, because 136 recurrences cannot support transcriptome-wide selection.

## 3. Computational approach

- Deposited RMA + ComBat arrays (585 arrays → 452 patients).
- Probe mapping: 54,675 → 43,163 → 21,815 genes.
- Fold-internal Hallmark scores.
- Unpenalised clinical Cox model (M1) vs M1 plus four ridge-penalised scores (M2).
- An optimism-corrected ΔC estimand with a paired nested bootstrap: 750 outer × 100 inner resamples, five Monte Carlo replicates, and a pre-specified precision rule.
- A hash-frozen pipeline, a verifier gate (60/60) and an independent recomputation (138 quantities).

Afterwards, two registered post-hoc programmes followed, each with independent review:
- X1–X9: limitation investigation;
- E1–E3: what the scores measure, with external replication.

## 4. What was newly investigated (E1–E3)

All outcome-free, registered before analysis (plan SHA-256 `c30dd5f4…`):
- **E1:** gene-level architecture of each score.
- **E2:** score relationships conditioned on estimated composition, genotype and subtype.
- **E3:** replication in three independent cohorts with genotype data. One cohort (GSE103479) was found to carry a probe-label defect and was excluded.

## 5. Strongest biological findings, with evidence class

1. **In bulk colorectal tumours, the Hallmark "EMT" score behaves as a stromal/ECM programme, not a tumour-cell EMT readout. This replicates in three independent cohorts across two platforms.**
   - Gene-disjoint correlation with fibroblast markers: ρ 0.943 (GSE39582), 0.950 (GSE41258), 0.939 (GSE143985), 0.979 (TCGA-COAD, RNA-seq).
   - Every cohort's EMT correlation exceeds its own random-gene-set baseline (200 random sets of the same size; upper 97.5 % bounds: GSE39582 0.05, GSE41258 −0.15, GSE143985 0.26, TCGA-COAD 0.76). In TCGA RNA-seq, random gene sets already correlate about 0.64 with the fibroblast estimate, so the TCGA ρ of 0.98 is less informative than its size suggests.
   - It holds within every CIT subtype (ρ 0.85–0.94).
   - *Class:* supported by independent external evidence (programme dominance); supported by computational estimates (fibroblast content).
2. **The p53 target-gene score tracks TP53 genotype.**
   - −0.45 SD in TP53-mutant tumours, −0.39 SD within pMMR (GSE39582).
   - **Replicated (pre-registered rule):** GSE41258 rank-biserial −0.22 (CI −0.43 to −0.01); TCGA-COAD −0.42 (−0.54 to −0.29), robust to composition adjustment and to hypermutation (−0.37 excluding hypermutated tumours); GSE143985 null (−0.06).
   - *Class:* directly demonstrated in GSE39582; supported by independent external evidence (two of three cohorts).
3. **Two apparent genotype effects separate once MMR status and estimated composition are accounted for.**
   - The p53-score elevation in BRAF-mutant and CIMP+ tumours disappears after adjustment for dMMR and estimated stromal/immune content (+0.89 → +0.03 SD).
   - The lower Wnt score in those tumours persists (−0.83 → −0.65 SD). It concentrates in canonical β-catenin targets (AXIN2, NKD1, TCF7, NOTCH1, PTCH1; rank-biserial −0.80).
   - In TCGA-COAD the raw BRAF effects point the same way (Wnt −0.42, p53 +0.23). However, 31 of 44 BRAF-mutant tumours are hypermutated, and after excluding hypermutated tumours or adjusting for mutation burden both effects vanish (Wnt −0.11, CI −0.44 to 0.21; p53 +0.12, CI −0.27 to 0.49). An independent re-check identified this, and the TCGA BRAF results are classified **confounded**, not replicated.
   - *Class:* exploratory association in GSE39582 only.

The p53 set appears to contain a stromal/immune-like gene group (CTSD, HMOX1, TGFB1; ρ 0.75 with stromal estimates) and a canonical-target group (DDB2, FDXR, MDM2, ZMAT3; strongest TP53 and MMR effects). The pre-registered stability criterion was **not met**, so this is a **hypothesis**.

## 6. Strongest methodological contributions

These are rigorous applications of established methods, not new algorithms:
- An optimism-corrected *difference* in C with paired nested resampling and an explicit Monte Carlo precision rule. The primary result stands on this estimand.
- **A reusable outcome-free audit of what a bulk pathway score measures:**
  - gene-level item–rest and leave-one-gene-out analysis;
  - fixed-rule modules with stability criteria;
  - composition-conditioned genotype effects;
  - random-gene-set nulls for conditioned correlations;
  - gene-disjoint marker tests;
  - pre-registered external replication with a fixed decision rule.
- **Data-quality detection:** positive-control co-expression QC identified an apparently misaligned public deposit (GSE103479). A naive replication would otherwise have reported spurious non-replication of the EMT–stroma relationship.
- Registration before analysis, independent re-implementation, and transparent recording of deviations throughout.

## 7. Predictive findings (unchanged)

- Corrected ΔC **+0.004316**; all five 95 % intervals include zero (about −0.009 to +0.05).
- Calibration slopes 0.833 / 0.862 indicate residual overfitting.
- No primary pathway hazard ratio excludes 1.
- Post-hoc: composition added no discrimination (corrected C 0.601 vs 0.604), and competing-risk estimates agree with cause-specific ones.
- **No clinical utility is claimed.**

## 8. Independent replication achieved

| Hypothesis (pre-registered) | Verdict (QC-passing cohorts) |
|---|---|
| R6 EMT score ~ fibroblast markers | **Replicated** (GSE41258, GSE143985, TCGA-COAD) |
| R1 p53 score lower in TP53-mutant tumours | **Replicated** (GSE41258, TCGA-COAD; GSE143985 null) |
| R3 Wnt lower in BRAF-mutant tumours | Not established: TCGA effect explained by hypermutation (registered rule alone would say partially replicated) |
| R4 p53 higher in BRAF-mutant tumours | Not established: TCGA effect explained by hypermutation |
| R2 E2F higher in TP53-mutant tumours | Not replicated (0 of 3) |
| R5 p53 higher in MSI-high tumours | Not replicated (GSE41258: +0.20, CI −0.04 to 0.43; TCGA MSI annotation too sparse) |

TCGA-COAD was added under a separately registered addendum (`docs/RESEARCH_EXTENSION_ADDENDUM_TCGA.md`). It covers all stages, whereas GSE39582 is stage II–III.

No external **predictive** validation was possible: no compatible cohort exists publicly (X8).

## 9. Remaining uncertainties

- Pathway *activity* was not measured.
- Tumour-cell-specific expression and tumour-cell EMT cannot be determined from bulk data.
- Composition is estimated computationally, not measured.
- The provenance of the ComBat `combat_mod` covariate is undocumented. ComBat protected it, which may inflate MMR-group differences.
- Baseline batch heterogeneity is real and of unknown source.
- The recurrence time unit is undocumented.
- All extension findings are post-hoc.

## 9b. A usable solution: the `pathway-audit` tool

The project's methods are packaged as an installable Python tool with a command-line interface (`tools/pathway-audit/`). Given any bulk expression matrix and gene sets, it reports:
- whether the data pass a positive-control quality gate;
- which genes drive each score;
- how much of each score is tissue composition;
- whether genotype effects survive adjustment;
- optionally, whether the scores add optimism-corrected discrimination to a clinical model.

**What has been checked:**
- 9 of 9 tests pass, including validation of the ESTIMATE re-implementation, in a fresh environment.
- Five-cohort case study (`reports/pathway_audit_case_study/`):
  - it flags EMT as composition-dominated in every usable cohort;
  - it fails GSE103479 at the quality gate;
  - its survival screen (+0.0048) is close to the frozen nested estimate (+0.0043).

**Not yet run:** the Dockerfile and the CI workflow, because Docker is not installed here and nothing has been pushed.

## 10. Contribution to computational biology

- A carefully established null result for adding pathway scores to clinical recurrence models, with an estimand designed for comparing models.
- Evidence, replicated in independent tumours, that a widely used "EMT" signature in bulk colorectal tissue is predominantly a stromal programme. This affects how such signatures should be interpreted in any bulk study.
- A demonstration that apparent genotype–pathway associations in bulk data can be explained by MMR status and tissue composition (p53 score), or survive them (Wnt targets).
- A reproducible, outcome-free audit framework for pathway scores, released as the `pathway-audit` tool, including data-quality controls that caught a defective public deposit.

## 11. Next logical study

- **Most informative:** spatial transcriptomics, or single-cell / dual-IHC (cytokeratin with vimentin/ZEB1; nuclear β-catenin; p21) on a subset of stage II–III tumours with matched bulk arrays. This would separate epithelial from stromal mesenchymal expression and test whether the Wnt-target and p53-target modules reflect tumour-cell activity.
- **Computational:** apply `pathway-audit` to further cohorts (CPTAC proteogenomic colon data would allow protein-level checks of the p53 and Wnt targets); obtain a corrected or raw-CEL version of GSE103479.
- **Predictive:** external validation of the frozen M1/M2 only in a cohort with about 340 or more recurrences and all five predictors.

## 12. Communication levels (evidence-grounded)

**Portfolio contribution (2–3 sentences).**
In 452 stage II–III colon cancers (GSE39582), I tested whether four predefined pathway-expression scores improve recurrence prediction beyond clinical variables. Using a nested-bootstrap estimate of the optimism-corrected gain, they did not (ΔC +0.004; all intervals include zero). Pre-registered follow-up work showed that the Hallmark EMT score in bulk tumours is dominated by a stromal/ECM programme, replicated in three independent cohorts (gene-disjoint ρ 0.94–0.98). I packaged the checks as `pathway-audit`, a tested command-line tool that flags composition-driven scores and broken datasets before they are misinterpreted.

**For a computational biology hiring manager.**
- **Design.** The primary study compared an unpenalised clinical Cox model with the same model plus four ridge-penalised Hallmark scores. The gene sets were fixed before any outcome was seen, because 136 events cannot support transcriptome-wide selection. Discrimination was estimated as an optimism-corrected ΔC with a paired nested bootstrap (750 × 100, five Monte Carlo replicates, pre-specified precision rule). The result was a small, non-reproducible gain (+0.0043; intervals −0.009 to +0.05), with calibration slopes below 1.
- **What the scores measure.** I then asked what the scores measure, in work registered before analysis and challenged by an independent re-implementation:
  - The EMT score correlates 0.94–0.98 with fibroblast markers even when no genes are shared, within every molecular subtype, and in three independent cohorts including TCGA RNA-seq. In bulk tissue it is a stromal-programme index.
  - The p53 target score is lower in TP53-mutant tumours (−0.45 SD), replicated in GSE41258 and TCGA-COAD.
  - Its apparent elevation in BRAF-mutant tumours disappears after adjusting for MMR status and estimated immune/stromal content. By contrast, lower expression of canonical β-catenin targets (AXIN2, NKD1) in those tumours persists.
- **Data quality.** A positive-control co-expression check showed that one public deposit has probe values misaligned to their labels. I excluded it rather than report spurious non-replication.
- **The tool.** I packaged all of this as `pathway-audit`, a tested Python command-line tool that produces this audit for any bulk dataset.
- **What is not claimed.** Pathway activity, tumour-cell-specific biology or clinical utility. Several module-level ideas failed their pre-set stability criteria and are reported as hypotheses.

**60–90 second interview version.**
"My project asked whether four classic colorectal-cancer pathway signatures improve recurrence prediction beyond stage and other clinical factors, in about 450 patients. I fixed the signatures before looking at outcomes and measured the gain with an optimism-corrected nested bootstrap, because apparent gains in small cohorts are mostly overfitting. The honest answer was no: the corrected gain was 0.004 in C, and the interval included zero.

Rather than stop there, I asked what these scores actually measure in bulk tumour tissue. The 'EMT' score turned out to be essentially a stromal or fibroblast signal. It correlates about 0.95 with fibroblast markers even when they share no genes, and that replicated in three independent cohorts, including TCGA. The p53 score tracks TP53 mutation as you'd expect. But its apparent link to BRAF mutation disappears once you account for mismatch-repair status and immune content, whereas lower Wnt-target expression in those tumours persists.

I registered every analysis in advance, had the key numbers independently re-implemented, and caught a public dataset whose probe values were misaligned. Then I turned the whole workflow into a command-line tool, `pathway-audit`, so anyone can run these checks on their own data. What I took from it is that a pathway score in bulk data is a readout of tissue composition and subtype as much as of the pathway, and you have to test that before interpreting it."

## Deviations from the registered extension plan

1. **E1 stability metric.** ARI proved degenerate for single-module sets (EMT median 0.0). A clusterwise Jaccard (Hennig 2007) was added post hoc as description. The registered ARI decision (H-E1 not evaluable) stands.
2. **E3 positive-control QC.** Added after GSE103479 gave implausible results. The cohort fails (median r 0.07). Verdicts are reported with and without it; the QC-passing verdicts are treated as primary. The defect was diagnosed independently as probe-label misalignment.
3. **E3 GSE41258 de-duplication.** 186 samples from 183 patients. Results were re-estimated with one tumour per patient; no verdict changed.
4. **Challenge tests** (`e_challenge.py`): partial-correlation null, gene-disjoint R6, Wnt-module dMMR adjustment. These were added after independent review; all are labelled.
5. **TCGA-COAD addendum.** Registered separately before analysis (hash in `reports/research_extensions/registration_tcga.json`). TCGA covers all stages.
6. **TCGA independent re-check (post-review).** Hypermutation sensitivity and cohort-specific random-gene-set baselines (`scripts/research_extensions/e3d_tcga_challenge.py` → `reports/research_extensions/tcga_challenge_tests.json`). R3/R4 reclassified as confounded; R1 robust; R6 holds against every cohort's baseline.
