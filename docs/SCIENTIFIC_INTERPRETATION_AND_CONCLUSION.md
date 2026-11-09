# Scientific interpretation and conclusion

**Date:** 8 November 2026.

**Evidence base:**
- the frozen primary analysis (`docs/PHASE9_CORRECTED_FINAL_REPORT.md`);
- the registered post-hoc exploratory analyses X1–X9 (`docs/ADDITIONAL_ANALYSES_METHODS_AND_RESULTS.md`);
- the limitation register (`docs/SCIENTIFIC_LIMITATION_RESOLUTION_REGISTER.md`).

Throughout, **biological plausibility** (textbook and literature) is kept separate from **what this cohort showed**.

## 1. Why EMT, E2F, Wnt/β-catenin and p53?

- With 136 recurrences, a transcriptome-wide signature would overfit. The project therefore fixed, before seeing outcomes, one representative Hallmark gene set per axis of colorectal-cancer biology:
  - **Wnt/β-catenin:** the initiating pathway in most colorectal cancers. APC or CTNNB1 is altered in about 80 %, and Wnt-pathway genes in 93 % (TCGA 2012).
  - **E2F targets:** the G1/S proliferation programme released when cyclin D–CDK4/6 phosphorylates RB.
  - **p53 pathway:** the DNA-damage checkpoint. TP53 is mutated in about 60 % of non-hypermutated tumours (TCGA 2012).
  - **EMT:** the mesenchymal programme linked to the poor-prognosis CMS4 subtype (Guinney 2015).
- Ten candidate sets were reduced to these four by correlation clustering, without outcome data.

## 2. What the scores represent in this cohort

Each score is the mean standardised mRNA of a gene set measured in bulk tumour tissue. The new analyses show what each one tracks.

| Score | Strongest observed correlates (exploratory) | Reading |
|---|---|---|
| **p53** | Lower in TP53-mutant tumours (−0.45 SD; holds within pMMR, −0.39 SD). Higher in the BRAF/CIMP+/CIN− cluster. Modest correlation with immune estimates (ρ 0.32–0.36) | Partly anchored to TP53 status, as expected if mutant p53 transactivates fewer targets. Part of its variation may come from immune/interferon-responsive genes in the set |
| **Wnt** | Lower in BRAF-mutant and CIMP+ tumours (−0.74 / −0.63 SD), higher in CIN+. Essentially unrelated to composition (\|ρ\| ≤ 0.18). Gene coherence very low (0.041) | Follows the known serrated vs conventional split. Weakly coordinated; it does not measure β-catenin activity |
| **E2F** | Higher in TP53-mutant and CIMP+ tumours. Lower with stromal content (ρ −0.32); higher with purity (ρ +0.30). The most coherent set (0.304) | A coherent proliferation programme, more evident in tumour-rich samples |
| **EMT** | ρ 0.95 with estimated fibroblast abundance and 0.91 with ESTIMATE stromal score; −0.80 with purity. Still 0.93 / 0.89 after removing all 29 marker genes it shares with stromal sets. Highest in the CIT C4 subtype | In this cohort the bulk EMT score behaves as a **stromal/fibroblast abundance index**. This agrees with reports that mesenchymal transcripts in bulk colorectal tumours are largely stromal (Calon 2015; Isella 2015) |

**Caveats on the table:**
- Every correlate comes from the same arrays.
- CIT subtypes were derived from these data, so the subtype associations are partly circular.
- Associations with mutations are links between expression and genotype, not measures of pathway activity.

## 3. Which associations are reproducible, and which are sensitive?

- **Robust within this cohort:**
  - p53 score with TP53 status (survives MMR restriction and batch-stratified permutation);
  - p53 score with MMR (survives batch stratification and holds within CIT C2);
  - EMT score with estimated stromal content (two independent estimators; survives removal of shared genes).
- **Sensitive to subtype composition:** the EMT, E2F and Wnt associations with MMR largely disappear within CIT C2, where 50 of 59 dMMR tumours lie.
- **Sensitive to batch:**
  - Batches differ in baseline recurrence hazard (frailty θ̂ 0.13, profile 95 % ≈ 0.03–0.50).
  - Pathway HRs are unchanged when this is modelled.
  - Leave-one-batch-out sign changes occur only for coefficients within 0.3 SE of zero.
- **Sensitive to composition and collinearity:** the EMT outcome association cannot be separated from estimated stromal content (VIF 12.8 / 19.0 with both in the model).
- **Sensitive to penalisation and MMR:**
  - E2F has HR 0.906 (0.807–1.016) in the frozen ridge model and 0.82 (0.69–0.98) unpenalised.
  - With dMMR adjustment it is 0.87 (0.72–1.05).
  - No pathway HR is a demonstrated prognostic association.
- **Five-set E2F, 0.862 (0.746–0.996):** sensitivity analysis only. It coincided with added collinearity and is not promoted.

## 4. Did the biology improve recurrence prediction beyond clinical variables?

**No reproducible improvement was demonstrated.**
- Apparent C rose from 0.6242 to 0.6382 (+0.0140).
- After optimism correction the gain was **+0.004316**.
- All five Monte Carlo replicate 95 % intervals include zero, spanning about −0.009 to +0.05.
- About 70 % of the apparent gain was optimism.

The exploratory analyses do not change this:
- Composition estimates add no discrimination beyond clinical variables (corrected C 0.601 vs 0.604).
- Composition-adjusted models differ by ≤ 0.013 in single-level corrected C, with no interval.
- The four-score block's unpenalised likelihood-ratio test is p = 0.056.
- Competing-risk modelling leaves every estimate essentially unchanged.

## 5. What the calibration evidence shows

- Corrected calibration slopes are 0.833 (M1) and 0.862 (M2), and 0.836 for the five-set model.
- Predicted risk differences are therefore too extreme in both models. Any clinical use would need shrinkage and recalibration.
- Calibration does not distinguish the pathway model from the clinical model.

## 6. Absolute risk and competing events

- Death without recurrence (46 patients; older, less often stage III, less often treated) is a genuine competing event.
- Recurrence cumulative incidence is 30.7 % at 60 time units overall, 22.9 % for stage II and 40.5 % for stage III (Aalen–Johansen).
- Naive Kaplan–Meier overstates it only slightly (31.3 %).
- Fine–Gray and cause-specific estimates agree.
- These absolute figures are **conditional on the time unit being months, which remains undocumented**.

## 7. What can be concluded about biological relevance

- The four scores capture **real, interpretable biological variation** in these tumours:
  - a TP53-linked target-gene programme;
  - a serrated-vs-conventional Wnt/p53 axis;
  - a coherent proliferation programme;
  - stromal content.
- They relate to molecular subtype and to MMR status.
- That variation **did not translate into demonstrable added recurrence discrimination** once clinical variables, including stage, were in the model and optimism was removed.
- This is consistent with prognostic information in bulk expression being concentrated in a few correlated axes, principally stromal content and the MSI/serrated axis. Those axes are only modestly prognostic within stage II–III after clinical adjustment, and are hard to estimate precisely with 136 events.

## 8. What cannot be concluded

- **That the pathways have no prognostic role.** The interval includes a modest gain (up to about +0.05).
- **That any pathway is active or inactive in tumour cells.** No protein, phosphorylation or functional data exist.
- **That tumour cells undergo EMT.** The bulk EMT score tracks stromal content; its cellular source cannot be determined from bulk data.
- **That cell composition was measured.** ESTIMATE and MCP-counter are computational estimates from the same arrays.
- **That the MMR associations are free of processing effects.** ComBat protected a covariate aligned with MMR, so group differences may be inflated.
- **That batch effects are resolved.** Baseline heterogeneity between batches is real and its source is unknown.
- **That the models transport to other populations.** No compatible independent cohort exists publicly (14 cohorts assessed; none eligible; about 340 events would be needed).
- **Absolute risk in calendar time.** The time unit is undocumented.

## 9. One-paragraph conclusion

In 452 patients with stage II–III colon cancer, four predefined Hallmark pathway-expression scores captured biologically coherent variation:
- a p53 target programme lower in TP53-mutant tumours;
- a Wnt/p53 axis separating serrated from conventional tumours;
- a proliferation programme;
- an EMT score that, in bulk tissue, behaves as an index of stromal content.

Adding the four scores to a clinical Cox model did not demonstrate a reproducible improvement in recurrence discrimination (optimism-corrected ΔC +0.004316; all intervals include zero). No pathway hazard ratio was a demonstrated association. Post-hoc analyses of competing risks, batch frailty, estimated cell composition and molecular anchoring sharpened the biological interpretation, but none changed this conclusion. External validation was not feasible with any public cohort identified.
