# Speaker notes: GSE39582 colon-cancer recurrence study

Generated from the same slide specification as the presentation (`portfolio/src/make_deck.py`). Slides 1–19 are the main talk; A1 and A2 are appendix slides.

## Slide 1. Do tumour pathway-expression scores improve colon-cancer recurrence prediction?

*Section:* GSE39582 recurrence study  
*Source:* Illustration: original conceptual drawing, not patient data.

Introduce yourself and the question in one sentence: after surgery for stage II–III colon cancer, do four predefined tumour pathway-expression scores add recurrence-discrimination information beyond standard clinical variables? Say up front that the answer is a carefully established null: the pathway scores did not demonstrate a reproducible improvement, and this does not establish that their true effect is zero. The illustration on the right is conceptual: a normal crypt beside tumour glands in stroma. It previews the main biological caveat of the talk, that a bulk tumour sample mixes many cell types.

## Slide 2. The clinical problem and the research question

*Section:* Clinical problem  
*Source:* Cohort counts: reports/figures/final/figure_values.json. Question: docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md.

Stage is the main clinical risk sorter after surgery, but it is imperfect, which motivates molecular markers. The trap is optimism: any model evaluated on the patients it was fitted to looks better than it will on new patients, and gene-expression signatures are especially exposed. So the design fixes every feature before outcomes are examined, compares the clinical model with the clinical-plus-pathway model, and reports the gain after optimism correction with its uncertainty. The endpoint is the cause-specific hazard of recurrence; deaths without recurrence are censored, so no absolute risk is reported. The CIT partition is internal, never external validation.

## Slide 3. Four biological axes behind the pathway scores

*Section:* Biological background  
*Source:* Established mechanisms (textbook background); gene sets: MSigDB Hallmark (Liberzon 2015).

This map shows how the four axes connect, using standard arrow conventions: pointed arrows mean activation, flat bars mean inhibition. Wnt signalling inhibits the destruction complex; the complex inhibits beta-catenin; beta-catenin with TCF/LEF switches on MYC and cyclin D. Cyclin D–CDK4/6 phosphorylates and inactivates RB; RB normally represses E2F, so RB inactivation releases E2F and its G1/S genes. p53 responds to DNA damage, induces p21, which blocks CDK activity, and induces MDM2, which degrades p53 in a negative feedback loop. On the right, TGF-beta drives SNAIL, ZEB and TWIST, which repress E-cadherin and induce mesenchymal genes; fibroblasts express many of the same genes. All of this is background. The study measured only one averaged mRNA score per axis.

## Slide 4. Wnt/β-catenin: β-catenin stability controls TCF target genes

*Section:* Pathway mechanism 1 of 4  
*Source:* Literature: TCGA 2012 (R6); Guinney 2015 (R1). Measured: phase9 execution and S-7 results; PHASE8C §3.1. Post-hoc: reports/exploratory/x1.

With APC intact, the destruction complex phosphorylates beta-catenin and it is degraded, so target genes stay off. In most colorectal cancers APC is lost or CTNNB1 is mutated; beta-catenin escapes degradation, enters the nucleus with TCF/LEF and drives targets such as MYC and AXIN2. The Wnt score is the mean of 41 target genes. Its primary hazard ratio is 1.016 (0.906–1.138), essentially 1, and its genes are loosely co-expressed (coherence 0.041). A near-universal genomic alteration does not mean the mRNA score varies in a prognostically useful way between patients, but that explanation was not tested. S-7 re-summarises the same arrays with representative probes, so it is not replication. A post-hoc analysis on the same data found the Wnt score lower in BRAF-mutant (-0.74 SD) and CIMP-positive tumours, matching the known serrated route; that is one correlated axis and an association, not a measure of β-catenin activity.

## Slide 5. E2F: RB phosphorylation releases the G1/S programme

*Section:* Pathway mechanism 2 of 4  
*Source:* Literature: Anjomshoaa 2009 (R4); Luo 2019 (R5). Measured: phase9_execution_results.json; phase9_s7_results.json; PHASE8C §3.1.

Growth signals induce cyclin D–CDK4/6, which phosphorylates RB. Phosphorylated RB releases E2F, which transcribes cyclin E, the MCM helicase, PCNA and other replication genes. p21 from the p53 pathway restrains this step. The E2F score is the most coherent of the four (coherence 0.304). Its primary hazard ratio, 0.906 (0.807–1.016), sits below 1 but the interval includes 1. The five-set model gives 0.862 (0.746–0.996), which excludes 1 only after the collinear PI3K score was added; that model was fixed in advance as a sensitivity analysis only, so it is not a finding. The literature on proliferation and outcome is mixed, so no direction is claimed.

## Slide 6. p53: a stress checkpoint held in check by MDM2

*Section:* Pathway mechanism 3 of 4  
*Source:* Literature: TCGA 2012 (R6). Measured: phase9 execution and S-7 results; phase8c/score_vs_annotation.csv. Post-hoc: reports/exploratory/x1, x4.

DNA damage activates ATM and ATR, which stabilise p53. p53 induces CDKN1A (p21) for cell-cycle arrest, BAX and PUMA for apoptosis, and GADD45 for repair. It also induces MDM2, which degrades p53: a negative feedback loop. The score is the mean of 194 target genes, so it is not a mutation or protein-function assay. Primary hazard ratio 0.995 (0.882–1.123). It has the strongest association with mismatch-repair status of the four scores (epsilon squared 0.1716, exploratory, n = 410), and it changes sign when some single batches are omitted. The direction and mechanism of the MMR association were not tested. Post-hoc, on the same data: the p53 score is lower in TP53-mutant tumours (-0.45 SD), and the difference holds within mismatch-repair-proficient tumours (-0.39 SD, batch-stratified permutation p 0.0003). This replicated in independent cohorts: GSE41258 -0.22 and TCGA-COAD -0.42 (rank-biserial). That links the target-gene score to genotype as expected; it is not a measure of p53 protein function.

## Slide 7. EMT: tumour-cell programme or stromal signal?

*Section:* Pathway mechanism 4 of 4  
*Source:* Literature: Guinney 2015 (R1); Calon 2015 (R2); Isella 2015 (R3). Measured: phase9 results. Post-hoc: reports/exploratory/x5, x6.

Two different biological sources produce the same mesenchymal transcripts. Tumour cells undergoing EMT lose E-cadherin and gain vimentin and fibronectin under SNAIL, ZEB and TWIST. Cancer-associated fibroblasts express those genes without any tumour-cell change. In bulk colorectal tumours, the mesenchymal-subtype signal is largely stromal according to Calon and Isella. The EMT hazard ratio, 1.070 (0.951–1.202), stays above 1 in every fit but is never distinguishable from 1. Cell composition was not measured. Post-hoc, it was estimated computationally from the same arrays: the EMT score correlates 0.95 with estimated fibroblast content, and 0.93 even after removing the 29 genes it shares with the stromal marker sets. That is consistent with a predominantly stromal source, as Calon and Isella reported, but it does not establish the source in either direction. In registered follow-up, the same gene-disjoint relationship replicated in GSE41258, GSE143985 and TCGA-COAD (ρ 0.95, 0.94, 0.98). An elevated bulk EMT-associated score does not establish tumour-cell EMT.

## Slide 8. What a bulk tumour sample measures

*Section:* What was measured  
*Source:* Conceptual illustration. Estimates: reports/exploratory/x5 (ESTIMATE v1.0.13; MCP-counter v1.2.0), post-hoc.

Every array is a piece of bulk tumour: tumour epithelium, fibroblasts, immune cells and vessels, all homogenised into one mRNA pool. A pathway score is the average of its member genes' standardised values in that pool, so it cannot tell which cells expressed the genes. Cell composition was not measured. In post-hoc analyses it was estimated with two published reference methods, ESTIMATE and MCP-counter, which agree with each other (stromal versus fibroblast estimate ρ 0.93). These are relative scores computed from the same arrays, not measured cell fractions, and they carry measurement error. All four pathway scores also associate with MMR status in outcome-free tests in 410 patients.

## Slide 9. GSE39582: from 585 arrays to 452 analysed patients

*Section:* Dataset and cohort  
*Source:* Figure 1 (reports/figures/final/Fig1_study_design_cohort.png). Censoring profile: appendix A2.

Walk the funnel: 585 arrays, the stage II–III primary cohort of 453, and the complete-case analysis population of 452. The model estimates the cause-specific hazard of recurrence. Deaths without recurrence are censored; appendix A2 shows those patients were older, less often stage III and less often treated, so censoring is informative. That is why no absolute risk or cumulative incidence is reported. The CIT partition shares processing history with the discovery set and is confounded with batch, so it is internal.

## Slide 10. From 54,675 probe sets to four pathway scores

*Section:* Preprocessing and features  
*Source:* figure_values.json (fig2.probes); config/pathways/m2_panel_4set.json and m2_panel_5set.json; PHASE6C §5.

The data were used as deposited, without re-normalisation. 8,893 probe sets were unmapped and 2,619 ambiguous; the remaining 43,163 were collapsed to 21,815 genes by the median. The four Hallmark sets contain 607 genes. Each score is the mean of member-gene z-scores, with standardisation re-estimated inside each training fold so held-out patients never shape the features. On the right: the ten candidate sets were highly collinear (condition number 77). One representative per correlation cluster was chosen without outcome data, which brings the four-set panel down to 2.95. The five-set panel, adding PI3K/AKT/mTOR, was fixed in advance as a sensitivity analysis only.

## Slide 11. Model architecture: clinical benchmark and molecular extension

*Section:* Model design  
*Source:* docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md §3; λ from phase9_execution_results.json.

M1 is the clinical benchmark and was never re-tuned. M2 adds the four scores with a ridge penalty on the molecular block only. Leaving the clinical block unpenalised keeps the comparison scale-equivariant and makes M1 exactly nested in M2, so any gain is attributable to the pathway block. Shrinking only that block makes the gain in C conservative. In M1, stage III versus II has hazard ratio 1.97 (1.30–2.98); chemotherapy is confounded by indication and never read causally. The five-set model is shown separately because it is sensitivity only and was never promoted.

## Slide 12. Validation: optimism-corrected nested bootstrap

*Section:* Validation design  
*Source:* docs/PHASE9_DELTAC_INTERVAL_PLAN.md §3; figure_values.json (fig8, figS1); verifier output.

Each of the 750 outer resamples refits both models and computes the apparent gain in C on that resample. Inside it, 100 inner resamples estimate optimism: each inner refit is scored on its own resample and on the outer resample, with M1 and M2 paired on the same resamples. Corrected gain equals apparent gain minus the difference in optimism, and the 95 % interval is the percentile interval of the corrected gains. The five runs differ only in random seed: they are Monte Carlo replicates of one analysis, not independent experiments. The run size met a pre-specified precision rule (largest endpoint Monte Carlo error 4.2 % of interval width, limit 5 %).

## Slide 13. Central result: no demonstrated gain in discrimination

*Section:* Central result  
*Source:* figure_values.json (fig7, fig9); reports/phase9/corrected_rerun_report/corrected_interval_results.json.

The pathway model's apparent lead was +0.0140 in C. Most of it was optimism: the corrected gain is +0.004316. The five run intervals are run 0: [-0.0086, +0.0519]; run 1: [-0.0090, +0.0469]; run 2: [-0.0093, +0.0484]; run 3: [-0.0091, +0.0518]; run 4: [-0.0091, +0.0494]. All include zero. A gain of about 0.004 in C is very small in practical terms, and the interval cannot exclude a small loss or a modest gain. This does not establish that the true effect is zero. Corrected calibration slopes of 0.833 and 0.862 show both models overfit.

## Slide 14. Pathway hazard ratios: every primary interval includes 1

*Section:* Pathway estimates  
*Source:* phase9_execution_results.json (M2 primary, M2 sensitivity); phase9_s7_results.json.

Hazard ratios are per training-fold standard deviation of each score, on a log axis with the no-effect line at 1. No primary-model interval excludes 1. EMT stays above 1 and E2F below 1 in every fit, which is compatible with parts of the literature but is not a demonstrated association. S-7 re-summarises the same arrays, so it is a robustness check on probe summarisation, not replication. The five-set E2F estimate 0.862 (0.746–0.996) is shown separately because the model is sensitivity only; its interval excluding 1 coincided with collinearity.

## Slide 15. Robustness: batch sensitivity, subtype association, calibration

*Section:* Robustness  
*Source:* phase9_execution_results.json (batch_sensitivity); score_vs_annotation.csv; figure_values.json (fig9).

Left: each point is one leave-one-batch-out refit, plotted against how many patients were omitted. Wnt changes sign once and p53 twice; EMT and E2F keep their sign. Because batches range from 14 to 103 patients, batch identity and sample size are mixed. Batch location effects did not survive stratification by the ComBat covariate combat_mod, which matches recorded MMR in all 410 samples but has unresolved provenance, so it is not described as an MMR variable. Batch effects were characterised, not resolved. Top right: all four scores associate with MMR status (exploratory); these associations cannot be compared with hazard ratios. Bottom right: corrected calibration slopes below 1 indicate residual overfitting.

## Slide 16. What this project established

*Section:* Contributions  
*Source:* tools/pathway-audit/; docs/RESEARCH_EXTENSION_REGISTER.md + TCGA addendum. All limitations: docs/SCIENTIFIC_LIMITATION_RESOLUTION_REGISTER.md.

Four completed contributions, strongest first. The follow-up analyses behind blocks 2 and 3 are post-hoc relative to the primary study: they were registered before they ran and used no outcome data, but none of this is confirmatory. Block 1, the tool: pathway-audit is a Python package and command-line tool that checks data quality, gene drivers and coherence, tissue composition, molecular associations with covariate adjustment, a random-gene-set null and an optional survival screen. Its 11 package tests pass, the built-in demo recovers all four planted truths, and the case study covers five colorectal cohorts, one of which, GSE103479, failed the quality gate. A Dockerfile and a GitHub Actions workflow are written but have not been run, so I do not claim they pass. Block 2: with the 6 shared genes removed, the bulk EMT score correlates 0.94, 0.95, 0.94 and 0.98 with fibroblast markers in GSE39582, GSE41258, GSE143985 and TCGA-COAD. Each is compared with its own random-gene-set baseline, whose upper 97.5 % bounds are 0.05, -0.15, 0.26 and 0.76; in TCGA random sets already reach a median of about 0.64, so that cohort is less informative than its size suggests. This is a bulk-tissue association with estimated stromal content. Fibroblast content is a computational estimate, not a measurement, and tumour-cell EMT is not established. Block 3: the p53 target-gene score is lower in TP53-mutant tumours. Discovery effect in GSE39582 -0.45 SD; replication as rank-biserial effects: GSE41258 -0.22, TCGA-COAD -0.42, and -0.37 after excluding hypermutated tumours; GSE143985 was not significant, so it replicated in two of three external cohorts. This is an expression association; it does not measure p53 protein activity or pathway function. E2F–TP53 and p53–MSI did not replicate, and the TCGA BRAF links were explained by hypermutation; these are not findings. Block 4: the pre-specified nested bootstrap gives a corrected gain in C of +0.004316, and all five intervals include zero, so no reproducible improvement in recurrence discrimination was demonstrated (reports/phase9/corrected_rerun_report/corrected_interval_results.json). That does not mean the pathways have no biological role. Still unresolved: tumour-cell EMT, direct pathway activation, external predictive validation (14 public cohorts assessed, none eligible) and the recurrence-time unit, provisionally months. Batch heterogeneity (frailty variance 0.13), competing risks, MMR association and computational composition are partially addressed or characterised; the full status of every limitation is in the scientific limitation resolution register.

## Slide 17. Decisions and audit trail

*Section:* Decisions and audit trail  
*Source:* docs/OWNER_DECISION_REGISTER.md; docs/PHASE9_CORRECTED_FINAL_REPORT.md §13–§15.

Each decision is recorded in the project's decision records, DC1 to DC7; describe the reasoning for each rather than claiming to have discovered every problem. DC1 fixed every analytic choice before outcomes; DC2 kept clinical terms unpenalised; DC3 kept the CIT partition internal; DC4 adjudicated the random-seed deviation and authorised a full re-run; DC5 completed the pre-specified 500 to 750 top-up; DC6 kept misleading results out of the conclusions; DC7 commissioned and reviewed the audits. The separate audit recomputed reported quantities from the recorded checkpoints; it did not refit models and is not an external replication.

## Slide 18. From finding to solution: replication and the pathway-audit tool

*Section:* From finding to solution  
*Source:* docs/RESEARCH_EXTENSION_REGISTER.md and TCGA addendum (registered before analysis); tools/pathway-audit/; reports/pathway_audit_case_study/.

Because the prediction answer was negative, the follow-up asked what these scores measure in bulk tissue, in analyses registered before they ran and using no outcome data. The clearest result: the EMT score tracks fibroblast markers even with the 6 shared genes removed, in GSE39582 and in three independent cohorts, including TCGA RNA-seq. One public cohort, GSE103479, failed a positive-control check (control gene pairs correlated at only 0.07); its values appear misaligned to probe labels, so it was excluded rather than reported as non-replication. I then packaged the checks into pathway-audit, a tested Python command-line tool: a data-quality gate, gene drivers, composition, genotype effects adjusted for composition and MMR, and an optional survival screen. On GSE39582 its survival screen gives a corrected gain of +0.0048, close to the frozen nested estimate. An independent re-check of TCGA showed that the apparent BRAF effects there are explained by hypermutation; with mutation burden as a covariate, the tool reports exactly that, while the TP53 link to the p53 score stays robust. Random gene sets also correlate with the fibroblast estimate in TCGA (about 0.64), so each cohort is compared with its own random baseline. Composition values are computational estimates, not measured cell fractions.

## Slide 19. Conclusion and the next scientific questions

*Section:* Conclusions  
*Source:* docs/SCIENTIFIC_INTERPRETATION_AND_CONCLUSION.md; docs/FINAL_RESEARCH_CONTRIBUTION.md; docs/EXTERNAL_COHORT_ELIGIBILITY_AUDIT.md.

This slide pulls the project together rather than repeating the evidence. The project asked whether four pathway scores improve recurrence prediction, and found no reproducible gain: corrected gain in C +0.004316, with every interval including zero. That does not establish that the pathways have no biological role. It then asked what the scores measure in bulk tumour tissue. In registered, outcome-free follow-up, the EMT score tracked estimated stromal content and the p53 score tracked TP53 status; other genotype links did not hold. It delivered pathway-audit so these checks can be run before a score is interpreted. The take-home is practical: audit bulk pathway scores for composition, genotype confounding and data quality before reading them as biology or adding them to prognostic models. Three questions follow, none of them performed here. First, whether EMT-associated expression sits in tumour cells or stroma, which needs spatial, single-cell or cell-resolved protein data. Second, whether the pathways are actually active, which needs functional or protein-level readouts. Third, whether any score adds recurrence information in a new cohort: that needs an independent stage II–III cohort with about 337 recurrences; the largest of the 14 public cohorts assessed had 124.

## Slide 20. References and figure sources

*Section:* Appendix A1  
*Source:* References verified in docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md.

Reference slide for questions. All eight references were checked against the published records. Every visual in the deck is an original drawing generated from code; numbers are read from the frozen outputs.

## Slide 21. Quality control and informative censoring

*Section:* Appendix A2  
*Source:* Figures 2 and S2 (reports/figures/final/).

Backup for questions on data quality. Quality control used soft flags only: 29 arrays were flagged and none removed, because there was no independent evidence of array failure. Processing batch is confounded with the discovery/CIT partition. Figure S2 profiles the deaths without recurrence that are censored under the cause-specific endpoint: they were older, less often stage III and less often treated, so censoring is informative and no absolute risk is reported.
