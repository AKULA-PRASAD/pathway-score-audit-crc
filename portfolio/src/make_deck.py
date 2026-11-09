"""Build the presentation (editable PPTX + PDF) and the speaker-notes file: 19 main slides + 2 appendix slides.

Usage: python portfolio/src/make_deck.py   (requires python-pptx; see portfolio/requirements-portfolio.txt)
Run portfolio/src/make_slide_visuals.py first; it draws the slide visuals from frozen values.
One slide specification drives both renderers, so the PPTX and PDF carry identical text.
Numbers come from values.py (frozen outputs, read-only). Fits no model.
"""
import os
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from pptx import Presentation  # noqa: E402
from pptx.dml.color import RGBColor  # noqa: E402
from pptx.util import Inches, Pt  # noqa: E402

import values  # noqa: E402

R = values.ROOT
FIG = R / "reports" / "figures" / "final"
VIS = R / "portfolio" / "diagrams" / "slides"
OUT = R / "portfolio" / "presentation"
W, H = 13.333, 7.5
INK, MUTED, PRIMARY, PANEL, GRID = "#1f2933", "#52606d", "#2a78d6", "#f5f7fa", "#e3e8ee"
PREVIEW = os.environ.get("DECK_PREVIEW_DIR")
AUTHOR, AFFIL, DATE = "Prasad Akula", "MS in Bioinformatics, Northeastern University", "November 2026"
plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "pdf.fonttype": 42, "text.color": INK})


def hr(t):
    return f"{t[0]:.3f} ({t[1]:.3f}–{t[2]:.3f})"


def build_slides(V):
    c, p = V["cohort"], V["probes"]
    m = {x["label"]: x for x in V["pathways"]}
    runs = "; ".join(f"run {i}: [{l:+.4f}, {u:+.4f}]" for i, (l, u) in enumerate(V["runs"]))
    d_app = V["C"]["M2_primary"]["apparent"] - V["C"]["M1"]["apparent"]
    st = V["stage_hr"]
    x = values.load_exploratory()
    E = values.load_extensions()
    SC = values.load_single_cell()
    S = []
    S.append(dict(kind="cover", section="GSE39582 recurrence study", image=VIS / "V01_cover_illustration.png",
        title="Do tumour pathway-expression scores improve colon-cancer recurrence prediction?",
        subtitle=f"An optimism-corrected test of four Hallmark pathway scores in {c['n']} stage II–III colon cancers (GSE39582), a four-cohort replication, and the pathway-audit tool",
        source="Illustration: original conceptual drawing, not patient data.",
        notes=("Introduce yourself and the question in one sentence: after surgery for stage II–III colon cancer, do four predefined tumour pathway-expression scores add recurrence-discrimination information beyond standard clinical variables? "
               "Say up front that the answer is a carefully established null: the pathway scores did not demonstrate a reproducible improvement, and this does not establish that their true effect is zero. "
               "The illustration on the right is conceptual: a normal crypt beside tumour glands in stroma. It previews the main biological caveat of the talk, that a bulk tumour sample mixes many cell types.")))
    S.append(dict(kind="visual", section="Clinical problem", image=VIS / "V02_problem_question.png",
        title="The clinical problem and the research question",
        source="Cohort counts: reports/figures/final/figure_values.json. Question: docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md.",
        notes=("Stage is the main clinical risk sorter after surgery, but it is imperfect, which motivates molecular markers. "
               "The trap is optimism: any model evaluated on the patients it was fitted to looks better than it will on new patients, and gene-expression signatures are especially exposed. "
               "So the design fixes every feature before outcomes are examined, compares the clinical model with the clinical-plus-pathway model, and reports the gain after optimism correction with its uncertainty. "
               "The endpoint is the cause-specific hazard of recurrence; deaths without recurrence are censored, so no absolute risk is reported. The CIT partition is internal, never external validation.")))
    S.append(dict(kind="visual", section="Biological background", image=VIS / "V03_biology_map.png",
        title="Four biological axes behind the pathway scores",
        subtitle=f"Measured here: one bulk mRNA score per axis (EMT {m['EMT']['genes']}, E2F {m['E2F targets']['genes']}, Wnt {m['Wnt/β-catenin']['genes']}, p53 {m['p53 pathway']['genes']} genes). Pathway activity itself was not measured.",
        source="Established mechanisms (textbook background); gene sets: MSigDB Hallmark (Liberzon 2015).",
        notes=("This map shows how the four axes connect, using standard arrow conventions: pointed arrows mean activation, flat bars mean inhibition. "
               "Wnt signalling inhibits the destruction complex; the complex inhibits beta-catenin; beta-catenin with TCF/LEF switches on MYC and cyclin D. "
               "Cyclin D–CDK4/6 phosphorylates and inactivates RB; RB normally represses E2F, so RB inactivation releases E2F and its G1/S genes. "
               "p53 responds to DNA damage, induces p21, which blocks CDK activity, and induces MDM2, which degrades p53 in a negative feedback loop. "
               "On the right, TGF-beta drives SNAIL, ZEB and TWIST, which repress E-cadherin and induce mesenchymal genes; fibroblasts express many of the same genes. "
               "All of this is background. The study measured only one averaged mRNA score per axis.")))
    S.append(dict(kind="visual", section="Pathway mechanism 1 of 4", image=VIS / "V04_wnt.png",
        title="Wnt/β-catenin: β-catenin stability controls TCF target genes",
        source="Literature: TCGA 2012 (R6); Guinney 2015 (R1). Measured: phase9 execution and S-7 results; PHASE8C §3.1. Post-hoc: reports/exploratory/x1.",
        notes=("With APC intact, the destruction complex phosphorylates beta-catenin and it is degraded, so target genes stay off. "
               "In most colorectal cancers APC is lost or CTNNB1 is mutated; beta-catenin escapes degradation, enters the nucleus with TCF/LEF and drives targets such as MYC and AXIN2. "
               f"The Wnt score is the mean of {m['Wnt/β-catenin']['genes']} target genes. Its primary hazard ratio is {hr(m['Wnt/β-catenin']['hr'])}, essentially 1, and its genes are loosely co-expressed (coherence {m['Wnt/β-catenin']['coherence']:.3f}). "
               "A near-universal genomic alteration does not mean the mRNA score varies in a prognostically useful way between patients, but that explanation was not tested. "
               "S-7 re-summarises the same arrays with representative probes, so it is not replication. "
               f"A post-hoc analysis on the same data found the Wnt score lower in BRAF-mutant ({x['braf_wnt_diff']:+.2f} SD) and CIMP-positive tumours, matching the known serrated route; that is one correlated axis and an association, not a measure of β-catenin activity.")))
    S.append(dict(kind="visual", section="Pathway mechanism 2 of 4", image=VIS / "V05_e2f.png",
        title="E2F: RB phosphorylation releases the G1/S programme",
        source="Literature: Anjomshoaa 2009 (R4); Luo 2019 (R5). Measured: phase9_execution_results.json; phase9_s7_results.json; PHASE8C §3.1.",
        notes=("Growth signals induce cyclin D–CDK4/6, which phosphorylates RB. Phosphorylated RB releases E2F, which transcribes cyclin E, the MCM helicase, PCNA and other replication genes. p21 from the p53 pathway restrains this step. "
               f"The E2F score is the most coherent of the four (coherence {m['E2F targets']['coherence']:.3f}). Its primary hazard ratio, {hr(m['E2F targets']['hr'])}, sits below 1 but the interval includes 1. "
               f"The five-set model gives {hr(V['five_set_e2f'])}, which excludes 1 only after the collinear PI3K score was added; that model was fixed in advance as a sensitivity analysis only, so it is not a finding. "
               "The literature on proliferation and outcome is mixed, so no direction is claimed.")))
    S.append(dict(kind="visual", section="Pathway mechanism 3 of 4", image=VIS / "V06_p53.png",
        title="p53: a stress checkpoint held in check by MDM2",
        source="Literature: TCGA 2012 (R6). Measured: phase9 execution and S-7 results; phase8c/score_vs_annotation.csv. Post-hoc: reports/exploratory/x1, x4.",
        notes=("DNA damage activates ATM and ATR, which stabilise p53. p53 induces CDKN1A (p21) for cell-cycle arrest, BAX and PUMA for apoptosis, and GADD45 for repair. It also induces MDM2, which degrades p53: a negative feedback loop. "
               f"The score is the mean of {m['p53 pathway']['genes']} target genes, so it is not a mutation or protein-function assay. Primary hazard ratio {hr(m['p53 pathway']['hr'])}. "
               f"It has the strongest association with mismatch-repair status of the four scores (epsilon squared {m['p53 pathway']['mmr_e2']:.4f}, exploratory, n = {m['p53 pathway']['mmr_n']}), and it changes sign when some single batches are omitted. "
               "The direction and mechanism of the MMR association were not tested. "
               f"Post-hoc, on the same data: the p53 score is lower in TP53-mutant tumours ({x['tp53_p53_diff']:+.2f} SD), and the difference holds within mismatch-repair-proficient tumours ({x['tp53_pmmr_diff']:+.2f} SD, batch-stratified permutation p {x['tp53_pmmr_p']:.4f}). "
               f"This replicated in independent cohorts: GSE41258 {E['r1_41258']:+.2f} and TCGA-COAD {E['r1_tcga']:+.2f} (rank-biserial). "
               "That links the target-gene score to genotype as expected; it is not a measure of p53 protein function.")))
    S.append(dict(kind="visual", section="Pathway mechanism 4 of 4", image=VIS / "V07_emt.png",
        title="EMT: tumour-cell programme or stromal signal?",
        source="Literature: Guinney 2015 (R1); Calon 2015 (R2); Isella 2015 (R3). Measured: phase9 results. Post-hoc: reports/exploratory/x5, x6.",
        notes=("Two different biological sources produce the same mesenchymal transcripts. Tumour cells undergoing EMT lose E-cadherin and gain vimentin and fibronectin under SNAIL, ZEB and TWIST. Cancer-associated fibroblasts express those genes without any tumour-cell change. "
               "In bulk colorectal tumours, the mesenchymal-subtype signal is largely stromal according to Calon and Isella. "
               f"The EMT hazard ratio, {hr(m['EMT']['hr'])}, stays above 1 in every fit but is never distinguishable from 1. "
               "Cell composition was not measured. Post-hoc, it was estimated computationally from the same arrays: the EMT score correlates "
               f"{x['emt_fib']:.2f} with estimated fibroblast content, and {x['emtex_fib']:.2f} even after removing the {x['emt_overlap']} genes it shares with the stromal marker sets. "
               "That is consistent with a predominantly stromal source, as Calon and Isella reported, but it does not establish the source in either direction. "
               f"In registered follow-up, the same gene-disjoint relationship replicated in GSE41258, GSE143985 and TCGA-COAD (ρ {E['disjoint']['GSE41258'][0]:.2f}, {E['disjoint']['GSE143985'][0]:.2f}, {E['disjoint']['TCGA-COAD'][0]:.2f}). "
               "An elevated bulk EMT-associated score does not establish tumour-cell EMT.")))
    S.append(dict(kind="visual", section="What was measured", image=VIS / "V08_bulk_tissue.png",
        title="What a bulk tumour sample measures",
        source="Conceptual illustration. Estimates: reports/exploratory/x5 (ESTIMATE v1.0.13; MCP-counter v1.2.0), post-hoc.",
        notes=("Every array is a piece of bulk tumour: tumour epithelium, fibroblasts, immune cells and vessels, all homogenised into one mRNA pool. "
               "A pathway score is the average of its member genes' standardised values in that pool, so it cannot tell which cells expressed the genes. "
               "Cell composition was not measured. In post-hoc analyses it was estimated with two published reference methods, ESTIMATE and MCP-counter, which agree with each other (stromal versus fibroblast estimate ρ 0.93). "
               "These are relative scores computed from the same arrays, not measured cell fractions, and they carry measurement error. "
               f"All four pathway scores also associate with MMR status in outcome-free tests in {m['EMT']['mmr_n']} patients.")))
    S.append(dict(kind="side", section="Dataset and cohort", image=FIG / "Fig1_study_design_cohort.png",
        title=f"GSE39582: from {c['arrays']} arrays to {c['n']} analysed patients",
        left=[f"{c['arrays']} Affymetrix HG-U133 Plus 2.0 arrays, deposited as RMA + ComBat (Marisa 2013).",
              f"{c['n']} stage II–III patients: {c['events']} recurrences, {c['censored']} censored.",
              f"Endpoint: cause-specific hazard of recurrence; {c['deaths_censored']} deaths without recurrence are censored (informative).",
              "Recurrence-time unit undocumented: rank-based measures only.",
              f"CIT partition ({c['cit']['n']} patients, {c['cit']['events']} events) is internal, never external validation."],
        source="Figure 1 (reports/figures/final/Fig1_study_design_cohort.png). Censoring profile: appendix A2.",
        notes=(f"Walk the funnel: {c['arrays']} arrays, the stage II–III primary cohort of {c['S3_n']}, and the complete-case analysis population of {c['n']}. "
               "The model estimates the cause-specific hazard of recurrence. Deaths without recurrence are censored; appendix A2 shows those patients were older, less often stage III and less often treated, so censoring is informative. "
               "That is why no absolute risk or cumulative incidence is reported. The CIT partition shares processing history with the discovery set and is confounded with batch, so it is internal.")))
    S.append(dict(kind="visual", section="Preprocessing and features", image=VIS / "V10_preprocessing_panel.png",
        title=f"From {p['total']:,} probe sets to four pathway scores",
        source="figure_values.json (fig2.probes); config/pathways/m2_panel_4set.json and m2_panel_5set.json; PHASE6C §5.",
        notes=(f"The data were used as deposited, without re-normalisation. {p['unmapped']:,} probe sets were unmapped and {p['ambiguous']:,} ambiguous; the remaining {p['single_gene']:,} were collapsed to {p['genes']:,} genes by the median. "
               f"The four Hallmark sets contain {V['panel']['union4']} genes. Each score is the mean of member-gene z-scores, with standardisation re-estimated inside each training fold so held-out patients never shape the features. "
               "On the right: the ten candidate sets were highly collinear (condition number 77). One representative per correlation cluster was chosen without outcome data, which brings the four-set panel down to 2.95. "
               "The five-set panel, adding PI3K/AKT/mTOR, was fixed in advance as a sensitivity analysis only.")))
    S.append(dict(kind="visual", section="Model design", image=VIS / "V11_model_architecture.png",
        title="Model architecture: clinical benchmark and molecular extension",
        source="docs/PHASE6B_MODELING_PROTOCOL_DRAFT.md §3; λ from phase9_execution_results.json.",
        notes=("M1 is the clinical benchmark and was never re-tuned. M2 adds the four scores with a ridge penalty on the molecular block only. "
               "Leaving the clinical block unpenalised keeps the comparison scale-equivariant and makes M1 exactly nested in M2, so any gain is attributable to the pathway block. Shrinking only that block makes the gain in C conservative. "
               f"In M1, stage III versus II has hazard ratio {st[0]:.2f} ({st[1]:.2f}–{st[2]:.2f}); chemotherapy is confounded by indication and never read causally. "
               "The five-set model is shown separately because it is sensitivity only and was never promoted.")))
    S.append(dict(kind="visual", section="Validation design", image=VIS / "V12_nested_bootstrap.png",
        title="Validation: optimism-corrected nested bootstrap",
        source="docs/PHASE9_DELTAC_INTERVAL_PLAN.md §3; figure_values.json (fig8, figS1); verifier output.",
        notes=(f"Each of the {V['B']['out']} outer resamples refits both models and computes the apparent gain in C on that resample. "
               f"Inside it, {V['B']['in']} inner resamples estimate optimism: each inner refit is scored on its own resample and on the outer resample, with M1 and M2 paired on the same resamples. "
               "Corrected gain equals apparent gain minus the difference in optimism, and the 95 % interval is the percentile interval of the corrected gains. "
               "The five runs differ only in random seed: they are Monte Carlo replicates of one analysis, not independent experiments. "
               f"The run size met a pre-specified precision rule (largest endpoint Monte Carlo error {V['mc']['rho_max']*100:.1f} % of interval width, limit {V['mc']['threshold']*100:.0f} %).")))
    S.append(dict(kind="visual", section="Central result", image=VIS / "V13_central_result.png",
        title="Central result: no demonstrated gain in discrimination",
        subtitle=f"Optimism-corrected gain in C: +{V['theta']:.6f}. All five 95 % intervals include zero, so no reproducible improvement was demonstrated.",
        source="figure_values.json (fig7, fig9); reports/phase9/corrected_rerun_report/corrected_interval_results.json.",
        notes=(f"The pathway model's apparent lead was +{d_app:.4f} in C. Most of it was optimism: the corrected gain is +{V['theta']:.6f}. "
               f"The five run intervals are {runs}. All include zero. "
               "A gain of about 0.004 in C is very small in practical terms, and the interval cannot exclude a small loss or a modest gain. This does not establish that the true effect is zero. "
               f"Corrected calibration slopes of {V['slopes']['M1']:.3f} and {V['slopes']['M2_primary']:.3f} show both models overfit.")))
    S.append(dict(kind="visual", section="Pathway estimates", image=VIS / "V14_forest.png",
        title="Pathway hazard ratios: every primary interval includes 1",
        source="phase9_execution_results.json (M2 primary, M2 sensitivity); phase9_s7_results.json.",
        notes=("Hazard ratios are per training-fold standard deviation of each score, on a log axis with the no-effect line at 1. "
               "No primary-model interval excludes 1. EMT stays above 1 and E2F below 1 in every fit, which is compatible with parts of the literature but is not a demonstrated association. "
               "S-7 re-summarises the same arrays, so it is a robustness check on probe summarisation, not replication. "
               f"The five-set E2F estimate {hr(V['five_set_e2f'])} is shown separately because the model is sensitivity only; its interval excluding 1 coincided with collinearity.")))
    S.append(dict(kind="visual", section="Robustness", image=VIS / "V15_robustness.png",
        title="Robustness: batch sensitivity, subtype association, calibration",
        source="phase9_execution_results.json (batch_sensitivity); score_vs_annotation.csv; figure_values.json (fig9).",
        notes=("Left: each point is one leave-one-batch-out refit, plotted against how many patients were omitted. Wnt changes sign once and p53 twice; EMT and E2F keep their sign. Because batches range from 14 to 103 patients, batch identity and sample size are mixed. "
               "Batch location effects did not survive stratification by the ComBat covariate combat_mod, which matches recorded MMR in all 410 samples but has unresolved provenance, so it is not described as an MMR variable. "
               "Batch effects were characterised, not resolved. "
               "Top right: all four scores associate with MMR status (exploratory); these associations cannot be compared with hazard ratios. Bottom right: corrected calibration slopes below 1 indicate residual overfitting.")))
    S.append(dict(kind="visual", section="Contributions", image=VIS / "V16_project_established.png",
        title="What this project established",
        source="tools/pathway-audit/; docs/RESEARCH_EXTENSION_REGISTER.md + TCGA addendum. All limitations: docs/SCIENTIFIC_LIMITATION_RESOLUTION_REGISTER.md.",
        notes=("Four completed contributions, strongest first. The follow-up analyses behind blocks 2 and 3 are post-hoc relative to the primary study: they were registered before they ran and used no outcome data, but none of this is confirmatory. "
               "Block 1, the tool: pathway-audit is a Python package and command-line tool that checks data quality, gene drivers and coherence, tissue composition, molecular associations with covariate adjustment, a random-gene-set null and an optional survival screen. "
               "It now also has a single-cell source check. Its 13 package tests pass, continuous integration on GitHub Actions passes on Python 3.10 and 3.12, the built-in demo recovers all four planted truths, and the case studies cover five bulk colorectal cohorts, one of which, GSE103479, failed the quality gate, plus one single-cell cohort. The Docker image builds and runs locally, with the tests and the demo passing inside it; it is not published to a registry. "
               f"Block 2: with the {E['n_shared_fib']} shared genes removed, the bulk EMT score correlates {E['disjoint']['GSE39582'][0]:.2f}, {E['disjoint']['GSE41258'][0]:.2f}, {E['disjoint']['GSE143985'][0]:.2f} and {E['disjoint']['TCGA-COAD'][0]:.2f} with fibroblast markers in GSE39582, GSE41258, GSE143985 and TCGA-COAD. "
               f"Each is compared with its own random-gene-set baseline, whose upper 97.5 % bounds are {E['baseline']['GSE39582']['p97_5']:.2f}, {E['baseline']['GSE41258']['p97_5']:.2f}, {E['baseline']['GSE143985']['p97_5']:.2f} and {E['baseline']['TCGA-COAD']['p97_5']:.2f}; in TCGA random sets already reach a median of about {E['baseline']['TCGA-COAD']['median']:.2f}, so that cohort is less informative than its size suggests. "
               "In bulk tissue that is an association with estimated stromal content. The registered single-cell extension then measured it directly in three public cohorts: "
               f"EMT-set transcripts are {SC['emt_median']['stromal']:.1f} times more concentrated in stromal cells than expression-matched random genes, and {SC['emt_median']['epithelial']:.2f} times in epithelial cells, and a pseudo-bulk EMT score tracks the counted stromal fraction (rho {SC['s4_rho']:.2f}). "
               f"Within malignant epithelial cells the EMT genes were not more coherent than matched random genes in any of {SC['s3_n']} patients, so a tumour-cell EMT programme was not detected; rare or partial states are not ruled out. "
               f"Block 3: the p53 target-gene score is lower in TP53-mutant tumours. Discovery effect in GSE39582 {x['tp53_p53_diff']:+.2f} SD; replication as rank-biserial effects: GSE41258 {E['r1_41258']:+.2f}, TCGA-COAD {E['r1_tcga']:+.2f}, and {E['tcga_hyper']['R1']['excluding_hypermutated']['rank_biserial'][0]:+.2f} after excluding hypermutated tumours; GSE143985 was not significant, so it replicated in two of three external cohorts. "
               "This is an expression association; it does not measure p53 protein activity or pathway function. E2F–TP53 and p53–MSI did not replicate, and the TCGA BRAF links were explained by hypermutation; these are not findings. "
               f"Block 4: the pre-specified nested bootstrap gives a corrected gain in C of +{V['theta']:.6f}, and all five intervals include zero, so no reproducible improvement in recurrence discrimination was demonstrated (reports/phase9/corrected_rerun_report/corrected_interval_results.json). That does not mean the pathways have no biological role. "
               f"Still unresolved: rare or partial tumour-cell EMT states, direct pathway activation, external predictive validation ({x['ext_assessed']} public cohorts assessed, none eligible) and the recurrence-time unit, provisionally months. "
               f"Batch heterogeneity (frailty variance {x['theta']:.2f}), competing risks, MMR association and computational composition are partially addressed or characterised; the full status of every limitation is in the scientific limitation resolution register.")))
    S.append(dict(kind="visual", section="Decisions and audit trail", image=VIS / "V17_decisions_audit.png",
        title="Decisions and audit trail",
        source="docs/OWNER_DECISION_REGISTER.md; docs/PHASE9_CORRECTED_FINAL_REPORT.md §13–§15.",
        notes=("Each decision is recorded in the project's decision records, DC1 to DC7; describe the reasoning for each rather than claiming to have discovered every problem. "
               "DC1 fixed every analytic choice before outcomes; DC2 kept clinical terms unpenalised; DC3 kept the CIT partition internal; DC4 adjudicated the random-seed deviation and authorised a full re-run; DC5 completed the pre-specified 500 to 750 top-up; DC6 kept misleading results out of the conclusions; DC7 commissioned and reviewed the audits. "
               "The separate audit recomputed reported quantities from the recorded checkpoints; it did not refit models and is not an external replication.")))
    S.append(dict(kind="visual", section="From finding to solution", image=VIS / "V19_solution_tool.png",
        title="From finding to solution: replication and the pathway-audit tool",
        source="docs/RESEARCH_EXTENSION_REGISTER.md and TCGA addendum (registered before analysis); tools/pathway-audit/; reports/pathway_audit_case_study/.",
        notes=("Because the prediction answer was negative, the follow-up asked what these scores measure in bulk tissue, in analyses registered before they ran and using no outcome data. "
               f"The clearest result: the EMT score tracks fibroblast markers even with the {E['n_shared_fib']} shared genes removed, in GSE39582 and in three independent cohorts, including TCGA RNA-seq. "
               f"One public cohort, GSE103479, failed a positive-control check (control gene pairs correlated at only {E['qc']['GSE103479']:.2f}); its values appear misaligned to probe labels, so it was excluded rather than reported as non-replication. "
               "I then packaged the checks into pathway-audit, a tested Python command-line tool: a data-quality gate, gene drivers, composition, genotype effects adjusted for composition and MMR, and an optional survival screen. "
               f"On GSE39582 its survival screen gives a corrected gain of +{E['survival_screen']:.4f}, close to the frozen nested estimate. "
               "An independent re-check of TCGA showed that the apparent BRAF effects there are explained by hypermutation; with mutation burden as a covariate, the tool reports exactly that, while the TP53 link to the p53 score stays robust. "
               "Random gene sets also correlate with the fibroblast estimate in TCGA (about 0.64), so each cohort is compared with its own random baseline. Composition values are computational estimates, not measured cell fractions. "
               f"The newest check, single-cell sources, ran on GSE132465 tumour cells and reproduced the registered single-cell analysis: EMT-set transcripts mainly stromal ({SC['emt_median']['stromal']:.1f} times the expression-matched share). The tool has 13 tests and passes continuous integration.")))
    S.append(dict(kind="visual", section="Conclusions", image=VIS / "V18_conclusion_questions.png",
        title="Conclusion and the next scientific questions",
        source="docs/SCIENTIFIC_INTERPRETATION_AND_CONCLUSION.md; docs/FINAL_RESEARCH_CONTRIBUTION.md; docs/EXTERNAL_COHORT_ELIGIBILITY_AUDIT.md.",
        notes=("This slide pulls the project together rather than repeating the evidence. "
               f"The project asked whether four pathway scores improve recurrence prediction, and found no reproducible gain: corrected gain in C +{V['theta']:.6f}, with every interval including zero. That does not establish that the pathways have no biological role. "
               "It then asked what the scores measure in bulk tumour tissue. In registered, outcome-free follow-up, the EMT score tracked estimated stromal content in four bulk cohorts, and single-cell data from three more cohorts confirmed that its transcripts come mainly from stromal cells; the p53 score tracked TP53 status; other genotype links did not hold. "
               "It delivered pathway-audit so these checks can be run before a score is interpreted. "
               "The take-home is practical: audit bulk pathway scores for composition, genotype confounding and data quality before reading them as biology or adding them to prognostic models. "
               "Three questions follow, none of them performed here. First, whether rare or partial tumour-cell EMT states exist that single-cell RNA missed, which needs spatial transcriptomics or multiplexed protein imaging at the invasive front. "
               "Second, whether the pathways are actually active, which needs functional or protein-level readouts. "
               f"Third, whether any score adds recurrence information in a new cohort: that needs an independent stage II–III cohort with about {x['ext_needed']} recurrences; the largest of the {x['ext_assessed']} public cohorts assessed had {x['ext_max_events']}.")))
    S.append(dict(kind="text", section="Appendix A1", title="References and figure sources",
        left_title="References",
        left=["R1 Guinney J et al. Nat Med 2015;21:1350–1356. doi:10.1038/nm.3967",
              "R2 Calon A et al. Nat Genet 2015;47:320–329. doi:10.1038/ng.3225",
              "R3 Isella C et al. Nat Genet 2015;47:312–319. doi:10.1038/ng.3224",
              "R4 Anjomshoaa A et al. Br J Cancer 2009;101:822–828. doi:10.1038/sj.bjc.6605229",
              "R5 Luo ZW et al. BMC Cancer 2019;19:123. doi:10.1186/s12885-019-5324-y",
              "R6 TCGA Network. Nature 2012;487:330–337. doi:10.1038/nature11252",
              "R7 Liberzon A et al. Cell Syst 2015;1:417–425. doi:10.1016/j.cels.2015.12.004",
              "R8 Marisa L et al. PLoS Med 2013;10:e1001453. doi:10.1371/journal.pmed.1001453"],
        right_title="Figure sources",
        right=["Slide visuals: portfolio/diagrams/slides/ (original vector drawings by portfolio/src/make_slide_visuals.py, numbers read from frozen and registered post-hoc outputs).",
               "Replication data: GEO GSE41258, GSE143985, GSE103479; TCGA-COAD via UCSC Xena. Tool: tools/pathway-audit/.",
               "Figures 1, 2 and S2: reports/figures/final/ (generated from frozen outputs by scripts/49).",
               "Mechanism panels are textbook background, labelled as established mechanism; no external images are used.",
               "Full results: docs/SCIENTIFIC_RESULTS_AND_FIGURES.md and docs/PHASE9_CORRECTED_FINAL_REPORT.md."],
        source="References verified in docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md.",
        notes="Reference slide for questions. All eight references were checked against the published records. Every visual in the deck is an original drawing generated from code; numbers are read from the frozen outputs."))
    S.append(dict(kind="two", section="Appendix A2", images=[FIG / "Fig2_preprocessing_qc.png", FIG / "FigS2_informative_censoring.png"],
        title="Quality control and informative censoring",
        captions=["Figure 2. Probe mapping, per-array QC (29 soft flags, none excluded) and batch confounded with partition.",
                  f"Figure S2. The {c['deaths_censored']} patients censored at death without recurrence differ from the rest: censoring is informative."],
        source="Figures 2 and S2 (reports/figures/final/).",
        notes=("Backup for questions on data quality. Quality control used soft flags only: 29 arrays were flagged and none removed, because there was no independent evidence of array failure. Processing batch is confounded with the discovery/CIT partition. "
               "Figure S2 profiles the deaths without recurrence that are censored under the cause-specific endpoint: they were older, less often stage III and less often treated, so censoring is informative and no absolute risk is reported.")))
    return S


# ------------------------------------------------------------------ layout geometry (inches, top-left origin)
FOOT_Y = 7.05


def alt_texts(V):
    m = {x["label"]: x for x in V["pathways"]}
    hrs = "; ".join(f"{x['label']} {hr(x['hr'])}" for x in V["pathways"])
    return {
        "V01_cover_illustration.png": "Conceptual vector illustration: a normal colon crypt lined with epithelial cells on the left, and crowded tumour glands surrounded by fibroblasts and immune cells on the right. Not patient data.",
        "V02_problem_question.png": f"Four steps: the clinical problem, the research question, the optimism trap and the approach. Below: {V['cohort']['n']} patients, {V['cohort']['events']} recurrences, 4 pathway scores, {V['probes']['total']:,} probe sets. Endpoint is the cause-specific hazard; validation is internal only.",
        "V03_biology_map.png": "Mechanism map. Wnt inhibits the APC-AXIN-GSK3 beta destruction complex, which inhibits beta-catenin; beta-catenin with TCF/LEF activates MYC and cyclin D; cyclin D-CDK4/6 phosphorylates and inactivates RB; RB represses E2F; E2F activates G1/S genes and proliferation. DNA damage activates p53, which induces p21 (blocking CDK activity) and MDM2 (which inhibits p53). TGF-beta activates SNAIL, ZEB and TWIST, which repress E-cadherin and induce vimentin and fibronectin; fibroblasts also express mesenchymal genes.",
        "V04_wnt.png": f"Two states of Wnt signalling. With APC intact, the destruction complex phosphorylates beta-catenin, which is degraded. With APC lost or CTNNB1 mutated, beta-catenin accumulates and enters the nucleus with TCF/LEF to transcribe MYC, AXIN2 and cyclin D. Evidence panel: literature, measured Wnt score HR {hr(m['Wnt/β-catenin']['hr'])}, and what was not tested.",
        "V05_e2f.png": f"E2F mechanism: growth signals induce cyclin D-CDK4/6, which phosphorylates RB and releases E2F; p21 inhibits CDK4/6; E2F transcribes cyclin E, MCM2-7, PCNA and CDC6, driving DNA replication and proliferation. Evidence panel: mixed literature, measured E2F HR {hr(m['E2F targets']['hr'])}, five-set estimate shown as sensitivity only.",
        "V06_p53.png": f"p53 mechanism: DNA damage activates ATM/ATR, which stabilise p53; p53 induces MDM2, which degrades p53; p53 induces CDKN1A (p21), BAX and PUMA, and GADD45 for arrest, apoptosis and repair. Evidence panel: TP53 mutation frequency, measured p53 HR {hr(m['p53 pathway']['hr'])} and MMR epsilon squared {m['p53 pathway']['mmr_e2']:.4f}.",
        "V07_emt.png": f"Two sources of mesenchymal transcripts feed one bulk EMT score: tumour cells undergoing EMT (E-cadherin down, vimentin and fibronectin up) and fibroblasts in the stroma. Banner: an elevated bulk EMT-associated score does not establish tumour-cell EMT. Measured EMT HR {hr(m['EMT']['hr'])}.",
        "V08_bulk_tissue.png": "Conceptual illustration of bulk tumour tissue: tumour glands, fibroblasts, immune cells and vessels. Tumour cells not separated; fibroblast, immune and endothelial content not measured but estimated computationally (ESTIMATE, MCP-counter). Flow: mixed tissue to one array to one score.",
        "Fig1_study_design_cohort.png": f"Flow diagram from {V['cohort']['arrays']} GSE39582 arrays to the {V['cohort']['n']}-patient analysis population with {V['cohort']['events']} recurrences and {V['cohort']['censored']} censored, with the discovery and internal CIT partitions and the cause-specific endpoint.",
        "V10_preprocessing_panel.png": f"Funnel from {V['probes']['total']:,} probe sets to {V['probes']['single_gene']:,} single-gene probes, {V['probes']['genes']:,} genes, {V['panel']['union4']} pathway genes and 4 scores per patient. Bar chart of condition numbers: 77 for ten candidate sets, {V['panel']['cond5']:.2f} for the five-set and {V['panel']['cond4']:.2f} for the four-set panel.",
        "V11_model_architecture.png": f"Three models: M1 clinical block unpenalised; M2 adds a ridge-penalised block of four pathway scores (lambda {V['lambda']['primary']:.3f}); the five-set model adds PI3K/AKT/mTOR and is sensitivity analysis only.",
        "V12_nested_bootstrap.png": f"Nested bootstrap workflow: {V['B']['out']} outer event-stratified resamples, each refitting M1 and M2; {V['B']['in']} inner resamples per outer resample estimate paired optimism; corrected gain equals apparent gain minus optimism difference; percentile interval; repeated in 5 runs differing only in random seed. Verifier 60/60.",
        "V13_central_result.png": f"Left: apparent C {V['C']['M1']['apparent']:.4f} and {V['C']['M2_primary']['apparent']:.4f} shrink to corrected {V['C']['M1']['corrected']:.4f} and {V['C']['M2_primary']['corrected']:.4f}. Right: five corrected-gain intervals all include zero, with point estimate +{V['theta']:.6f}.",
        "V14_forest.png": f"Forest plot of pathway hazard ratios on a log axis. Primary: {hrs}. S-7 estimates beside each. All primary intervals include 1. Five-set E2F {hr(V['five_set_e2f'])} shown in a separate sensitivity-only box.",
        "V15_robustness.png": "Leave-one-batch-out log hazard ratios for each pathway against omitted batch size (EMT 0, E2F 0, Wnt 1, p53 2 sign changes); MMR association bars (p53 strongest); corrected calibration slopes 0.833, 0.862 and 0.836, all below 1.",
        "V16_project_established.png": "Four contribution blocks and a strip of open questions. 1, pathway-audit framework for bulk and single-cell data: six-step workflow (data-quality gate, gene drivers and coherence, tissue composition, molecular associations, random-set null, survival screen); 13 of 13 package tests pass; CI passes on Python 3.10 and 3.12; demo recovers 4 of 4 planted truths; 5 bulk and 1 single-cell case-study cohorts; Docker image builds and runs locally. 2, EMT score tracks stroma in 4 bulk and 3 single-cell cohorts: dot plot of gene-disjoint Spearman rho 0.94, 0.95, 0.94 and 0.98 above grey cohort-specific random-gene-set bands; bars of single-cell EMT-set transcript enrichment: stromal 5.08, endothelial 1.59, immune 0.92, epithelial 0.45 times the expression-matched share; tumour-cell EMT programme not detected, not ruled out. 3, molecular anchoring: p53 score lower in TP53-mutant tumours; rank-biserial effects GSE41258 -0.22, GSE143985 -0.06 not significant, TCGA-COAD -0.42, TCGA non-hypermutated -0.37; not a measure of p53 protein activity; E2F-TP53, p53-MSI and TCGA BRAF links not established. 4, predictive evaluation: corrected delta C +0.004316 with five intervals all including zero. Unresolved: rare or partial tumour-cell EMT, direct pathway activation, external predictive validation, recurrence-time unit.",
        "V17_decisions_audit.png": "Nine-stage project timeline with decisions DC1 to DC7 attached to their stages, and verification facts: verifier 60 of 60, independent recomputation within 2.2e-16 for 138 quantities, 1,097 tests.",
        "V19_solution_tool.png": "Left: bar chart of gene-disjoint Spearman correlation between the EMT score and fibroblast markers, 0.94 in GSE39582, 0.95 in GSE41258, 0.94 in GSE143985 and 0.98 in TCGA-COAD; GSE103479 failed the quality gate. Right: the six checks of the pathway-audit tool, including the single-cell source check, and its case-study results in five bulk cohorts and one single-cell cohort.",
        "V18_conclusion_questions.png": "Left: a five-step flow. Asked whether four pathway scores improve recurrence prediction; found no reproducible gain, corrected delta C +0.004316 with every interval including zero; asked what the scores measure in bulk tissue; found the EMT score tracks stroma in 4 bulk and 3 single-cell cohorts and the p53 score tracks TP53 status; delivered pathway-audit. Take-home box: audit bulk pathway scores for composition, genotype confounding and data quality before interpreting them. Right: three next scientific questions, not performed, each with how and why: rare or partial tumour-cell EMT states missed by single-cell RNA; pathway activity; predictive value in an independent cohort with about 337 recurrences.",
        "Fig2_preprocessing_qc.png": "Four panels: probe mapping to 21,815 genes; per-array medians by batch; 29 soft QC flags with none excluded; batch composition confounded with partition.",
        "FigS2_informative_censoring.png": "Profile of patients censored at death without recurrence compared with the others, showing differences in age, stage and treatment: censoring is informative.",
    }


def img_box(path, x, y, w, h):
    im = mpimg.imread(path)
    ar = im.shape[1] / im.shape[0]
    if w / h > ar:
        iw, ih = h * ar, h
    else:
        iw, ih = w, w / ar
    return im, x + (w - iw) / 2, y + (h - ih) / 2, iw, ih


def layout(s):
    """Return (text blocks, image boxes) in inches."""
    k = s["kind"]
    if k == "visual":
        top = 1.42 if s.get("subtitle") else 1.15
        return [], [(0.55, top, 12.2, FOOT_Y - 0.1 - top)]
    if k == "side":
        return [("left", 0.55, 1.3, 4.1, 5.5)], [(4.8, 1.1, 8.2, 5.85)]
    if k == "two":
        return [], [(0.55, 1.15, 9.9, 2.8), (0.55, 4.1, 9.9, 2.8)]
    if k == "cover":
        return [], [(7.6, 0.85, 5.3, 5.3)]
    return [("left", 0.55, 1.3, 6.4, 5.5), ("right", 7.2, 1.3, 5.6, 5.5)], []


def images(s):
    return s.get("images") or ([s["image"]] if "image" in s else [])


# ------------------------------------------------------------------ PDF renderer
def wrap(t, w, fs):
    return textwrap.fill(t, max(10, int(w * 72 / (fs * 0.5))))


def pdf_slide(pdf, s, n, total, alt):
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    blocks, boxes = layout(s)
    if s["kind"] == "cover":
        ax.text(0.7, 1.1, s["section"].upper(), fontsize=12, color=PRIMARY, fontweight="bold", va="top")
        t = wrap(s["title"], 6.4, 32)
        ax.text(0.7, 1.55, t, fontsize=32, fontweight="bold", va="top", linespacing=1.15)
        y = 1.55 + (t.count("\n") + 1) * 32 * 1.15 / 72 + 0.35
        st = wrap(s["subtitle"], 6.4, 15)
        ax.text(0.7, y, st, fontsize=15, color=MUTED, va="top", linespacing=1.3)
        ax.plot([0.7, 2.2], [5.45, 5.45], color=PRIMARY, lw=3)
        ax.text(0.7, 5.7, AUTHOR, fontsize=18, fontweight="bold", va="top")
        ax.text(0.7, 6.1, AFFIL, fontsize=14, color=MUTED, va="top")
        ax.text(0.7, 6.45, DATE, fontsize=14, color=MUTED, va="top")
    else:
        ax.text(0.55, 0.22, s["section"].upper(), fontsize=10, color=PRIMARY, fontweight="bold", va="top")
        ax.text(0.55, 0.45, s["title"], fontsize=24, fontweight="bold", va="top")
        if s.get("subtitle"):
            ax.text(0.55, 1.0, wrap(s["subtitle"], 12.2, 13), fontsize=13, color=MUTED, va="top")
    for side, x, y, w, h in blocks:
        head = s.get(f"{side}_title")
        yy = y
        if head:
            ax.text(x, yy, head, fontsize=15, fontweight="bold", va="top")
            yy += 0.45
        fs = 11.5 if s["kind"] == "text" else 13
        for it in s[side]:
            t = wrap(it, w - 0.3, fs)
            ax.text(x, yy, "•", fontsize=fs, color=PRIMARY, va="top")
            ax.text(x + 0.25, yy, t, fontsize=fs, va="top", linespacing=1.3)
            yy += (t.count("\n") + 1) * fs * 1.3 / 72 + 0.18
    for i, (path, b) in enumerate(zip(images(s), boxes)):
        im, ix, iy, iw, ih = img_box(path, *b)
        ax.imshow(im, extent=(ix, ix + iw, iy + ih, iy), interpolation="lanczos")
        if s["kind"] == "two":
            ax.text(10.75, b[1] + 0.2, wrap(s["captions"][i], 2.1, 12), fontsize=12, color=INK, va="top", linespacing=1.3)
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.text(0.55, FOOT_Y + 0.12, s["source"], fontsize=8.5, color=MUTED, va="center")
    ax.text(W - 0.45, FOOT_Y + 0.12, f"{n} / {total}", fontsize=9, color=MUTED, va="center", ha="right")
    pdf.savefig(fig)
    if PREVIEW:
        fig.savefig(Path(PREVIEW) / f"slide{n:02d}.png", dpi=60)
    plt.close(fig)


# ------------------------------------------------------------------ PPTX renderer
def rgb(h):
    return RGBColor.from_string(h.lstrip("#"))


def add_text(slide, x, y, w, h, paras, size, bold=False, color=INK, bullets=False):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    for i, t in enumerate(paras):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        para.text = ("• " if bullets else "") + t
        para.space_after = Pt(6 if bullets else 0)
        for r_ in para.runs:
            r_.font.size = Pt(size)
            r_.font.bold = bold
            r_.font.name = "Arial"
            r_.font.color.rgb = rgb(color)
    return tb


def style_title(ttl, x, y, w, h, size):
    ttl.left, ttl.top, ttl.width, ttl.height = Inches(x), Inches(y), Inches(w), Inches(h)
    tf = ttl.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = 1  # top
    for para in tf.paragraphs:
        para.alignment = 1  # left
        for r_ in para.runs:
            r_.font.size = Pt(size)
            r_.font.bold = True
            r_.font.name = "Arial"
            r_.font.color.rgb = rgb(INK)


def pptx_slide(prs, s, n, total, alt):
    sl = prs.slides.add_slide(prs.slide_layouts[5])  # "Title Only": a real, visible title placeholder
    blocks, boxes = layout(s)
    ttl = sl.shapes.title
    ttl.text = s["title"]
    if s["kind"] == "cover":
        add_text(sl, 0.7, 0.95, 6.4, 0.4, [s["section"].upper()], 12, True, PRIMARY)
        style_title(ttl, 0.7, 1.4, 6.6, 2.4, 32)
        add_text(sl, 0.7, 3.75, 6.4, 1.4, [s["subtitle"]], 15, color=MUTED)
        bar = sl.shapes.add_shape(1, Inches(0.7), Inches(5.43), Inches(1.5), Inches(0.04))
        bar.fill.solid()
        bar.fill.fore_color.rgb = rgb(PRIMARY)
        bar.line.fill.background()
        add_text(sl, 0.7, 5.6, 6.4, 0.4, [AUTHOR], 18, True)
        add_text(sl, 0.7, 6.0, 6.4, 0.75, [AFFIL, DATE], 14, color=MUTED)
    else:
        add_text(sl, 0.55, 0.12, 12.2, 0.3, [s["section"].upper()], 10, True, PRIMARY)
        style_title(ttl, 0.55, 0.38, 12.2, 0.6, 24)
        if s.get("subtitle"):
            add_text(sl, 0.55, 0.93, 12.2, 0.45, [s["subtitle"]], 13, color=MUTED)
    for side, x, y, w, h in blocks:
        yy = y
        head = s.get(f"{side}_title")
        if head:
            add_text(sl, x, yy - 0.05, w, 0.4, [head], 15, True)
            yy += 0.45
        add_text(sl, x, yy, w, h - (yy - y), s[side], 11.5 if s["kind"] == "text" else 13, bullets=True)
    for i, (path, b) in enumerate(zip(images(s), boxes)):
        _, ix, iy, iw, ih = img_box(path, *b)
        pic = sl.shapes.add_picture(str(path), Inches(ix), Inches(iy), Inches(iw), Inches(ih))
        pic._element.nvPicPr.cNvPr.set("descr", alt[Path(path).name])
        pic.name = Path(path).stem
        if s["kind"] == "two":
            add_text(sl, 10.7, b[1] + 0.1, 2.15, 2.5, [s["captions"][i]], 12)
    add_text(sl, 0.55, FOOT_Y, 11.0, 0.3, [s["source"]], 8.5, color=MUTED)
    add_text(sl, W - 1.45, FOOT_Y, 1.0, 0.3, [f"{n} / {total}"], 9, color=MUTED)
    sl.notes_slide.notes_text_frame.text = s["notes"]


def main():
    V = values.load()
    S = build_slides(V)
    alt = alt_texts(V)
    OUT.mkdir(parents=True, exist_ok=True)
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    total = len(S)
    with PdfPages(OUT / "GSE39582_recurrence_portfolio.pdf", metadata={"Title": "GSE39582 colon-cancer recurrence study", "Author": AUTHOR, "CreationDate": None}) as pdf:
        for i, s in enumerate(S, 1):
            pdf_slide(pdf, s, i, total, alt)
            pptx_slide(prs, s, i, total, alt)
    prs.core_properties.title = "GSE39582 colon-cancer recurrence study"
    prs.core_properties.author = AUTHOR
    prs.save(OUT / "GSE39582_recurrence_portfolio.pptx")
    lines = ["# Speaker notes: GSE39582 colon-cancer recurrence study", "",
             "Generated from the same slide specification as the presentation (`portfolio/src/make_deck.py`). "
             "Slides 1–19 are the main talk; A1 and A2 are appendix slides.", ""]
    for i, s in enumerate(S, 1):
        lines += [f"## Slide {i}. {s['title']}", "", f"*Section:* {s['section']}  ", f"*Source:* {s['source']}", "", s["notes"], ""]
    (R / "portfolio" / "SPEAKER_NOTES.md").write_text("\n".join(lines))
    print("slides:", total)


if __name__ == "__main__":
    main()
