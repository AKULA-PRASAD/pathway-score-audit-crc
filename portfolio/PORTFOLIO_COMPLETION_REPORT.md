# Portfolio completion report

**Date:** 8 November 2026 (final).

**Scope:** readiness of the GSE39582 recurrence study as an entry-level US computational-biology and bioinformatics/ML portfolio.

**Integrity:**
- No frozen analysis output, configuration, checkpoint, protocol or scientific figure was changed.
- No model or simulation was re-run.
- Nothing was staged, committed, pushed, released or made public.

Status key (one status per requirement):

| Status | Meaning |
|---|---|
| **VERIFIED** | Checked automatically against frozen sources, or by an executed command, with a recorded result |
| **VISUALLY INSPECTED** | Rendered output was looked at |
| **STRUCTURAL ONLY** | File structure checked; not rendered in its target application |
| **BLOCKED** | Attempted, but could not be performed in this environment |
| **NOT PERFORMED** | Not attempted |
| **OWNER DECISION** | Outside the scope of verification; requires the owner |

## 1. Final acceptance matrix

### 1.1 Scientific integrity

| ID | Requirement | Status | Evidence |
|---|---|---|---|
| S1 | Cohort, exclusions, endpoint and censoring | VERIFIED | QC checks against `reports/figures/final/figure_values.json` and the execution JSON: 452 / 136 / 316; 46 censored deaths; CIT 99 / 29 |
| S2 | Probe mapping, gene sets and scoring | VERIFIED | QC against `fig2.probes` and `config/pathways/m2_panel_*.json`: 54,675 → 43,163 → 21,815; 198 / 195 / 41 / 194; union 607 |
| S3 | Model specifications | VERIFIED | QC against frozen M1 terms; ridge on molecular block only |
| S4 | Corrected ΔC, five intervals, calibration slopes | VERIFIED | QC against `corrected_interval_results.json` and `figure_values.json`: +0.004316; all five intervals include zero; 0.833 / 0.862 / 0.836 |
| S5 | Primary, S-7 and five-set hazard ratios; batch-stratified and leave-one-batch-out results | VERIFIED | QC against execution and S-7 JSON; per-batch values from `batch_sensitivity` |
| S6 | MMR associations and `combat_mod` caveat | VERIFIED | QC against `reports/phase8c/score_vs_annotation.csv` and Phase 8B/8C text |
| S7 | External-validation feasibility | VERIFIED | QC against Phase 4/5 text: 48 / 18 / about 338; 129 overlapping patient pairs |
| S8 | Biological interpretation and citations | VERIFIED | References R1–R8 verified (`docs/BIOLOGICAL_EVIDENCE_AND_CITATION_AUDIT.md`); QC forbidden-claim scan |
| S9 | All ten unresolved limitations disclosed | VERIFIED | QC required-statement checks in slides, README, research report and web page |
| S10 | Post-hoc limitation investigation (X1–X9) | VERIFIED | Registered plan (hash-checked); 22/22 validation checks; independent re-implementation matched all key quantities after one estimator correction; `docs/SCIENTIFIC_LIMITATION_RESOLUTION_REGISTER.md`. No limitation classified RESOLVED; primary conclusion unchanged |

### 1.2 Interactive portfolio

| ID | Requirement | Status | Evidence |
|---|---|---|---|
| W1 | Published source retrieved and compared with the repository source | VERIFIED | Round 7: live version confirmed equal to the pre-edit source before republishing; updated source published as version 9 at the same URL. Earlier: published file downloaded from the artifact. Identical to `portfolio/web/index.html` except the platform wrapper (3 lines) and one platform meta tag |
| W2 | Automated interaction, tooltip, keyboard and overflow tests on the identical source | VERIFIED | `portfolio/src/web_check.py`, run locally in headless Chrome on the retrieved published file: 74/74 at 1280, 820 and 500 px; no script errors |
| W3 | Every number on the page | VERIFIED | All 137 decimal tokens in the retrieved source matched to frozen outputs or documented quotations. The rest are styling constants, DOIs and version strings |
| W4 | Layout at desktop, tablet and narrow widths | VISUALLY INSPECTED | Full-page screenshots of the source reviewed; earlier layout defects fixed |
| W5 | Colour contrast in both themes | VERIFIED | Computed WCAG ratios; all checked text and status pairs ≥ 4.5:1 |
| W6 | Interaction testing of the signed-in published page itself | NOT PERFORMED | The published page requires a signed-in session; W2 covers the identical source |
| W7 | Real screen-reader use (VoiceOver or NVDA) | NOT PERFORMED | Semantics were checked programmatically only (within W2) |
| W8 | Physical phones narrower than 500 px | NOT PERFORMED | Headless Chrome's minimum tested width is 500 px |

### 1.3 Presentation and accessibility

| ID | Requirement | Status | Evidence |
|---|---|---|---|
| P1 | Slide count | VERIFIED | 21 slides (19 main plus 2 appendix); one slide created per specification entry (QC) |
| P2 | One meaningful, visible title per slide | VERIFIED | 21/21 title placeholders on the slide, each carrying the slide's own title (QC) |
| P3 | Descriptive alt text | VERIFIED | 21/21 pictures, written from verified values (QC) |
| P4 | Speaker notes | VERIFIED | 21/21 slides; notes match the slide specification in order (QC) |
| P5 | PDF content and layout | VISUALLY INSPECTED | All 20 pages reviewed after the redesign; slides 1, 6, 7, 16, 18 and 19 re-inspected after the round-7 update. Residual: Figure 1 (slide 9) and Figure 2 (A2) keep small labels; slide 16 body text is about 9–10 pt at slide size |
| P6 | PPTX rendered in PowerPoint or Keynote | BLOCKED | No owner inspection result recorded. The round-7 rebuild (21 slides; slides 1, 6, 7, 16, 18, 19 changed or new) voids any earlier check. Manual checklist in section 3 |
| P7 | PPTX structure as a whole | STRUCTURAL ONLY | P1–P4 inspected in the file; not rendered (see P6) |
| P8 | PDF accessibility tagging | NOT PERFORMED | The PDF is untagged; the PPTX is the accessible version |

### 1.4 Documentation and reproducibility

| ID | Requirement | Status | Evidence |
|---|---|---|---|
| D1 | Research report | VERIFIED | `portfolio/RESEARCH_REPORT.md`; key numbers checked by QC |
| D2 | Recruiter summary | VERIFIED | `portfolio/README.md`, "In plain terms" |
| D3 | Interview preparation | VERIFIED | `portfolio/INTERVIEW_PREPARATION.md`; key numbers checked by QC |
| D4 | Skills evidence and role fit | VERIFIED | `portfolio/SKILLS_EVIDENCE_AND_ROLE_FIT.md` |
| D5 | Technical README and reproducibility guide | VERIFIED | Repository `README.md`; `portfolio/REPRODUCIBILITY_GUIDE.md` with each command marked TESTED or DOCUMENTED, NOT RE-RUN |
| D6 | Tested reproduction commands | VERIFIED | Raw-data checksums OK; archive manifest 52/52; independent audit reproduced in about 7 s; pytest 1,097 passed; portfolio rebuilt from a fresh environment (diagrams byte-identical, QC passed). Slide visuals added in round 5 rebuild byte-identical in the existing environment; not yet rebuilt from a fresh one |
| D7 | Fresh GEO download and full analysis re-run | NOT PERFORMED | Network download and many compute hours; documented in the guide but not executed |
| D8 | Scientific figures and captions | VERIFIED | `reports/figures/final/` (11 figures, captions, values, sources), unchanged; diagrams D1–D5; slide visuals V01–V18 (labels and frozen values checked in their SVG text by QC) |
| D9 | Asset inventory | VERIFIED | `portfolio/ASSET_INVENTORY.md`: paths, sizes and hashes recomputed. The generator refuses duplicates, missing deliverables, unlabelled files and missing sources; re-running it gives identical output |
| D10 | No tool-attribution terms in portfolio-facing text | VERIFIED | QC scan of slides, notes, README, documents and web page. The QC script itself contains the search patterns |

### 1.5 Repository and release integrity

| ID | Requirement | Status | Evidence |
|---|---|---|---|
| R1 | Frozen outputs, figures and archive unchanged | VERIFIED | 345 protected files and 25 figure files match recorded hashes; archive manifest 52/52 |
| R2 | GitHub Release unchanged | VERIFIED | Tag `phase9-corrected-campaign-2026-11-08`; asset digest `sha256:0812b3ba…2aa2f` |
| R3 | Repository state unchanged | VERIFIED | HEAD `bf3c456`, level with the remote; nothing staged; visibility PRIVATE |
| R4 | Code licence | DONE | MIT chosen by the owner (8 Nov 2026); `LICENSE` at the root and in `tools/pathway-audit/` |
| R5 | Per-sample data governance before any public release | OWNER DECISION | Open items in `docs/PUBLIC_RELEASE_REPRODUCIBILITY_AUDIT.md` |

**Final readiness: PORTFOLIO READY WITH DISCLOSED LIMITATIONS.** Full readiness is withheld because P6 is blocked, and W6, W7, W8, P8 and D7 were not performed.

## 2. Owner actions
1. Open the PPTX in PowerPoint or Keynote and work through the checklist in section 3 (resolves P6).
2. Share the web page from its Share menu before sending the link.
3. Decide R4 and R5 before making the repository public.

## 3. Manual PPTX checklist
Open `portfolio/presentation/GSE39582_recurrence_portfolio.pptx` beside `portfolio/presentation/GSE39582_recurrence_portfolio.pdf`.

**Every slide:**
- the title is visible at the top and not overlapped;
- no text runs past its box;
- Arial is shown, or a substitution was accepted;
- the source line is bottom left and "n / 21" is bottom right;
- the layout matches the PDF page.

**Slide by slide:**

| Slide | What to check |
|---|---|
| 1 Cover | Title wraps to at most three lines inside the left half; three-line subtitle (mentions the replication and tool), "Prasad Akula", "MS in Bioinformatics, Northeastern University" and "November 2026" visible; illustration on the right with its "Conceptual illustration" line |
| 2 Clinical problem | Four boxes fully inside the slide; four numbers and the endpoint line visible |
| 3 Biological axes | Subtitle on one line; arrowheads and flat inhibition bars distinguishable; "nucleus" label not overlapping arrows |
| 4–7 Mechanisms | Each shows the four evidence labels; pathway text panel does not overlap the diagram; on slide 7 the banner sentence is fully visible |
| 8 Bulk tissue | "NOT SEPARATED" and three "NOT MEASURED" statuses readable; "Conceptual illustration" line visible |
| 9 Cohort | Five bullets do not overlap Figure 1; note whether Figure 1's small labels are legible on the projector |
| 10 Preprocessing | Funnel numbers and the three condition-number bars (77, 6.95, 2.95) readable |
| 11 Model architecture | Three model boxes; λ = 0.351 visible; three bullets below |
| 12 Nested bootstrap | "750", "100", "paired", "5 runs" and "differ only in random seed" visible; both footnote lines inside the slide |
| 13 Central result | Subtitle shows +0.004316; θ̂ shows the hat over theta; conclusion sentence inside the blue box |
| 14 Forest plot | Eight interval labels do not overlap the dashed sensitivity box |
| 15 Robustness | Four batch panels, MMR bars and calibration bars; footnote on two lines |
| 16 What this project established | Four contribution blocks (tool, EMT replication, p53 anchoring, predictive evaluation) inside the slide; six workflow chips not clipped; dot-plot values 0.94 / 0.95 / 0.94 / 0.98 and p53 effects readable; "What remains unresolved" strip with four items visible above the footer |
| 17 Decisions | Four decision boxes do not overlap; DC1–DC7 all visible |
| 18 From finding to solution (new) | Four bars with values 0.94 / 0.95 / 0.94 / 0.98; GSE103479 "failed quality gate" line; five tool steps and four case-study bullets inside the slide |
| 19 Conclusion and next questions | Five-step Asked / Found / Delivered flow and the take-home box on the left; three question cards (Q1–Q3) with How / Why lines and the "NOT TESTED / HYPOTHESIS" label on the right |
| A1 References | Eight references and four source bullets readable |
| A2 QC and censoring | Both figures and their captions visible |

**Speaker notes:** with View → Notes, each slide matches the same-numbered section of `portfolio/SPEAKER_NOTES.md`.

**Accessibility checker:** Review → Check Accessibility should report no missing-title or missing-alt-text errors. Record any reading-order warnings.

## 4. Acceptance history (superseded; for the record only)
- **Round 1, completion audit.** Missing deliverables created (research report, interview preparation, reproducibility guide, skills and role fit, web source copy, web test runner). Web contrast, chart-button semantics and two wording gaps corrected. QC extended to cover all deliverables.
- **Round 2, acceptance.** The published source was compared with the repository and every number on it audited. Missing PPTX slide titles and file-name alt text were found and corrected; the titles on diagram slides were first placed above the visible area.
- **Round 3, closure.** Off-slide titles replaced by visible section titles; PPTX and PDF rebuilt; PDF re-inspected; the inventory generator hardened against duplicates and omissions.
- **Round 4, documentation consistency.** This report restructured into five groups with unique IDs, and the live-page wording separated into W1, W2 and W6.
- **Round 9, contribution-led slides.** Slide 16 redesigned as "What this project established": four contribution blocks (pathway-audit framework, replicated EMT–stroma association with cohort-specific baselines, TP53–p53 anchoring, pre-specified predictive evaluation) and a compact "What remains unresolved" strip; full limitation detail moved to the speaker notes and the limitation register. Slide 19 rebuilt as "Conclusion and the next scientific questions" so it no longer repeats slide 16. No scientific value changed; frozen outputs untouched. P6 remains blocked.
- **Round 8, portfolio hardening.** Independent re-check of TCGA (all numbers matched; R3/R4 reclassified as confounded by hypermutation; cohort-specific random-gene-set baselines replace a single-cohort value); tool gains a built-in demo, self-contained HTML reports, user-supplied composition and burden adjustment (tests 11/11); repository front page, CITATION.cff and evidence-backed resume text; slide 18 and web page corrected and republished (version 10).
- **Round 7, solution update.** Registered research extensions E1–E3 with TCGA-COAD addendum; `pathway-audit` tool (tests 9/9) and five-cohort case study; new slide 18 and updated slides 1, 6, 7, 16, 19; web page gains a Replication & tool section, updated evidence matrix and roadmap (web tests 74/74 at three widths) and was republished.
- **Round 6, limitation investigation.** Registered post-hoc analyses X1–X9 (molecular anchoring, competing risks, batch frailty, MMR provenance, ESTIMATE / MCP-counter composition, composition-adjusted models, external-cohort audit, time-unit search), independent review, four new documents and register; slides 4, 6, 7, 8, 16 and 18 updated; QC extended to 380 checks. Frozen outputs unchanged.
- **Round 5, presentation redesign.** The deck was rebuilt as 18 main and 2 appendix slides: a cover, mechanism-level slides for Wnt/β-catenin, E2F, p53 and EMT with evidence labels, and new vector visuals for the cohort funnel, model architecture, nested bootstrap, central result, forest plot, robustness and evidence matrix (`portfolio/src/make_slide_visuals.py`). QC extended to the new visuals; all 20 PDF pages inspected. No scientific value changed. P6 remains blocked.
