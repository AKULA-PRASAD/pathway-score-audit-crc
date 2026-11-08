# Biological evidence and citation audit

**Date:** 8 November 2026. **Type:** documentation and evidence audit only. No model was fitted and no analysis was run. The change set in section 5 was approved by the owner and applied on 8 November 2026 (wording, citations and one figure caption only).

**Files inspected:**
- `docs/SCIENTIFIC_RESULTS_AND_FIGURES.md` (SRF)
- `docs/PHASE9_CORRECTED_FINAL_REPORT.md` (FR)
- `reports/figures/final/FIGURE_CAPTIONS.md` (CAP)
- `README.md`
- for consistency: `docs/SCIENTIFIC_RESULTS_FIGURE_AUDIT.md` (SFA)

**Reference verification method:**
- Each reference was retrieved on 8 November 2026 from Europe PMC (REST API: title, authors, journal, volume, pages, DOI, PMID, PMCID and abstract).
- Where the abstract did not contain the specific statement, the open-access full text was retrieved (NCBI BioC, PMC) and the exact sentence located.
- Nothing below is cited from memory. One candidate (Samuels et al., *Science* 2004, doi:10.1126/science.1096502) was verified for metadata only. Its abstract is not available through the API, so it is **not used** to support any number.

---

## 1. Verified references

| # | Reference | Identifiers | What was verified (exact support) |
|---|---|---|---|
| R1 | Guinney J, Dienstmann R, Wang X, et al. The consensus molecular subtypes of colorectal cancer. *Nat Med* 2015;21(11):1350–1356 | doi:[10.1038/nm.3967](https://doi.org/10.1038/nm.3967); PMID 26457759; PMC4636487 | Abstract: "CMS4 (mesenchymal, 23%), prominent transforming growth factor-β activation, stromal invasion and angiogenesis"; "CMS2 (canonical, 37%), epithelial, marked WNT and MYC signaling activation". Full text: "CMS4 tumors displayed worse overall survival … and relapse-free survival". Full text: **GSE39582 is one of the 18 data sets** used to derive the CMS |
| R2 | Calon A, Lonardo E, Berenguer-Llergo A, et al. Stromal gene expression defines poor-prognosis subtypes in colorectal cancer. *Nat Genet* 2015;47(4):320–329 | doi:[10.1038/ng.3225](https://doi.org/10.1038/ng.3225); PMID 25706628 | Abstract: the predictive power of poor-prognosis classifications "arises from genes expressed by stromal cells rather than epithelial tumor cells"; stromal markers "associate robustly with disease relapse" |
| R3 | Isella C, Terrasi A, Bellomo SE, et al. Stromal contribution to the colorectal cancer transcriptome. *Nat Genet* 2015;47(4):312–319 | doi:[10.1038/ng.3224](https://doi.org/10.1038/ng.3224); PMID 25706627 | Abstract: mRNA levels of stem/serrated/mesenchymal-subtype genes "were mostly due to stromal expression"; high CAF-signature expression "was associated with poor prognosis in untreated CRC" |
| R4 | Anjomshoaa A, Nasri S, Humar B, et al. Slow proliferation as a biological feature of colorectal cancer metastasis. *Br J Cancer* 2009;101(5):822–828 | doi:[10.1038/sj.bjc.6605229](https://doi.org/10.1038/sj.bjc.6605229); PMID 19654572; PMC2736847 | Abstract (73 primary CRCs, 27 liver metastases, gene proliferation signature): "relapsed primaries or those with established metastases had significantly lower proliferative activity than CRCs that were non-metastatic and did not relapse". Small study |
| R5 | Luo ZW, Zhu MG, Zhang ZQ, et al. Increased expression of Ki-67 is a poor prognostic marker for colorectal cancer patients: a meta analysis. *BMC Cancer* 2019;19:123 | doi:[10.1186/s12885-019-5324-y](https://doi.org/10.1186/s12885-019-5324-y); PMID 30727976 | Abstract (34 studies, 6,180 patients): high Ki-67 associated with poorer OS (HR 1.54) and DFS (HR 1.43). It **contradicts** a simple "high proliferation = better outcome" statement. Ki-67 protein staining is not an E2F mRNA score |
| R6 | Cancer Genome Atlas Network. Comprehensive molecular characterization of human colon and rectal cancer. *Nature* 2012;487(7407):330–337 | doi:[10.1038/nature11252](https://doi.org/10.1038/nature11252); PMID 22810696; PMC3401966 | Full text: "the WNT signalling pathway was altered in 93% of all tumours, including biallelic inactivation of APC … or activating mutations of CTNNB1 in ∼80% of cases"; TP53 mutated in "60 versus 20%" of non-hypermutated vs hypermutated tumours; "mutually exclusive mutations in PIK3R1 and PIK3CA as well as deletions in PTEN in 2%, 15% and 4% of non-hypermutated tumours" |
| R7 | Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst* 2015;1(6):417–425 | doi:[10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004); PMID 26771021; PMC4707969 | Abstract: each hallmark is a refined gene set "that conveys a specific biological state or process and displays coherent expression". It defines the gene sets; it makes no prognostic claim |
| R8 | Marisa L, de Reyniès A, Duval A, et al. Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value. *PLoS Med* 2013;10(5):e1001453 | doi:[10.1371/journal.pmed.1001453](https://doi.org/10.1371/journal.pmed.1001453); PMID 23700391; PMC3660251 | Source of GSE39582: 566 samples meeting RNA quality, a discovery subset of 443, six subtypes |

Stable record links: `https://europepmc.org/article/MED/<PMID>`.

---

## 2. Claim-by-claim evidence table

Support codes: **D** = this project's data; **L** = external literature (verified, section 1); **H** = hypothesis or inference; **M** = methodological statement (true by definition of the measurement).

| ID | Exact claim (location) | Support | Evidence | Verdict | Narrower wording needed? |
|---|---|---|---|---|---|
| B1 | "Mesenchymal and stromal expression programmes characterise the poor-prognosis CMS4 consensus subtype (Guinney et al., *Nat Med* 2015)" (FR §9.2; SRF §11: "stromal / mesenchymal programmes; CMS4 subtype") | L | R1: CMS4 is mesenchymal, with stromal invasion, and has worse OS and RFS | **Supported** | Yes. Add that **GSE39582 contributed to the CMS derivation (R1)**, so this context is **not independent** of the present cohort |
| B2 | "In bulk tissue, EMT gene expression also reflects stromal and fibroblast content" (FR §9.2; SRF §11 "Cannot be concluded: tumour-cell EMT vs stroma") | L (indirect) | R2, R3 show that mesenchymal/poor-prognosis subtype transcripts in bulk CRC are mostly stromal. Neither tested the Hallmark EMT set specifically | **Partially supported** | Yes: "mesenchymal / EMT-like expression in bulk colorectal tumours largely reflects stromal (fibroblast) content (R2, R3); this was not tested for the Hallmark EMT score here" |
| B3 | "stromal-rich tumours may recur more" (SRF §11, Reasonable hypothesis, EMT) | L + H | R2 (stromal markers and relapse), R3 (CAF signature and poor prognosis) | Supported as a **literature-based hypothesis**. This dataset shows no association (HR 1.07, CI includes 1) | Add citations R2, R3; keep labelled as a hypothesis |
| B4 | "higher proliferation signatures have been reported in several studies to accompany *better* outcome, unlike many other cancers" (FR §9.2); "higher proliferation reported with better colorectal-cancer outcome in several studies" (SRF §11) | L, **mixed** | R4 supports lower proliferation in relapsed primaries (one small gene-signature study). R5 (meta-analysis, Ki-67) finds high Ki-67 associated with **worse** survival. "Unlike many other cancers" was not verified | **Overstated / partially supported** | **Yes.** "Evidence is mixed: one gene-expression study found lower proliferation in relapsed primary tumours (R4), whereas a meta-analysis of Ki-67 staining associated high proliferation with poorer survival (R5)." Delete "several studies" and "unlike many other cancers" |
| B5 | "lower proliferation may accompany higher recurrence" (SRF §11 hypothesis, E2F); "E2F (< 1) … biologically plausible (proliferation …)" (README) | L (mixed) + H | R4 for; R5 against (different measure) | Hypothesis only; plausibility is **contested** | Yes: "consistent with some, but not all, prior reports (R4 vs R5)" |
| B6 | "Pathway activation, usually through APC loss, is nearly universal in colorectal cancer" (FR §9.2); "activation near-universal in colorectal cancer" (SRF §11) | L | R6: Wnt pathway **genomically altered** in 93%; APC/CTNNB1 in about 80% | Supported **for genomic alteration**. "Activation" is a functional claim and should be narrowed | Yes: "Wnt-pathway genes are genomically altered in about 93 % of colorectal cancers (APC or CTNNB1 in about 80 %; R6)" |
| B7 | "…which limits between-tumour contrast" (FR §9.2) | H, **contradicted in part** | R1: Wnt/MYC signalling activation is a *distinguishing* feature of CMS2, so Wnt-related expression does vary between tumours | **Unsupported as stated** | **Remove.** Do not use "limited contrast" to explain the null Wnt coefficient |
| B8 | "none, a null is expected" (SFA §6, Wnt) | H, unsupported | as B7 | **Unsupported** | Remove (planning document) |
| B9 | "The gene set captures p53 transcriptional targets. TP53 mutation is not equivalent to low target expression" (FR §9.2; SRF §11) | M + L | R7 defines Hallmark sets as coherent expression programmes; R6 gives TP53 mutation in 60 % (non-hypermutated) | Methodologically sound | Add: "TP53 is mutated in about 60 % of non-hypermutated colorectal cancers (R6); the score is not a mutation readout, and TP53 status was not analysed here" |
| B10 | "PIK3CA mutations occur in a minority of colorectal cancers" (FR §9.2); "PIK3CA alteration in a minority of tumours" (SRF §11) | L | R6: PIK3CA mutated in 15 % of non-hypermutated tumours | Supported | Add the figure and citation |
| B11 | "Proliferation (E2F) and EMT programmes are plausible recurrence correlates" (FR §14) | L + H | B1–B5 | EMT/stroma plausibility supported; proliferation contested | Yes: "EMT/stromal programmes are plausible recurrence correlates (R1–R3); evidence for proliferation is mixed (R4, R5)" |
| B12 | "ideally with orthogonal measures such as stromal content, mutation data or protein assays" (FR §14) | H (recommendation) | — | Appropriate | No |
| B13 | Score = mean standardised mRNA; not activation, protein activity, signalling flux or mutation status (FR §4, §9.2; SRF §5; CAP Fig 4; README) | M | definition | Correct | No |
| B14 | "No pathway shows an independent association" (primary) and five-set E2F as a sensitivity result (FR §9, §9.1; SRF §7; CAP Fig 6; README) | D | EX JSON; ER §5 | Correct | No |
| B15 | Hallmark gene-set source (FR §4; README "Liberzon et al., *Cell Syst* 2015") | L | R7 | Supported | Add DOI |
| B16 | Dataset source (README "Marisa L. et al., *PLoS Med* 2013;10:e1001453") | L | R8 | Supported | Add DOI |

No document claims pathway activation, causality, protein activity or biomarker validity. All such statements are already negated.

---

## 3. Scientific-wording audit

| ID | Issue | Locations | Status | Proposed wording |
|---|---|---|---|---|
| W1 | "Independently reproduced" / "INDEPENDENTLY REPRODUCED" could be read as external replication | FR line 7 ("Independent reproduction"), §Summary row "Verification", §7.2 ("Classification: INDEPENDENTLY REPRODUCED"), §14 ("independently reproduced 95 % percentile intervals"), §14 status item 2; CAP Fig 8 ("and independently reproduced"); README status line 6 ("independent checkpoint-level audit") | **Needs clarification** (SRF §9 already says "by a separate implementation within the project") | At first use in each document: "independently reproduced — an independent re-implementation of the calculations within this project from the raw checkpoints, not an external replication by another group" |
| W2 | "five independent runs" | README (findings table, method sketch); FR Summary row and §6 "five independent runs", §14 "Five independent, verified …" | **Needs clarification** (SRF and CAP already say "differ only in random seed") | "five runs that differ only in random seed (Monte Carlo replicates of one analysis, not independent cohorts or experiments)" |
| W3 | Cause-specific hazard vs absolute probability / cumulative incidence | FR §2, §13; SRF §2, §12; CAP Fig 1, Fig 5, S2; README | Correct | — |
| W4 | HR association vs predictive discrimination | SRF §6; CAP Fig 5 | **Missing in FR** | Add to FR §9: "Hazard ratios describe adjusted association; added predictive discrimination is assessed by ΔC (§8), and the two need not agree." |
| W5 | Plausibility vs demonstrated pathway activity; mRNA vs protein, activation or mutation | all documents | Correct | — (B4–B7 narrow the plausibility statements) |
| W6 | Five-set E2F vs primary four-set | FR, SRF, CAP, README | Correct | — |
| W7 | S-7 vs independent replication | FR §11, SRF §7, CAP Fig 6 | Correct | — |
| W8 | Historical single-layer estimates vs corrected uncertainty | FR §9 heading, CAP Fig 6 and Fig 9 | **Missing in SRF §7, §8, §10 headings** | Add "(historical single-layer estimates; no nested intervals)" to SRF §7, §8 and §10 headings |
| W9 | Internal CIT partition vs external validation | FR §13, SRF §2 and §12, CAP Fig 1, README | Correct | — |

---

## 4. Unsupported or partially supported claims (summary)

1. **B4 (proliferation and better outcome, "several studies", "unlike many other cancers"):** overstated. Evidence is mixed (R4 vs R5).
2. **B7 (Wnt activation "limits between-tumour contrast"):** unsupported and partly contradicted by R1 (CMS2). Remove.
3. **B8 (SFA: Wnt "a null is expected"):** unsupported. Remove.
4. **B6 ("activation nearly universal"):** narrow to "genomically altered in about 93 %" (R6).
5. **B2 (EMT score reflects stroma):** supported for mesenchymal subtype signatures (R2, R3), not tested for the Hallmark EMT score. Narrow.
6. **B1 (CMS4 prognosis):** supported, but **not independent of GSE39582** (R1 included it). Disclose.

---

## 5. Approved documentation-only change set (applied 8 November 2026)

The owner approved groups A–E and optional changes O1–O3 on 8 November 2026, and they were applied on the same day. None changes a number, model, endpoint, interval, coefficient or ΔC conclusion. Codes: **U** corrects an unsupported claim; **N** narrows an interpretation; **P** methodological precision or disclosure.

**Corrections to the first draft of this section.** Checking the draft against the documents found five issues, all fixed below:
1. The draft labelled SRF §8 "historical single-layer estimates". Those corrected C values define θ̂ and are not labelled historical in FR §7.3, so the heading now says "full-data point estimates that define θ̂".
2. The draft cited R7 for "p53 transcriptional targets". R7 defines the Hallmark collection in general and was not verified for that statement, so the citation was dropped.
3. The draft FR EMT row did not narrow "EMT gene expression also reflects stromal and fibroblast content" (B2). F7 does.
4. The draft applied W1 in FR §7.2 only; line 7, the Summary verification row and §14 item 2 are now included.
5. The SFA Wnt row also says "near-universal activation" (B6), so it is annotated too.

### A. `docs/SCIENTIFIC_RESULTS_AND_FIGURES.md`
- **A1 (P, W8)** §7 heading: append "; historical single-layer estimates" inside the parentheses.
- **A2 (P)** §8 heading: "(Figure 7; full-data point estimates that define θ̂; the interval for ΔC is in §9)".
- **A3 (P, W8)** §10 heading: "(Figure 9; historical single-layer estimates, descriptive)".
- **A4 (P, W1)** §9: "…within the project (`reports/phase9/independent_audit/`), not an external replication by another group or in another cohort."
- **A5** §11 table, context and hypothesis cells only (number columns unchanged); header "Plausible context (literature; not tested here)":
  - EMT context (N, P; R1–R3): "mesenchymal/stromal expression characterises the CMS4 subtype, which had worse relapse-free survival (R1; GSE39582 contributed to the CMS derivation, so this is not independent evidence); mesenchymal-subtype transcripts in bulk colorectal tumours are mostly stromal (R2, R3), not tested for the Hallmark EMT score". Hypothesis: "stromal-rich tumours may recur more (R2, R3)".
  - E2F context (U; R4, R5): "proliferation and cell-cycle genes; evidence on proliferation and outcome is mixed: lower proliferation-signature expression in relapsed primary tumours (R4), but high Ki-67 associated with poorer survival in a meta-analysis (R5)". Hypothesis (N): "lower proliferation may accompany higher recurrence (consistent with R4, not with R5)".
  - Wnt context (U, N; R6, R1): "Wnt-pathway genes genomically altered in about 93 % of colorectal cancers (R6); genomic alteration is not activation, and marked WNT/MYC activation characterises CMS2 specifically (R1), so Wnt-related expression varies between tumours".
  - p53 context (P; R6): "p53 transcriptional targets, not TP53 mutation status; TP53 is mutated in about 60 % of non-hypermutated tumours (R6), and TP53 status was not analysed here".
  - PI3K context (P; R6): "PIK3CA mutated in about 15 % of non-hypermutated tumours (R6); gene expression is not pathway activity".
- **A6 (P)** new "§14. References" (R1–R8).
- **O1 (N)** §11 closing: "That makes them **hypotheses** consistent with part of the prior literature (EMT/stroma: R1–R3; proliferation: mixed, R4 vs R5)." §13 Biology bullet: the same qualifier replaces "biologically plausible hypotheses".

### B. `docs/PHASE9_CORRECTED_FINAL_REPORT.md`
- **F0 (P)** header revision note: wording and citations only; no number, model, endpoint or conclusion changed.
- **F1 (P, W1)** line 7: "Here "independent" means a separate re-implementation within this project, run on the same raw checkpoints; it is not an external replication by another group or in another cohort."
- **F2 (P, W2)** Summary row label: "five runs differing only in random seed".
- **F3 (P, W1)** Summary verification row: "an independent re-implementation within the project reproduced all values (not external replication)".
- **F4 (P, W2)** §6: "**five runs that differ only in random seed** (Monte Carlo replicates of one analysis of one cohort, not independent cohorts or experiments)".
- **F5 (P, W1)** §7.2: "This classification means the calculations were re-implemented separately within this project and reproduced from the same raw checkpoints. It is not an external replication: no other group, cohort or dataset was involved."
- **F6 (P, W4)** §9: "**Association, not discrimination.** Hazard ratios describe adjusted association. Added predictive discrimination is assessed by ΔC (§8), and the two need not agree."
- **F7** §9.2 intro cites "references R1–R8 at the end of this report"; Background cells:
  - EMT (N, P): CMS4 with worse relapse-free survival (R1); GSE39582 used to derive the CMS, so not independent evidence; mesenchymal-subtype transcripts mostly stromal (R2, R3); not tested for the Hallmark EMT score.
  - E2F (U): mixed evidence (R4 vs R5); neither measured E2F-target mRNA. "several studies" and "unlike many other cancers" removed.
  - Wnt (U, N): genomically altered in about 93 %, APC/CTNNB1 in about 80 % (R6); alteration is not activation; Wnt-related expression varies (CMS2; R1). "which limits between-tumour contrast" removed.
  - p53 (P): TP53 mutated in about 60 % of non-hypermutated tumours (R6); TP53 status not analysed.
  - PI3K (P): PIK3CA mutated in about 15 % of non-hypermutated tumours (R6).
- **O1 (N)** §9.2 closing: "**hypotheses consistent with part of the prior literature** (EMT/stroma: R1–R3; proliferation: mixed, R4 vs R5)".
- **F8 (P, W1, W2)** §14 Observed: "Five verified 95 % percentile intervals, from runs that differ only in random seed and reproduced by an independent re-implementation within the project, …".
- **F9 (N, B11)** §14: "EMT/stromal programmes are plausible recurrence correlates (R1–R3); evidence on proliferation (E2F) is mixed (R4, R5) (§9.2)."
- **F10 (P, W1)** §14 item 2: "resolved by the verified nested bootstrap, reproduced by an independent re-implementation within the project (§8)."
- **F11 (P)** new "References" section (R1–R8).

### C. `README.md`
- **M1 (P, W1)** status: "an independent re-implementation within this project reproduced every reported number from the raw checkpoints (not an external replication)".
- **M2, M3 (P, W2)** findings table and method sketch: "five runs differing only in random seed".
- **M4 (U, N)** interpretation bullet: directions not demonstrated; EMT consistent with Guinney 2015 (whose derivation included GSE39582), Calon 2015, Isella 2015; Hallmark EMT vs stroma not tested; proliferation evidence mixed (Anjomshoaa 2009; Luo 2019); pointer to SRF §14.
- **M5, M6 (P)** DOIs for Marisa 2013 (R8) and Liberzon 2015 (R7).
- **O2 (P, W4)** "(hazard ratios measure association; added discrimination is the ΔC above)".
- **O3** status line corrected: branch pushed and private GitHub Release `phase9-corrected-campaign-2026-11-08` published; licence and public-release review pending.

### D. `docs/SCIENTIFIC_RESULTS_FIGURE_AUDIT.md` (planning record)
- Wnt row annotated, not deleted: B6 note on "near-universal activation"; "a null is expected" marked as withdrawn (B8).

### Figure 8 caption (`scripts/49`, regenerated `reports/figures/final/`)
- "and independently reproduced." → "and reproduced by an independent re-implementation within this project from the same checkpoints (not external replication)."
- Only string literals changed. The original output was backed up and hash-verified before regeneration; all 22 images and `figure_values.json` regenerated byte-identical; only the Fig 8 caption line and the script hash in `figure_sources.json` changed.

---

## 6. Unresolved evidence gaps

1. **Hallmark EMT vs stroma.** No reference tests whether the **Hallmark EMT** score itself is stroma-driven in GSE39582. R2 and R3 concern subtype signatures. Settling it would need deconvolution or stromal-content estimates, which would be a new analysis and is not proposed.
2. **E2F direction.** The literature conflicts (R4 vs R5), and they measure different things (a gene signature vs Ki-67 protein). No E2F-Hallmark-specific prognostic study was verified.
3. **CMS circularity.** CMS prognostic evidence includes GSE39582 (R1), so external corroboration would need studies that exclude this cohort. None was verified here.
4. **TP53 and PIK3CA.** Mutation status is available in the GEO metadata (`tp53.mutation`, etc.) but was not analysed. Any link between scores and mutations would be a new analysis.
5. **Unverified candidate.** Samuels et al. 2004 is metadata-verified only, and is not used.
