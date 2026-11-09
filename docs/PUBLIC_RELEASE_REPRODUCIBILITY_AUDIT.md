# Public-release and reproducibility audit

**Date:** 8 November 2026. **Type:** audit and recommendation only. Nothing was deleted, redacted, moved, committed, pushed or published. Repository visibility, the existing GitHub Release and its asset are unchanged. This document is local and uncommitted.

**Scope:**
- the whole tracked repository (487 files) and its complete Git history, including the local-only `main` branch;
- the GitHub Release `phase9-corrected-campaign-2026-11-08` and its asset;
- the ignored local data directories, for completeness.

**Final decision: READY AFTER SPECIFIED REMEDIATIONS** (section 10).

---

## 1. Executive assessment

The scientific content is defensible and conservatively stated. Integrity is verified end to end. No credentials, secrets or non-public patient data were found. Public release is nevertheless **not yet responsible** without remediation.

**Blocking items:**
1. **The repository has no licence.** Public visibility without one grants no reuse rights.
2. **Per-sample tables need a data-governance decision before visibility changes.** The repository and its pushed history contain derived per-sample tables of clinical, outcome, mutation and expression data from GSE39582. The source is already public and de-identified, but some tables are new derived inferences about individual samples, notably the expression-based sex-discordance review. Making the repository public publishes the **entire history**. Any later removal would then require a history rewrite or a new repository, so this decision must come first.
3. **Third-party attribution is missing.** No notice covers:
   - the MSigDB Hallmark GMT files redistributed under CC BY 4.0;
   - the GEO, ArrayExpress and PLoS-supplement derived content.
4. **The reproduction instructions overstate portability.** By design, the frozen launcher, verifier and report tool refuse to run outside the original interpreter and virtual environment. Only the independent audit (`scripts/48`) is portable.
5. **Two documents are stale about the archive.** They state the archive has not been built or uploaded; it has been, privately.

Several medium and low presentation items could also be challenged by reviewers (section 8).

---

## 2. Repository and release integrity (verified 8 Nov 2026)

| Item | Result |
|---|---|
| Branch / HEAD | `phase9-corrected-pipeline` @ `1146afade28a87fd62dfc1f42f778288b1298453` |
| Remote branch / remote HEAD | `1146afa…` (identical); working tree clean |
| Tag | `phase9-corrected-campaign-2026-11-08` → `1146afa…` (commit) |
| Release asset | `phase9_campaign_artifacts.tar`, 579,716,096 bytes; GitHub-reported digest `sha256:0812b3ba…c2aa2f`, equal to the local build |
| Archive manifest | `reports/phase9/ARCHIVE_MANIFEST_SHA256.txt` SHA-256 `c45bb806…6b6f`; 52/52 entries verify against the files on disk |
| Checkpoints | five `runs/*.jsonl` on disk, hashes equal to their end records, ignored, never committed |
| Repository visibility | PRIVATE |
| Git history | 14 commits on the pushed branch, all authored and committed by "Prasad Akula <akula.pra@northeastern.edu>". **No tool-attribution trailers in pushed history.** Two commits with a `Co-Authored-By` tool trailer (`6fda652`, `1b8c1c2`) exist **only on the local, unpushed `main` branch** |
| Secrets scan (all history) | 0 hits for GitHub, AWS, Slack, OpenAI-style or private-key patterns; no `.env`, `.pem` or `.key` files tracked |

---

## 3. Data-disclosure and privacy audit

Classification codes: (1) safe to publish; (2) publishable with context or documentation; (3) requires removal or redaction; (4) requires legal, licensing or data-governance review; (5) must remain private.

| # | Finding | Location | Class | Notes |
|---|---|---|---|---|
| D1 | Per-sample GEO clinical and outcome metadata (585 samples × 37 columns: sex, age, TNM, location, chemotherapy, RFS/OS event and time, MMR/CIMP/CIN, TP53/KRAS/BRAF mutation detail, CIT subtype) | `metadata/GSE39582_sample_metadata.csv` (pushed since commit `1cd2bba`) | **4** | A direct mirror of the public GEO series characteristics. It is de-identified (GSM accessions, no names or dates of birth), but combines age, sex, stage, mutations and survival per sample. GEO terms could not be re-verified this session (the page returned a captcha). The recalled position is that NCBI places no restrictions on use but submitters may claim rights. Needs owner and institution sign-off, plus a provenance notice |
| D2 | Derived per-sample cohort and endpoint tables | `reports/phase1/cohort_manifest.csv` (585 × 26), `excluded_records.csv` (8 × 53), `zero_time_records.csv` (8 × 16), `five_year_reclassified.csv` (6 × 8) | 2 | Subsets and recodings of D1 |
| D3 | Table S1 values from the source publication | `reports/phase1/zero_time_vs_tableS1.csv` (8 × 9) | 2 | Values from Marisa et al., *PLoS Med* 2013 supplement (PLOS articles are CC BY). The raw `.xls` is ignored and not redistributed. Needs attribution |
| D4 | Expression-derived per-sample QC, PCA and distribution metrics | `reports/phase7c/*` (453-row files), `reports/phase8a/sample_*`, `pca_scores_gene_level.csv`, `distributional_heterogeneity.csv` | 2 | Technical summaries of public arrays. Not identifying beyond the public GSM |
| D5 | Full 453 × 453 sample correlation matrix | `reports/phase8a/sample_correlation_pearson.csv` | 2 | Derived; no clinical content. Large but under limits |
| D6 | Per-sample pathway scores (exploratory, full cohort) | `reports/phase8a/pathway_scores_*`, `reports/phase8c/scores_*` | 2 | Derived from public expression; labelled EXPLORATORY |
| D7 | **Expression-inferred sex vs recorded sex, per sample**, including 8 discordant samples with QC metrics | `reports/phase8a/sex_positive_control.csv` (453 × 8), `reports/phase8b/sex_discordant_eight_review.csv` (8 × 15) | **4** | A new derived inference about individuals, not present in GEO. Discordance can reflect sample mix-up, but it could also be read as a statement about a person's sex or gender. Recommend an aggregate-only replacement (counts and summary) in any public view, subject to owner and data-governance decision |
| D8 | Per-sample data for **other** cohorts (GSE14333/GSE17536 overlap: age, sex, stage) | `reports/phase4/overlap_GSE14333_GSE17536_*.csv` (133 and 129 rows × 9) | 2 / 4 | Public GEO data from other submitters. Same considerations as D1 |
| D9 | Individual GSM accessions cited in prose | `docs/PHASE1…`, `PHASE7C…`, `PHASE8A…`, `PHASE8B…`, `config/analysis/QC_DECISIONS_RECORD.json`, etc. (≤ 38 IDs each) | 1 / 2 | Public accessions used to document QC decisions |
| D10 | Checkpoint archive content | release asset | 1 (data) | No sample IDs: resamples are stored as position-indexed multiplicity vectors and hashes. No local paths. Releases follow repository visibility, so **making the repository public makes this asset public** |
| D11 | Local filesystem paths and username (`/Users/kalki/…`, `/Users/kalki/miniforge3/…`) | 27 tracked files: the frozen freeze record and A-19/M2 freeze records, `QC_SPECIFICATION.json`, Phase 2–9 summary JSONs, `final_project/audit/pytest_project_regression.log`, `final_project/reproducibility/environment.json`, five historical script docstrings (scripts 33–37), three closure docs | 2 | Not secret; reveals a username and directory layout. **Frozen and hash-pinned files cannot be edited** without breaking verification. Disclose rather than redact |
| D12 | Institutional email in commit metadata | all commits | 2 | Standard Git practice; owner's choice whether to keep it in a public history |
| D13 | Tool-attribution trailers | local-only `main` (`6fda652`, `1b8c1c2`) | 5 (keep unpushed) | Never push `main`. The owner may delete the local branch or keep it local |
| D14 | Raw and processed data, external downloads (SOFT files, Parquet bundle, SDRF, Table S1 `.xls`) | `data/` (ignored) | 5 / not in repository | Correctly excluded |
| D15 | Secrets and credentials | — | — | None found |

---

## 4. Per-sample data inventory (tracked files with ≥ 5 GSM identifiers)

"Public?" means the underlying values are already public at the source. "Aggregate-only OK?" means an aggregate replacement would still support the reported analysis.

| File | Format | Rows × cols | Variable types | Identification risk | Public? | Redistribution | Aggregate-only OK? |
|---|---|---|---|---|---|---|---|
| `metadata/GSE39582_sample_metadata.csv` | CSV | 585 × 37 | clinical, outcome, mutation, subtype | indirect (quasi-identifier combination on a public accession) | yes (GEO) | review (D1) | no: it is the reproducible metadata snapshot; re-downloadable from GEO by `scripts/01` |
| `reports/phase1/cohort_manifest.csv` | CSV | 585 × 26 | endpoint, stage, cohort flags | indirect | derived from public | with context | partly; regenerable |
| `reports/phase1/excluded_records.csv` | CSV | 8 × 53 | full clinical rows | indirect | yes | with context | yes |
| `reports/phase1/zero_time_records.csv` | CSV | 8 × 16 | clinical, endpoint | indirect | yes | with context | yes |
| `reports/phase1/zero_time_vs_tableS1.csv` | CSV | 8 × 9 | endpoint vs publication | indirect | yes (CC BY supplement) | with attribution | yes |
| `reports/phase1/five_year_reclassified.csv` | CSV | 6 × 8 | endpoint | indirect | yes | with context | yes |
| `reports/phase4/overlap_GSE14333_GSE17536_all.csv`, `_unambiguous.csv` | CSV | 133 × 9, 129 × 9 | age, sex, stage (other cohorts) | indirect | yes (GEO) | with context | yes (counts suffice) |
| `reports/phase7c/sample_qc_metrics.csv`, `degradation_ratios.csv`, `pca_scores.csv`, `investigation_summary.csv` | CSV | 453 × 10, 453 × 3, 453 × 8, 29 × 7 | technical QC, batch, subtype | low | derived | yes | yes for public view; per-sample needed to audit QC decisions |
| `reports/phase8a/sample_distribution_summary.csv`, `distributional_heterogeneity.csv`, `sample_mean_correlation.csv`, `pca_scores_gene_level.csv` | CSV | 453 × 20, 453 × 6, 453 × 2, 453 × 24 | expression summaries, batch | low | derived | yes | yes |
| `reports/phase8a/sample_correlation_pearson.csv` | CSV | 453 × 454 | expression correlation | low | derived | yes | yes |
| `reports/phase8a/pathway_scores_*_EXPLORATORY_full_cohort.csv`; `reports/phase8c/scores_*_reproduced_full_cohort.csv` | CSV | 453 × 5–8 | pathway scores, batch | low | derived | yes | partly (scores are reproducible from the bundle) |
| `reports/phase8a/sex_positive_control.csv` | CSV | 453 × 8 | inferred vs recorded sex | **sensitive derived inference** | no (new inference) | **review (D7)** | **yes, recommended** |
| `reports/phase8b/sex_discordant_eight_review.csv` | CSV | 8 × 15 | inferred sex discordance, QC | **sensitive derived inference** | no | **review (D7)** | **yes, recommended** |
| `reports/phase8b/mlh1_axis_by_mmr_group.csv` | CSV | 453 × 4 | expression axis, MMR | low | derived | yes | yes |
| `reports/phase8c/score_influence_loo.csv` | CSV | 5 × 6 | per-sample influence | low | derived | yes | yes |
| JSON summaries (`phase1_findings.json`, `qc_summary.json`, `investigation_records.json`, `phase8a/8b/8c_summary.json`, `cohort_summary.json`, `QC_DECISIONS_RECORD.json`, test fixture) | JSON | — | IDs embedded in findings | low | derived | yes | n/a |
| Phase docs listing GSMs | MD | — | IDs in prose | low | yes | yes | n/a |
| **No** model predictions or per-patient risk scores are stored in tracked files. Bootstrap and validation records (checkpoints) contain no sample IDs | | | | | | | |

**Removal is not reversible after publication.** Every file above is already in the pushed branch history. If the owner decides any of them must not be public, the options are:
- (a) rewrite history before visibility changes, which needs explicit approval and invalidates the commit hashes cited in the documents and the release tag; or
- (b) publish a **new, curated public repository or snapshot** while keeping this repository private as the full audit record. **Option (b) is recommended**, because it preserves the hash-anchored evidence trail.

---

## 5. Reproducibility assessment

| Category | What | Basis |
|---|---|---|
| **Fully reproducible from the repository alone** | Code review; the synthetic-data unit and integration tests (most of the 1,085; some need local data); documentation and decision trail | `README.md` §Environment; `tests/` |
| **Reproducible with documented external inputs** | The Phase 1–9 pipeline (scripts 01–39) and θ̂: GEO SOFT files (`scripts/01`, `09`; SHA-256 in `data/raw/CHECKSUMS.sha256`); ArrayExpress SDRF (auto-fetched by `src/external_annotation.py`; needed for ComBat batch labels); Table S1 `.xls` (auto-fetched by `scripts/15`). The **independent audit** (`scripts/48`, standard library only) **given the checkpoint archive**. The archive is currently in a *private* release, so external researchers cannot obtain it yet | `scripts/01`, `09`, `13`, `16–23`, `38`; `scripts/48`; `docs/LARGE_ARTIFACT_ARCHIVE.md` |
| **Reproducible only with the original execution environment** | The frozen launcher (`scripts/42`), the post-run verifier (`scripts/43`), the gated report tool (`scripts/46`), and bit-identical re-execution of the corrected campaign. `verify_code_identity` compares `ENV_KEYS` = interpreter path and SHA-256 (`/Users/kalki/miniforge3/bin/python3.13`), virtual-environment prefix (`/Users/kalki/coad_transportability/.venv`), installed distributions, NumPy BLAS/LAPACK build (Accelerate) and platform, so on another machine it is expected to refuse. This is reasoned from code; it was not tested on another machine | `src/phase9_corrected/orchestration.py:1230–1231, 1457–1462`; freeze record `environment` |
| **Not independently reproducible** | The depositors' RMA + ComBat preprocessing (the project uses the deposited matrix; CEL-level reprocessing not done). Independent investigator replication: the "independent audit" is an independent *implementation* by the same project, not an independent team | `docs/PHASE2_GENE_REPRESENTATION.md`; IA §1 |

**Further reproducibility observations:**
- **Determinism.** Seeds are fixed and documented, both the per-replicate `SeedSequence` design and the Phase 9 seeds 20261102–20261106. Bit-identical results are demonstrated only on the original machine (V5). A different BLAS or CPU may change floating-point results slightly; this is not quantified.
- **OS assumptions.** The scripts assume macOS arm64. The launcher needs system `jq` and `shasum`, and the durability code uses `F_FULLFSYNC` (macOS-specific). Linux portability is untested.
- **Execution order.** The README gives "scripts 13, 16–23 in order" without the exact command list or the external fetch steps. No single driver script or Makefile exists.
- **Absolute paths in executable code:** none. Only docstrings in scripts 33–37 contain absolute example paths.
- **Tests:** 1,085 passed on the original machine. Passing tests establish software conformance, not reproducibility elsewhere.

---

## 6. Scientific presentation review

**Accurate and appropriately restrained:**
- the research question;
- the cohort (452 patients / 136 events) and the endpoint;
- the M1 and four-set M2 definitions;
- the separation of primary, five-set and S-7 analyses;
- the corrected ΔC and its five intervals;
- the paired-interval distinction (FR §8.1);
- the cause-specific vs CIF vs composite distinction (FR §2);
- batch instability, informative censoring and overfitting;
- an explicit "not demonstrated / not equivalence" conclusion with no utility or causal claims.

**Language a reviewer could reasonably challenge:**

| # | Text | Location | Concern | Suggested change |
|---|---|---|---|---|
| P1 | "how to evaluate an omics signature honestly rather than report its apparent performance" | `docs/FINAL_PROJECT_HANDOVER.md` §5 (portfolio description) | Implies others are dishonest; self-congratulatory | "…how optimism correction and nested resampling change the apparent value of an omics signature" |
| P2 | "unusually strong computational controls" | `docs/PHASE9_SCIENTIFIC_INTERPRETATION_REVIEW.md` §7 | Unsupported comparative claim | "extensive computational controls" |
| P3 | "INDEPENDENTLY REPRODUCED" (headline) | README, FR, IA | Could be read as independent investigators | Add "independent re-implementation by the project (not an external replication)" at first use |
| P4 | "five independent runs" | README, FR | They are independent Monte Carlo replicates of one analysis, not independent studies | Add "(differing only in random seed)" at first use |
| P5 | Biology background without citations (proliferation paradox, APC ubiquity, PIK3CA frequency) | FR §9.2 | A supervisor will ask for sources | Add 2–3 primary citations, or mark "general background, uncited" explicitly |
| P6 | Hazard ratios with "95 % Wald CI" from a ridge-penalised fit | FR §9, README | Penalised Wald intervals are not standard inference | Already caveated in FR §9; add "(approximate)" in the README table |
| P7 | Protocol files named `*_DRAFT.md` although approved; two "Decision N" numbering systems | `docs/PHASE6B_…_DRAFT.md`, `PHASE9_MODELING_PROTOCOL_DRAFT.md`; register; FR §15 | Confusing to a reader | A short "document status" note in the README (files are hash-pinned and cannot be renamed) |
| P8 | Archive described as "not built or uploaded" | `docs/LARGE_ARTIFACT_ARCHIVE.md` §Status; `docs/release/RELEASE_READINESS.md` §1.0 | Stale since the private release | Update to record the release, its tag and digest |
| P9 | "Without re-running, from the archived checkpoints: run the verifier, the report tool and the independent audit" | `README.md` §Reproduction | The verifier and report tool refuse outside the original environment (section 5) | State that only `scripts/48` is portable; the verifier and report tool reproduce on the original environment only |
| P10 | Historical provisional package (`reports/phase9/final_project/`), including a presentation with superseded "PROVISIONAL" numbers | `reports/phase9/final_project/` | A casual reader may cite superseded numbers | Keep, but add a "superseded" notice in a public README or exclude it from a curated public snapshot |

The negative or inconclusive result is presented as a legitimate finding (FR §14; README "Main findings").

---

## 7. Licensing and third-party restrictions

**Current state:** no licence file and no third-party notice. Copyright statements appear only inside third-party package metadata, not in repository files.

| Component | Terms | Status in repository | Obligation |
|---|---|---|---|
| Original code (src, scripts, tests) | none chosen | — | Owner chooses a licence |
| Documentation and reports | none chosen | — | Owner may choose a separate content licence (for example CC BY 4.0) |
| GSE39582 data (GEO) and derived tables | NCBI GEO disclaimer (recalled: no NCBI restriction; submitters may claim rights). Could not be re-verified this session. Associated publication: Marisa et al., *PLoS Med* 2013 (CC BY) | metadata mirror plus derived per-sample tables | Cite the source. Data-governance decision on per-sample tables (D1, D7). **Re-verify the GEO terms at the source** |
| GSE14333 / GSE17536 (GEO) | as above | per-sample overlap tables | Cite; same decision |
| ArrayExpress SDRF E-GEOD-39582 | EMBL-EBI terms of use | not redistributed (ignored); derived batch labels appear in per-sample reports | Cite |
| Marisa 2013 Table S1 | CC BY (PLOS) | `.xls` not redistributed; 8-row derived comparison table | Attribution |
| GPL570 platform annotation | GEO platform record (Affymetrix/Thermo Fisher annotation) | derived probe-mapping tables in `reports/phase2/` | Cite; **licensing review** of redistributing derived annotation tables |
| **MSigDB Hallmark gene sets v2026.2.Hs** | CC BY 4.0 (MSigDB terms page, checked 8 Nov 2026: applies to v2022.1 and later; the exceptions listed are KEGG, BioCarta and STKE, and Hallmark is not among them) | **full GMT files redistributed**: `reports/phase6c/h.all.v2026.2.Hs.{entrez,symbols}.gmt`; panel JSONs | **Attribution required** (MSigDB; Liberzon et al., *Cell Syst* 2015) and indication of changes |
| Python dependencies | numpy (BSD-3 and others), pandas (BSD-3), scipy (BSD-3), lifelines (MIT), pyarrow (Apache-2.0), pytest (MIT), GEOparse (BSD), xlrd (BSD), matplotlib (PSF-style), scikit-learn (BSD-3), tqdm (MPL-2.0 AND MIT) | **not redistributed**; only pinned in `requirements*.txt` | None for listing; no copyleft code is vendored |
| Checkpoint archive | project-generated; no third-party data rows (no IDs, only resample multiplicities and statistics) | private release asset | Covered by whatever licence the owner chooses for project outputs |

**Code-licence options** (for original code only; third-party data and annotations are not relicensed by any of these):

| Licence | Reuse | Patent grant | Copyleft | Fit |
|---|---|---|---|---|
| MIT | very permissive | no explicit grant | no | simplest; common for portfolios |
| BSD-3-Clause | very permissive, plus a no-endorsement clause | no explicit grant | no | matches the NumPy, pandas and SciPy ecosystem |
| Apache-2.0 | permissive | **explicit patent grant** | no | stronger legal clarity; longer |
| GPL-3.0 | derivatives must stay GPL | yes | **strong** | limits reuse; unusual for a portfolio analysis |

Content (docs and reports) could use **CC BY 4.0**, consistent with the MSigDB and PLOS terms. **No licence has been selected.**

---

## 8. Prioritised findings

| ID | Severity | Location | Evidence | Why it matters | Remediation | Changes science? | Owner approval |
|---|---|---|---|---|---|---|---|
| F1 | **Critical (gating)** | `metadata/GSE39582_sample_metadata.csv`; `reports/phase8a/sex_positive_control.csv`; `reports/phase8b/sex_discordant_eight_review.csv`; other per-sample tables (section 4) | already in pushed history | Publication is irreversible and includes history. Derived sex-discordance inferences are new per-individual statements | Owner and institution data-governance decision **before any visibility change**. Recommended: keep this repository private as the audit record, and publish a curated public snapshot with aggregate-only replacements for D7 (and optionally D1/D2, with a GEO download instruction instead) | packaging only | **Yes** |
| F2 | **High** | repository root | no `LICENSE` | Public code without a licence cannot be reused legally | Owner selects a code licence (and optionally CC BY 4.0 for content) | packaging | **Yes** |
| F3 | **High** | `reports/phase6c/*.gmt`; panel JSONs; per-sample GEO-derived tables; Table S1 comparison | no attribution notice | CC BY 4.0 and source-citation obligations | Add a `NOTICE` / data-sources section: MSigDB (Liberzon 2015, CC BY 4.0, v2026.2.Hs), GEO GSE39582 / GSE14333 / GSE17536, ArrayExpress E-GEOD-39582, Marisa 2013 Table S1 | packaging | Yes (wording) |
| F4 | **High** | `README.md` §Reproduction ("Without re-running…") | `ENV_KEYS` include interpreter path, hash and prefix (`orchestration.py:1230`) | Overstates portability; a reviewer following it would hit a refusal | Reword: `scripts/48` is portable; the verifier and report tool run only in the original environment; describe what an external user can verify | packaging | Yes |
| F5 | **High** | release visibility coupling | releases follow repository visibility | A visibility change would publish the 580 MB archive and history at once | Decide visibility together with F1; for a curated public snapshot, decide whether to attach the archive there | packaging | **Yes** |
| F6 | Medium | `docs/LARGE_ARTIFACT_ARCHIVE.md` §Status; `docs/release/RELEASE_READINESS.md` §1.0; `docs/FINAL_PROJECT_HANDOVER.md` | "not yet built or uploaded" | Stale after the private release | Record the release tag, URL, asset size and SHA-256 | packaging | Yes (documentation commit) |
| F7 | Medium | `README.md` §Reproduction step 1 | "scripts 13, 16–23 in order"; external fetches undocumented; SDRF and Table S1 not in `CHECKSUMS.sha256` | Execution order and external inputs not fully specified | Add an exact ordered command list and an external-inputs table (URL, purpose, where its hash is recorded) | packaging | Yes |
| F8 | Medium | GEO terms | disclaimer page not retrievable this session | Redistribution basis asserted from memory | Re-verify the GEO and EBI terms at the source and cite them in the notice | packaging | Yes |
| F9 | Medium | GPL570-derived mapping tables (`reports/phase2/*`) | derived from vendor annotation | Redistribution terms unconfirmed | Licensing review, or exclude from a public snapshot (regenerable by `scripts/09`/`17`) | packaging | Yes |
| F10 | Medium | P1–P4 (section 6) | quoted text | Reviewer credibility | Wording changes | wording only (no numbers) | Yes |
| F11 | Low | P5–P7, P10 | section 6 | Clarity | Citations, status notes, superseded notice | wording only | Yes |
| F12 | Low | 27 files with `/Users/kalki/…` (D11) | `git grep` | Username exposure; frozen files cannot change | Accept and disclose; optionally exclude non-frozen historical logs from a curated snapshot | packaging | Yes |
| F13 | Low | commit metadata email (D12) | `git log` | Personal or institutional email public | Owner's choice; a history rewrite is not recommended | packaging | Yes |
| F14 | Low (process) | local `main` with tool-attribution trailers (D13) | `6fda652`, `1b8c1c2` | Must never be pushed | Keep local or delete the local branch; never push `main` | none | Yes (to delete) |
| F15 | Low | OS portability (`F_FULLFSYNC`, macOS `shasum`/`jq`, Accelerate BLAS) | code and freeze record | External users on Linux may fail | Document "tested on macOS arm64 only" | packaging | Yes |

No finding requires changing any scientific number, model, endpoint, cohort, protocol, configuration or checkpoint.

---

## 9. Recommended remediation sequence

1. **F1 + F5. Data-governance and visibility decision.** Choose between:
   - (a) a full public repository with per-sample tables (after institution sign-off); or
   - (b) **recommended:** keep this repository private as the authoritative audit record, and create a curated public snapshot. The snapshot would contain code, documentation, aggregate results and archive instructions, with D7 (and optionally D1, D2, D8 and the GPL570-derived tables) replaced by aggregates or download instructions.
2. **F2.** Select a code licence, plus a content licence if wanted.
3. **F3 + F8 + F9.** Re-verify third-party terms at the source and add the NOTICE / data-sources file.
4. **F4 + F6 + F7 + F15.** Documentation commit: portability statement, release record, exact reproduction commands, external-input table, tested-platform note.
5. **F10 + F11.** Wording corrections P1–P7 and P10.
6. **Final pre-publication check.** Re-run the secret and path scan and the test suite on the exact tree to be published, and confirm that the frozen hashes are unchanged.
7. **Only then:** a separate owner authorisation to change visibility or publish the snapshot. Do not push `main` (F14).

---

## 10. Go / no-go

**READY AFTER SPECIFIED REMEDIATIONS.**

The science, integrity and audit trail are defensible as they stand. Public release should wait for:
- the data-governance decision on per-sample tables (F1), decided together with visibility and the release asset (F5);
- a code licence (F2);
- third-party attribution and re-verified terms (F3, F8, F9);
- corrected reproduction and portability statements (F4, F6, F7).

A clean test suite and a verified private release do not establish privacy compliance, licensing compatibility or external reproducibility. This audit makes no change and authorises nothing.
