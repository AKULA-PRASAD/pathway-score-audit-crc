# Project completion inventory

**Prepared:** 7 November 2026; **final update 8 November 2026**. **Scope:** the whole repository (Phases 0–9), as it stands in the working tree.

**Closure documents (8 November 2026):** `docs/FINAL_PROJECT_HANDOVER.md` (consolidated handover), `docs/OWNER_DECISION_REGISTER.md` (open owner decisions D1–D5), `reports/phase9/independent_audit/` (independent audit: INDEPENDENTLY REPRODUCED).

**Owner approvals (8 November 2026):** decisions D1–D5 approved as recommended and the local commit plan authorised (`docs/OWNER_DECISION_REGISTER.md`, "Owner approvals recorded"). Push, release, archive upload and the licence remain open; items in section 5 that refer to D1–D5 are resolved accordingly.

**Scientific status:** the corrected campaign is complete. Five runs × 750 outer replicates passed the frozen verification at 750 (PASS, 60 of 60). Reporting was permitted, and the final interval is computed and reported in `docs/PHASE9_CORRECTED_FINAL_REPORT.md`. What remains are owner decisions on version control and archiving (section 5).
**Status of this file:** a review document for the project owner. Section 1 gives the repository state after the approved local commits on 8 November 2026. Sections 2–5 are retained as the pre-commit inventory on which the commit plan was based. Nothing has been deleted or moved, and nothing has been pushed.

## 1. Repository state

| Item | State |
|---|---|
| Branch | `phase9-corrected-pipeline`, tracking `origin/phase9-corrected-pipeline` (private GitHub repository) |
| Local commits not pushed | **10** on top of `origin` (`6bb3f26`): `c432eee` (Path C sizing-decision writer, 7 Nov 2026; kept unchanged, owner decision D1), then the eight approved commit groups (G1 `7e0df40`, G2 `1c3bcaf`, G3 `5169b9f`, G4 `6300e39`, G5 `f35346d`, G6 `4ae202c`, G7 `3924526`, G8 `97ed2a0`), then a documentation-only follow-up commit that records this state. Nothing has been pushed; the remote is unchanged |
| Local `main` | still at `1b8c1c2` (original Phase 0 commit, with the earlier co-author trailer); not pushed, not changed |
| Working tree | clean after the commits: no modified or untracked files. The five checkpoints `reports/phase9/corrected_rerun/runs/run0–4.jsonl` remain on disk, unchanged, and excluded from Git by `.gitignore` (owner decision D2) |
| Tracked files | 487 after commit G8 (88 before the commit plan) |

## 2. What was missing from version control (pre-commit inventory)

**Status (8 November 2026):** everything in this section except the five checkpoints has since been committed in groups G1–G8. The section is retained as the record the commit plan was based on.

Every file the corrected pipeline imports is tracked: no tracked module imports an untracked one, and all 27 frozen files are tracked. What is missing is the **evidence and history of Phases 1–9**, which exists only in the working tree.

### 2.1 Recommended for commit (small, text, needed for reproducibility and audit)

| Group | Files | Why it matters |
|---|---|---|
| Approved Phase 9 freeze | `config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json` | Trust anchor for every corrected run; its SHA-256 (`618f20eb…a790`) is recorded in the configuration. Without it a clone cannot run the launcher. |
| Phase 9 configuration change | the 4-line diff in `PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json` | Records Decision 13 approval and Decision 14 authorisation; the committed version still shows them pending |
| Other approved configuration | `config/analysis/PHASE6B_AMENDMENT_RECORD.json`, `QC_DECISIONS_RECORD.json`, `QC_SPECIFICATION.json`; `config/pathways/m2_candidate_panel.json`, `m2_panel_5set.json` | Decision records for the modelling amendment, QC and the pathway panels (the five-set panel is the frozen sensitivity panel) |
| Phase 1–9 documents | 25 files in `docs/` (Phase 1 to Phase 9 reports, plans and audits) | The protocol trail. Several are cited by tracked files (for example the README and the freeze procedure) |
| Phase 1–9 scripts | `scripts/15`–`39` (25 files) | Produce every Phase 1–9 artifact, including the Phase 9 primary execution (`38`) that is the source of θ̂ |
| Supporting modules | `src/aggregation.py`, `bundle.py`, `cohort.py`, `endpoint.py`, `external_annotation.py`, `gene_mapping.py`, `m2_amendment.py`, `precision_sim.py`, `provenance.py` | Imported by scripts 15–39 and by the Phase 3–9 tests |
| Phase 3–9 tests | 17 files plus `tests/fixtures/premove_scorer_reference.json` | Verification of the earlier phases |
| Final results and report | `docs/PHASE9_CORRECTED_FINAL_REPORT.md`; `docs/FINAL_PROJECT_HANDOVER.md`; `docs/OWNER_DECISION_REGISTER.md`; `reports/phase9/independent_audit/` and `scripts/48_phase9_independent_audit.py`; `reports/phase9/corrected_rerun_report/` (interval results JSON and tables, 70 KB); `reports/phase9/ARCHIVE_MANIFEST_SHA256.txt` | The verified final result and the hash manifest for the archive |
| New, this session | `scripts/46_phase9_corrected_report.py`, `tests/test_phase9_corrected_report_tool.py`, `scripts/47_phase9_theta_alignment_check.py`, `docs/PHASE9_THETA_ALIGNMENT_AUDIT.md`, `docs/LARGE_ARTIFACT_ARCHIVE.md`, this file | Final interval computation (gated on a PASS verification) and its tests; the θ̂ alignment audit and its read-only check; the archive proposal |
| Phase 1–9 summary outputs | `reports/phase1` … `reports/phase8d`, `reports/phase9/*.json|csv` (decision register, protocol change log, readiness checklist, summaries) | Small tables, JSON and figures (about 25 MB in total, the largest a 6 MB CSV) |

### 2.2 Needs an owner decision before commit (large or derived)

| Path | Size | Issue |
|---|---|---|
| `reports/phase9/corrected_rerun/runs/run*.jsonl` | 113 MB each at B_out = 750 (final), about 565 MB in total | **GitHub rejects files over 100 MB.** Options: Git LFS; a release asset or archive outside Git; or keep out of Git and rely on the hashes recorded in the end records and the verification report. |
| `reports/phase9/corrected_rerun/` (halts, gate, logs, manifest, sizing, verification at 500 and 750) | small | The audit trail of the corrected campaign. Recommended for commit together with whichever option is chosen for the checkpoints (`docs/LARGE_ARTIFACT_ARCHIVE.md`; manifest already generated). |
| `reports/phase9/dc_interval/` | 14 MB | Historical, seed-defective nested-bootstrap checkpoints. Pinned evidence for the deviation record and for V4. |
| `reports/phase9/execution/`, `execution_s7/` | small | Historical Phase 9 primary and S-7 results. Pinned inputs of the corrected manifest (θ̂ source). |
| `reports/phase9/final_project/` | 2.7 MB | The earlier provisional report package for the historical analysis. Its headline numbers are **not** final findings (see section 4). Contains a `__pycache__` directory, which `.gitignore` already excludes. |
| `data/processed/` | ignored | Derived from the public GEO series by scripts 13 and 23; pinned by hash in the corrected manifest. Correctly excluded. |

No credentials, tokens or private data were found. The tracked `metadata/GSE39582_sample_metadata.csv` is public GEO sample metadata.

## 3. Tests

**Current status: all tests pass.** Re-run at project closure, 8 November 2026: 1,085 passed, 0 failed, 0 skipped, 42 warnings (3 min 6 s). Earlier final run, 7 November 2026 (`.venv/bin/python -B -m pytest -q -p no:cacheprovider tests/`): **1,085 passed, 0 failed, 0 skipped** (42 warnings: NumPy runtime warnings inside deliberate non-finite tests, and the fork-in-a-threaded-process deprecation notice from the multiprocessing integration tests). No test is failing, skipped or excluded.

History: an earlier full run that day reported 3 failures. All three are closed:

| Test | Cause | Resolution |
|---|---|---|
| `test_phase9_corrected_d11_remediation.py::test_e1_control_stale_pyc_is_read_without_prefix` | Not a defect. The first full run set `PYTHONPYCACHEPREFIX`, which disables the in-tree bytecode cache that this negative control deliberately relies on. | Passes when the suite is run normally. No change. |
| `test_phase9_corrected_freeze_prepare.py::test_prepare_never_touches_the_real_repository` | Obsolete expectation: it asserted that the real freeze record is absent, the state before Decision 13. | The expectation now checks the same intent against the approved state: the real freeze record exists and still hashes to the approved Decision 13 anchor, so no test or prepare run created or altered it. |
| `test_phase9_corrected_nonestimable_stop.py::test_real_freeze_is_absent_and_verification_fails_closed` (renamed `…_is_the_approved_anchor_and_verification_fails_closed`) | Obsolete expectation, as above; with the record now present, the wrong-anchor refusal message changed from "missing or unreadable" to "does not match the approved Decision 13 hash". | The expectation now asserts that the installed record matches the approved anchor and that a configuration anchored to any other hash is still refused. The "not approved" refusal check is unchanged. |

Exactly two test functions were edited (the second and third rows); the first needed no change. No frozen code, configuration, protocol or result was changed.

Some tests read the locally held study inputs and historical outputs, which are not in Git, so they cannot pass on a fresh clone.

## 4. Documentation consistency

| Finding | Where | Recommended handling |
|---|---|---|
| The adjudication record is still headed "DRAFT — PENDING PROJECT-OWNER APPROVAL" with unticked boxes, while the configuration records owner approvals that cite it | `reports/phase9/final_project/documentation/PHASE9_METHODOLOGICAL_CLARIFICATION_AND_OWNER_ADJUDICATION_DRAFT.md` vs `config/…/PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json` (`approvals`) | Keep the draft unchanged as historical evidence; the approvals are recorded in the configuration with dates and references. A short dated approval note would remove the ambiguity. |
| The CV-seed deviation record is still headed "DRAFT — NOT YET ADJUDICATED" while `cv_seed_deviation_adjudicated` is approved | `reports/phase9/final_project/documentation/CV_SEED_DEVIATION_DRAFT_v2.md` | As above |
| The earlier final report and completion report describe the **historical as-run** nested bootstrap, labelled PROVISIONAL | `reports/phase9/final_project/documentation/PHASE9_FINAL_REPORT.md`, `PROJECT_COMPLETION_REPORT.md` | Keep as the historical record. The corrected report supersedes them for all findings. |
| The README described the pre-run status, a 453-sample cohort and 5- and 9-parameter models | `README.md` | Corrected on 8 November 2026: 452-patient A-19 population, 6 and 10 coefficients, verified result and limitations |
| The procedure documented a missing sizing tool | `docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md` section 7 | Resolved by `scripts/45` (commit `c432eee`) |

## 5. Unresolved decisions and technical debt

1. **Decision 12** (`lambda_at_failure`): the configuration entry is null and unapproved, although Option C was chosen and implemented (frozen `fitting.py`). It has no effect on this campaign, which had 0 fit failures, and the pipeline does not read it. Editing the configuration would break its hash binding, so record it separately (owner decision D3).
2. **Commit `c432eee`**: kept unchanged (owner decision D1, 8 Nov 2026). Not pushed.
3. **Large checkpoints**: choose how they are archived (section 2.2). The SHA-256 manifest is generated; the tar, the release and the `.gitignore` entry await approval.
4. **Two obsolete test expectations**: updated to the approved post-freeze state, and the suite passes in full (section 3). Only the review of the diff remains.
5. **Scripts 45 and 46 sit outside the Decision 13 freeze** by design. A future freeze could include them and give the launcher `size` and `report` commands.
6. **θ̂ alignment: audited and aligned** (`docs/PHASE9_THETA_ALIGNMENT_AUDIT.md`). θ̂'s authority is the plan (v2.1 line 6 and section 3.1), not "Decision 3": the owner's recorded Decision 3 is B_in,max and the RNG structural argument, and the θ̂ retention proposal (Decision 3 of the adjudication record) remains a DRAFT. θ̂ and the corrected outer replicates are the same functional on the same code and inputs, and the apparent component reproduces bit for bit. The differences are inner Monte Carlo precision and CV random-number assignment, plus one unmeasured apparent-CV link (strict analytical bound 1.6×10⁻³; typical size about 2×10⁻⁵). No additional analysis is needed. The owner may record an explicit approval of the disclosure wording. Re-estimating θ̂ under corrected streams would be a sensitivity analysis only and is not authorised.
7. **Calibration, S-7 and the five-set panel** were not re-run under the corrected streams. Their values come from the historical Phase 9 primary and S-7 runs and are reported only as such (final report, section 8).
8. **Not started, by instruction:** the batch-perturbation sensitivity analysis listed as a candidate in Phase 8D.
9. **Unexplained bootstrap shift.** A7 and A8 are triggered in all runs: the outer distribution sits above θ̂, opposite to the plan's section 7.2 anticipation. This is reported as a warning, as the plan requires. An investigation of its cause would be a new analysis and is not authorised.
10. **Draft adjudication records.** The θ̂ disclosure proposal (draft Decision 3) and the reporting-format proposal (draft Decision 9) were followed but never formally approved. Their numbering differs from the owner's recorded decisions.
