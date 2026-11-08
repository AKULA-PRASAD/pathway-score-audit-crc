# Release-readiness package

**Prepared:** 8 November 2026, from inspection of the working tree at HEAD `c432eee314c00dc174ded6c95e976eb41dbbe66d`. **Updated the same day** after the owner approved decisions A–E and the eight commit groups were committed locally (section 1.0). Nothing has been pushed, tagged, uploaded, moved or deleted, and the repository is **not publicly released**. No release version or tag is proposed; that is the owner's choice.

## 1.0 Repository state after the approved local commits (fresh inspection, 8 Nov 2026)

| Item | State |
|---|---|
| Branch | `phase9-corrected-pipeline`, tracking `origin/phase9-corrected-pipeline` |
| Remote tip | `6bb3f2673bed09fb1f516cb7f84a1d538abccce9`: **unchanged, nothing pushed** |
| Local commits on top of the remote | `c432eee` (kept, decision A), then G1 `7e0df408`, G2 `1c3bcaf`, G3 `5169b9f`, G4 `6300e39`, G5 `f35346d`, G6 `4ae202c`, G7 `3924526`, G8 `97ed2a04` (HEAD `97ed2a048a506e08bbdf30a3ace8e4620507173b` before this record), then one documentation-only follow-up commit containing this update |
| Working tree | clean: no modified or untracked files |
| Tracked files | 487 after G8 |
| Checkpoints | `reports/phase9/corrected_rerun/runs/run0–4.jsonl` on disk and unchanged; excluded by `.gitignore`; none in Git; archive manifest 52/52 OK |
| Largest tracked file | 6.3 MB (`reports/phase2/probe_statistics.csv`) |
| Still pending (separate authorisation) | push; building and publishing the checkpoint archive as a GitHub Release; code licence; public-release review, including per-sample tables |

## 1. Pre-commit state (as inspected before the commits; retained as the record the commit plan was based on)

| Item | State (checked on 8 Nov 2026) |
|---|---|
| Branch / HEAD | `phase9-corrected-pipeline` @ `c432eee`. Parent and `origin` tip `6bb3f267`. Ahead 1, behind 0 |
| Modified, tracked (6) | `.gitignore`, `README.md`, `config/analysis/PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json`, `requirements-analysis.txt`, `tests/test_phase9_corrected_freeze_prepare.py`, `tests/test_phase9_corrected_nonestimable_stop.py` |
| Untracked | 405 files. Exact list by commit group in `docs/release/COMMIT_PLAN_FILES.txt` |
| Frozen code | 27/27 files match `PHASE9_CORRECTED_FREEZE_RECORD.json` (`618f20eb…a790`). That hash equals the configuration anchor. The driver hash matches |
| Input manifest | 29/29 entries match (`1ff24387…3f27`) |
| Checkpoints | 5/5 match their end records and the archive manifest; never modified after the top-up |
| Archive manifest | `shasum -a 256 -c`: 52/52 OK |
| Frozen verification | `L00750.report.json` present: PASS 60/0/0/0, reporting permitted |
| Independent audit | present: INDEPENDENTLY REPRODUCED |
| Tests | **1,085 passed, 0 failed, 0 skipped, 42 warnings** (3 min 9 s, 8 Nov 2026, `-rs` reports no skips). The suite's temporary `scratch_p9c_*` directory was removed by the suite itself |

### Discrepancies found during this inspection, and their handling
1. **θ̂ RNG-link "bound".** Earlier text described the ≈ 2×10⁻⁵ figure as a bound. It is a typical-size (five-SD) figure. The strict analytical worst case is 4/2500 = 1.6×10⁻³. **Corrected** in the alignment audit, final report, inventory and register. Neither figure is an empirical measurement.
2. **Untracked count.** An earlier document said 407; the actual count was 401 before the closure documents were added, and 405 immediately before the commits. **Corrected.** All of those files except the five checkpoints are now committed.
3. **"Every pathway CI includes 1."** This is true for the primary four-set model only. The five-set sensitivity model's E2F interval is 0.746–0.996 (`docs/PHASE9_EXECUTION_REPORT.md` §5). **Corrected** after owner approval of G1–G13 (8 Nov 2026) in the final report, README and handover.
4. **Run-3 lower endpoint rounding** (−0.0091, not −0.0092). Already corrected on 8 Nov. Every quoted interval now matches the results file.

## 2. Package contents

| # | Item | Path | SHA-256 (first 16, current) |
|---|---|---|---|
| 1 | Final scientific report | `docs/PHASE9_CORRECTED_FINAL_REPORT.md` | `dc96b61029f2e670` |
| 2 | Independent audit report and results | `reports/phase9/independent_audit/INDEPENDENT_AUDIT_REPORT.md`, `independent_audit_results.json` | `abdea3df82285d8e`, `b67560a824182c45` |
| 3 | Final project handover | `docs/FINAL_PROJECT_HANDOVER.md` | `a71a2fbf236f2a37` |
| 4 | README | `README.md` | `e2862c0efbb35afa` |
| 5 | Reproduction instructions | `README.md` §Reproduction; `docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md` | `c0ac2982cfd9ec8b` (procedure) |
| 6 | Project completion inventory | `docs/PROJECT_COMPLETION_INVENTORY.md` | `10734aff4496bbe2` |
| 7 | Owner decision register | `docs/OWNER_DECISION_REGISTER.md` | `847ac3be997f3cba` |
| 8 | Archive manifest and restoration | `reports/phase9/ARCHIVE_MANIFEST_SHA256.txt`; `docs/LARGE_ARTIFACT_ARCHIVE.md` | `c45bb806ce844346`, `174102f48d2740f7` |
| 9 | Final test report | section 1 above | — |
| 10 | Commit plan (exact files) | `docs/release/COMMIT_PLAN_FILES.txt` | `29fec1bc71bdbb97` |
| — | Scientific interpretation review (gaps G1–G13, applied) | `docs/PHASE9_SCIENTIFIC_INTERPRETATION_REVIEW.md` | `6b502b79008ca264` |
| — | θ̂ alignment audit | `docs/PHASE9_THETA_ALIGNMENT_AUDIT.md` | `8d6e319a3833f8db` |
| — | Frozen verification at 750 | `reports/phase9/corrected_rerun/verification/L00750.report.json` | `9558aaa5dba467fb` |
| — | Official interval results | `reports/phase9/corrected_rerun_report/corrected_interval_results.json` | `1c504b137b177f2e` |
| — | Freeze record / interval plan v2.1 / modelling protocol v3.1 | `config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json`, `docs/PHASE9_DELTAC_INTERVAL_PLAN.md`, `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md` | `618f20eb26d41423`, `d27c52a29930a338`, `ab56708919f3746f` |

All hashes in this table were recomputed on 8 November 2026 after the G1–G13 corrections and the post-commit documentation update. This file is not listed, because a file cannot record its own hash.

**Superseded documentation hashes** (values in the first version of this table, before the G1–G13 corrections and the post-commit update; kept for traceability):

| File | Previous | Current |
|---|---|---|
| `docs/PHASE9_CORRECTED_FINAL_REPORT.md` | `11a7f5136010643d` | `dc96b61029f2e670` |
| `docs/FINAL_PROJECT_HANDOVER.md` | `76d3668415bc6660` | `a71a2fbf236f2a37` |
| `README.md` | `af61f3d7324bd02d` | `e2862c0efbb35afa` |
| `docs/PROJECT_COMPLETION_INVENTORY.md` | `4d74ca10544a7fa0` | `10734aff4496bbe2` |
| `docs/OWNER_DECISION_REGISTER.md` | `5102dc76f859833c` | `847ac3be997f3cba` |
| `docs/LARGE_ARTIFACT_ARCHIVE.md` | `02d32f40e1916b5a` | `174102f48d2740f7` |
| `docs/release/COMMIT_PLAN_FILES.txt` | `5ae8ca46015b8fed` | `29fec1bc71bdbb97` |
| `docs/PHASE9_SCIENTIFIC_INTERPRETATION_REVIEW.md` | `3f7a4b5eb9f02c55` | `6b502b79008ca264` |

Unchanged and not touched: the frozen code and freeze record, input manifest, configuration anchor, verification report, interval results, independent audit, θ̂ alignment audit, protocols, archive manifest and checkpoints.

## 3. Owner decisions A–E (summary; full evidence in `docs/OWNER_DECISION_REGISTER.md`)

| | Decision | Recommendation |
|---|---|---|
| A | Commit `c432eee` (3 files: `scripts/45`, its test, procedure doc). No frozen artifact touched. The working-tree copies are identical, and the sizing file it produced is part of the audit trail. Made before the no-commit rule | Keep it; push with the closure commits after approval. No reset, amend or revert |
| B | Archive of the five checkpoints (sizes and hashes in `docs/LARGE_ARTIFACT_ARCHIVE.md`) | GitHub Release asset plus the committed manifest; `.gitignore` the `runs/*.jsonl` and `*.lock` entries |
| C | Decision 12. Intended to resolve how `lambda_at_failure` is recorded when a fit fails (plan line 270). Option C (pass-through capture) is documented as the owner's choice in frozen `src/phase9_corrected/fitting.py` ("Owner Decision 12 (Option C, 2026-11-06)"); implemented and verified (D12-lambda PASS, 0 failed fits) | Proposed separate register entry (dated on the day of approval, **not backdated**): *"Decision 12 (lambda_at_failure): Option C — capture the original selected λ through a pass-through around the frozen `tune_lambda`, without additional calls — as documented in frozen `src/phase9_corrected/fitting.py` (header dated 2026-11-06). Implemented, tested and verified (verifier check D12-lambda PASS in all five runs; 0 failed fits, so no value was recorded). The configuration field `D12_lambda_at_failure` remains null because editing the configuration would break its hash binding to the completed campaign. Recorded on [date of owner approval]."* |
| D | Untracked Phase 1–8 and other files (classification in section 4) | Commit in groups G1–G8; never commit X or temporary files |
| E | θ̂ disclosure | Approve the wording below, in a separate dated note |

### Decision E: recommended disclosure wording
> *The point estimate θ̂ = +0.004316342804030704 is the pre-specified optimism-corrected ΔC on the full data (plan v2.1, line 6 and section 3.1). It was computed once by the Phase 9 primary execution from 2,500 event-stratified optimism resamples (5 seeds × 500) and was not re-estimated under the corrected random-number streams. Each outer replicate of the corrected nested bootstrap applies the same estimator with 100 inner resamples, so θ̂ and the interval share one functional and differ only in inner Monte Carlo precision; θ̂'s apparent component was reproduced bit-for-bit by the corrected engine. θ̂'s Monte Carlo standard error (≈ 0.00019) excludes the cross-validation randomness on the full data that selected λ. One historical random-number link — the apparent-fit cross-validation on the full data starts from the same generator state as one of the 2,500 optimism resamples — was not measured. Analytically its effect on θ̂ is at most 1.6 × 10⁻³ (strict worst case), and about 2 × 10⁻⁵ if that replicate is typical; these are calculations, not measurements. The primary percentile intervals do not use θ̂; the basic interval, z₀, the bias-corrected diagnostic, A7 and A8 do.*

The wording is checked against `docs/PHASE9_THETA_ALIGNMENT_AUDIT.md` §1, §4.1, §4.3 and §6, and against `scripts/38_phase9_execute.py` lines 90–92, 117–145.

## 4. Classification of untracked and modified files

| Class | Contents | Count | Disposition |
|---|---|---:|---|
| 1. Required project source | `src/` (9 modules), scripts 15–39 and 46–48, `tests/` (18 files plus fixture), configuration records, approved freeze record | 63 | Commit (G1, G3, G5, G7) |
| 2. Scientific documentation | `docs/PHASE1…PHASE9*`, closure documents, README | 40 | Commit (G4, G5, G8) |
| 3. Reproducibility and audit evidence | `reports/phase9/corrected_rerun/` minus the checkpoints, verification, sizing, gate, halts, manifests, independent audit, historical `dc_interval/`, `h5_spotcheck/` | ≈ 60 | Commit (G5, G6, G7) |
| 4. Generated outputs | `reports/phase1`–`phase8d` (152 files, 24 MB, largest 6.3 MB), `reports/phase9` historical outputs and the `final_project/` package | ≈ 230 | Commit (G4, G5). All under GitHub's 50 MB warning |
| 5. Temporary files | `reports/phase9/final_project/analysis/__pycache__/` (already ignored by `.gitignore`); test scratch directories (self-removed) | — | Never commit (already excluded) |
| 6. Large data and checkpoints | `reports/phase9/corrected_rerun/runs/run0–4.jsonl` (113 MB each); `data/raw`, `data/processed`, `data/external` (ignored) | 5 + ignored | Archive (B); never in Git |
| 7. Never publish | Ignored data directories (`data/external/marisa2013_tableS1.xls` is a journal supplement); any credentials (none found) | — | Already excluded |

**Per-sample content.** About 38 report files contain per-sample rows keyed by GEO accession (for example `reports/phase1/cohort_manifest.csv`, `reports/phase7c/sample_qc_metrics.csv`, `reports/phase8a/pathway_scores_*`). They are derived from **public, de-identified** GEO data, which already carry the same identifiers and clinical fields. Publishing them is consistent with GEO's open access. **The owner should confirm comfort with redistributing derived per-sample tables**, including the sex-discordance review `reports/phase8b/sex_discordant_eight_review.csv`.

**Local paths.** `/Users/kalki/…` occurs in 15 report files, 4 configuration records (including the frozen freeze record's environment fingerprint) and 5 historical script docstrings. They reveal a username and directory layout; they are not secrets. **They cannot be edited without breaking frozen or historical hashes**, so they are disclosed rather than changed.

**Secrets:** none (token, key and private-key patterns scanned). **Personal email:** appears only in Git commit metadata, which is standard.

## 5. Repository organisation review

**Recommendation: no file moves.**
- Every script resolves paths relative to the repository root, and many artifacts are pinned by path in manifests and freeze records. These include `config/…` in the freeze record and configuration anchor, and `reports/phase9/corrected_rerun/…` in the end records and verification.
- Renaming `config/` to `configs/`, or moving `reports/` or `scripts/`, would break the launcher's hash checks, the input manifest and the historical audit trail, and would invalidate the verification.
- No `notebooks/` directory exists, and none is needed.

| Observation | Recommendation |
|---|---|
| `results/probe_correlation_audit/` (tracked, Phase 0) alongside `reports/` | Keep; referenced by scripts 11–12 and the Phase 0 changelog |
| Protocols named `*_DRAFT.md` although approved (`PHASE6B_MODELING_PROTOCOL_DRAFT.md`, `PHASE9_MODELING_PROTOCOL_DRAFT.md`) | Keep the names; the Phase 9 protocol is pinned in the input manifest. The approval status is recorded in the config and docs |
| `reports/phase9/final_project/` (provisional historical report package) and its `PROJECT_COMPLETION_REPORT.md` / `PHASE9_FINAL_REPORT.md` | Keep unchanged as historical evidence. Supersession is stated in the final report header, README and inventory |
| Draft adjudication record and CV-seed deviation record still headed DRAFT | Keep unchanged; formal approval via separate notes (decisions C and E) |
| `docs/PHASE9_EXECUTION_REPORT.md` §15 says "Ready for final interpretation? No" | Historical; its three blockers are now resolved or reported (gap G11 proposes an explicit closure note in the final report) |

## 6. GitHub publication readiness

| Check | Result |
|---|---|
| README covers question, motivation, data, preprocessing, scoring, models, validation, results, audit, interpretation, limitations, reproduction, environment, layout, data access | Yes. Restrained conclusion; no robust-improvement, utility, absolute-risk, external-validation, causal or generalisation claim |
| File-size limits | Only the five checkpoints exceed 100 MB (excluded); largest committed file 6.3 MB |
| Broken references | None in the deliverable documents (shorthand such as `scripts/46` refers to numbered scripts) |
| Licence | None chosen (owner decision; required before public release) |
| Wording fixes G1–G13 | Applied after owner approval (8 Nov 2026) |

## 7. Proposed Git commit plan (not executed)

Exact files per group: `docs/release/COMMIT_PLAN_FILES.txt`. Groups apply in order, after decision A, by pushing `c432eee` along with them.

| Group | Purpose | Files | New / modified | Depends on | Sensitive or generated risk | State afterwards |
|---|---|---:|---|---|---|---|
| G1 | Approved freeze record, owner config record, environment, ignore rules (plus the B entry if approved) | 4 | 1 new, 3 modified | — | none | a clone can run the launcher's identity checks |
| G2 | Two test expectations updated to the post-freeze state | 2 | modified | G1 (tests read the freeze record) | none | test suite consistent with approved state |
| G3 | Phase 1–8 source, config records, tests, fixture | 52 | new | — | none | earlier-phase code and tests tracked |
| G4 | Phase 1–8 docs and outputs | 172 (24.5 MB) | new | G3 | per-sample derived public data (owner confirmation) | protocol trail and evidence tracked |
| G5 | Phase 9 primary, S-7 and historical as-run evidence | 119 (18.1 MB) | new | G3 | local paths in historical JSON | historical, superseded results preserved |
| G6 | Corrected campaign evidence without checkpoints, plus archive manifest | 38 | new | G1 | none | audit trail tracked; checkpoints pinned by hash |
| G7 | Report, θ̂ check, independent audit tools and outputs | 8 | new | G6 | none | final numbers and their audit tracked |
| G8 | Final documentation (README, report, audits, review, archive, inventory, register, handover, release package) | 11 | 1 modified, 10 new | G5–G7 | none | portfolio-ready documentation |
| X | Five checkpoints | 5 (565 MB) | — | — | GitHub limit | **not committed**; archived per decision B |

Commit messages should be professional English, with no tool attribution.

## 8. Final acceptance checklist

| Criterion | Status |
|---|---|
| Frozen verification PASS preserved | ✅ |
| Independent numerical audit preserved and satisfactory | ✅ INDEPENDENTLY REPRODUCED |
| Numerical outputs match the approved checkpoints | ✅ (138 quantities, 2.2×10⁻¹⁶) |
| Frozen code, input manifest and configuration hashes documented and re-verified | ✅ |
| Test suite current | ✅ 1,085 / 0 / 0 |
| Archive instructions match the files and manifest | ✅ 52/52 |
| Report, README, inventory and handover agree on all numbers | ✅ |
| Final report interpretation complete | ✅ G1–G13 applied after owner approval (8 Nov 2026) |
| Owner decisions A–E | ✅ approved 8 Nov 2026 (recorded in the decision register) |
| Licence | ⏳ open |
| Local commits | ✅ G1–G8 committed 8 Nov 2026 (`7e0df408`…`97ed2a04`), plus a documentation-only follow-up |
| Push, GitHub Release, archive upload | ⏳ not authorised |
