# Owner decision register — project closure

**Prepared:** 8 November 2026. These decisions belong to the project owner. None has been taken; nothing has been committed, pushed, amended, reset, published or uploaded.

**Cross-cutting constraint.**
- The configuration hash `bf939a4e…6a56` covers the whole run configuration except `n_workers`, and every checkpoint header, invocation record and verification report is bound to it.
- **Any edit to `config/analysis/PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json` would break that binding** for the completed campaign.
- Decisions that need recording should therefore go in a separate dated record (D3–D5 below), not in the configuration.

## D1. Unapproved local commit `c432eee`

| | |
|---|---|
| Decision | Keep, push, or leave unpushed and replace by an owner-made commit. History may not be rewritten |
| Current state | `c432eee314c0…` on `phase9-corrected-pipeline`, parent `6bb3f267…` (the pushed tip). Ahead of `origin` by 1; **not pushed**. Author Prasad Akula, 7 Nov 2026 12:19 −0500. Professional message with no attribution trailer |
| Contents | `scripts/45_phase9_sizing_decision.py` (+143), `tests/test_phase9_corrected_sizing_writer.py` (+74), `docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md` (+21/−4) |
| Working-tree implications | The three files are byte-identical in the working tree; nothing was modified after the commit. No frozen file, configuration or result is in the commit. `scripts/45` (SHA-256 `cf4e3702…`) is the tool that wrote `sizing/level00500.json`, the decision the owner-authorised top-up relied on. The frozen verifier re-derived that decision independently (sizing at 500: PASS) |
| Violation? | It was made **before** the instruction "Do not commit, push, publish, or modify Git history" was given. It breaches no rule in force at the time, but it was not explicitly approved. It has been disclosed since |
| Recommended option | **Keep the commit as is and push it with the rest of the closure work**, once the owner approves. It is the provenance of a file in the audit trail. Removing it would need a history rewrite (prohibited), or a revert that leaves the sizing file without its tool in history |
| Consequences | Keep: the history reflects what happened. Revert commit: the tool disappears from the tip while its output remains; not recommended |
| Owner approval required | **Yes**, for any push and for the disposition |

## D2. Archive destination for the five checkpoints, and the `.gitignore` entry

| | |
|---|---|
| Decision | Where the five 113 MB checkpoints live (Release asset, Zenodo or Git LFS), and whether to add the documented `.gitignore` entry |
| Current state | Checkpoints in place, unchanged. Manifest `reports/phase9/ARCHIVE_MANIFEST_SHA256.txt` (52 files, all OK). No tar built, nothing uploaded, `.gitignore` unchanged |
| Evidence | `docs/LARGE_ARTIFACT_ARCHIVE.md`: exact sizes, hashes, directory structure, restore and verify commands. GitHub rejects files over 100 MB |
| Recommended option | **GitHub Release asset** (deterministic uncompressed tar) plus the committed manifest, and ignore only `reports/phase9/corrected_rerun/runs/*.jsonl` and `*.lock`. Optionally add Zenodo for a DOI |
| Consequences | Release asset: free, outside history, verifiable through the manifest. Git LFS: a seamless clone but a 1 GB quota, and each clone downloads about 0.6 GB |
| Owner approval required | **Yes**: upload, destination and `.gitignore` |

## D3. Blank Decision 12 entry (`D12_lambda_at_failure`)

| | |
|---|---|
| Decision | How to record the Decision 12 outcome in the decision log |
| Current state | The configuration entry has `value: null`, `approved: false`. Its template lists only `accept_not_recorded` and `approve_additional_tune_lambda_call` |
| Intended decision | The owner chose **Option C** (2026-11-06): capture the original λ at failure through a pass-through around the frozen `tune_lambda`, without extra calls. This is recorded in `src/phase9_corrected/fitting.py` (header "Owner Decision 12 (Option C, 2026-11-06)"), implemented, tested, independently audited before Decision 13, and inside the freeze |
| Evidence for completion | Verifier check `D12-lambda`: PASS in all 5 runs, over 153,000 fits per run, with 0 failed fits, so no λ-at-failure value ever had to be recorded. The pipeline does not read the configuration entry |
| Recommended option | **Do not edit the configuration** (cross-cutting constraint). Record a dated owner note (for example `docs/OWNER_DECISION_12_RECORD.md`) stating "Option C approved 2026-11-06; implemented in frozen fitting.py; no failures occurred". Reference it from the inventory |
| Consequences | Completes the decision log with no effect on the campaign. Editing the configuration instead would invalidate the configuration binding of every checkpoint and verification report |
| Owner approval required | **Yes**: the wording and the confirmation that Option C was the decision |

## D4. Committing the Phase 1–8 (and other untracked) files

| | |
|---|---|
| Decision | Whether to commit the untracked project history, and in which grouping |
| Current state | 401 untracked files and 6 modified tracked files (`git status`). Includes 25 docs, scripts 15–39 and 46–48, 9 `src` modules, 17 test files, configuration records, the approved freeze record, `reports/phase1`–`phase8d` (about 25 MB) and `reports/phase9` |
| Evidence | `docs/PROJECT_COMPLETION_INVENTORY.md` section 2. Without the freeze record a clone cannot run the launcher. The tracked README and protocol cite several untracked docs. No credentials found. Some historical JSON and docstrings contain the local path `/Users/kalki/…`, which is not sensitive but reveals a username; the frozen and historical files cannot be edited without breaking their hashes |
| Recommended option | **Commit everything except the five checkpoints**, in a few logical commits: (1) the approved freeze record and configuration change; (2) Phase 1–8 docs, scripts, `src`, tests and reports; (3) Phase 9 historical outputs; (4) the corrected campaign evidence, final report, independent audit, archive manifest and README |
| Consequences | The repository becomes reproducible and auditable from a clone, plus the archive. Not committing leaves the public record without its protocol trail |
| Owner approval required | **Yes**, for all commits and pushes |

## D5. Approval of the θ̂ disclosure wording (and the reporting format)

| | |
|---|---|
| Decision | Whether to approve formally two proposals of the draft adjudication record that were followed: Decision 3 (retain θ̂ with disclosures) and Decision 9 (per-run intervals primary; pooled draws descriptive only) |
| Current state | Both are unticked in `reports/phase9/final_project/documentation/PHASE9_METHODOLOGICAL_CLARIFICATION_AND_OWNER_ADJUDICATION_DRAFT.md`. Their numbering differs from the owner's recorded decisions (where Decision 3 = B_in,max and Decision 9 = verification gate) |
| Evidence | `docs/PHASE9_THETA_ALIGNMENT_AUDIT.md`: θ̂ is the plan-specified estimator (plan v2.1, line 6 and section 3.1), the same functional the interval bootstraps; its apparent component reproduces bit for bit; one unmeasured RNG link has a strict analytical bound of 1.6×10⁻³ (typical size about 2×10⁻⁵). Section 6 of that audit holds the five disclosures used in the final report. The plan itself designates the per-run percentile interval as primary and no pooled interval |
| Recommended option | **Approve both, in a separate dated owner note.** Keep the draft unchanged as historical evidence |
| Consequences | Closes the governance gap; no change to any number. Rejecting would require only rewording, because the percentile interval does not use θ̂ |
| Owner approval required | **Yes** |

## Minor, optional

- **Code licence:** the repository has no licence file. Choose one if the repository is to be public. Owner decision.
- **Career dossier artifact:** shows the pre-result status. It can be refreshed from the final report on request.

---

## Owner approvals recorded — 8 November 2026

On 8 November 2026 the project owner approved decisions D1–D5 as recommended above and authorised the local commit plan (`docs/release/COMMIT_PLAN_FILES.txt`). Pushing, GitHub Release creation, archive upload, publication and the choice of a code licence were **not** authorised and remain open.

| Decision | Recorded outcome |
|---|---|
| D1 (A) | Commit `c432eee` kept unchanged in history; no amend, reset, rebase, revert or cherry-pick. Push not yet authorised |
| D2 (B) | Archive destination: GitHub Release asset with the committed SHA-256 manifest. `.gitignore` now excludes `reports/phase9/corrected_rerun/runs/*.jsonl` and `*.lock`. The archive has not been built or uploaded |
| D3 (C) | Decision 12 record, below |
| D4 (D) | Untracked project files committed locally in groups G1–G8; the five checkpoints excluded |
| D5 (E) | θ̂ disclosure wording, below; the reporting format (per-run percentile intervals primary, pooled draws descriptive only) is approved |

**Decision 12 (lambda_at_failure), recorded 8 November 2026.** Option C: capture the original selected λ through a pass-through around the frozen `tune_lambda`, without additional calls, as documented in frozen `src/phase9_corrected/fitting.py` (header dated 2026-11-06). Implemented, tested and verified: verifier check D12-lambda PASS in all five runs, with 0 failed fits, so no value was recorded. The configuration field `D12_lambda_at_failure` remains null, because editing the configuration would break its hash binding to the completed campaign. This record is dated on the day of approval and is not backdated.

**Approved θ̂ disclosure, recorded 8 November 2026.**

> The point estimate θ̂ = +0.004316342804030704 is the pre-specified optimism-corrected ΔC on the full data (plan v2.1, line 6 and section 3.1). It was computed once by the Phase 9 primary execution from 2,500 event-stratified optimism resamples (5 seeds × 500), and was not re-estimated under the corrected random-number streams.
>
> Each outer replicate of the corrected nested bootstrap applies the same estimator with 100 inner resamples. θ̂ and the interval therefore share one functional and differ only in inner Monte Carlo precision. θ̂'s apparent component was reproduced bit for bit by the corrected engine.
>
> θ̂'s Monte Carlo standard error (≈ 0.00019) excludes the cross-validation randomness on the full data that selected λ.
>
> One historical random-number link was not measured: the apparent-fit cross-validation on the full data starts from the same generator state as one of the 2,500 optimism resamples. Analytically its effect on θ̂ is at most 1.6 × 10⁻³ (strict worst case), and about 1.9 × 10⁻⁵ if that replicate is typical. These are calculations, not measurements.
>
> The primary percentile intervals do not use θ̂. The basic interval, z₀, the bias-corrected diagnostic, A7 and A8 do.
