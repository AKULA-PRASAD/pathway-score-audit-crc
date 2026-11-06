# Phase 9 corrected pipeline — Decision 13 freeze procedure

**Status:** procedure only. It does not establish a freeze, approve Decision 13, or authorise any scientific execution (Decision 14).

## 1. Components

| Component | Path | Role |
|---|---|---|
| Scientific driver | `scripts/41_phase9_corrected_driver.py` | `validate` / `run` (unchanged by Decision 13) |
| Launcher | `scripts/42_phase9_corrected_launch.sh` | pre-Python identity checks; starts the driver (`validate`, `run`) or the runner (`verify`) |
| Post-run runner | `scripts/43_phase9_corrected_verify.py` → `src/phase9_corrected/postrun.py` | Stage A, Stage B, sizing recomputation, V1–V5, D12 statuses, verdict |
| Freeze preparation | `scripts/44_phase9_freeze_prepare.sh` | non-scientific bootstrap of the first freeze candidate |
| Capture entry point | `scripts/44_phase9_freeze_capture.py` | capture-only: environment and code inventory |
| Freeze record | `config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json` | absent until the owner approves a reviewed candidate |
| Owner anchor | owner configuration, `owner_decisions_outside_schema.D13_driver_freeze` | `sha256` (driver), `freeze_record_sha256` (freeze bytes), `approved`, `date`, `reference` |

## 2. Freeze record (schema `p9c-freeze-2`)

| Field | Content | Checked at launch |
|---|---|---|
| `schema` | `p9c-freeze-2` | yes |
| `code` | `{relative path: SHA-256}` for every `src/phase9_corrected/*.py` and every file in `orchestration.REQUIRED_FILES` | yes, in full (launcher `shasum -c`; Python re-hash; package file set must match exactly) |
| `environment` | blocking fields (below) | yes, exact equality |
| `environment_record_only` | `/usr/bin/sw_vers` output, `/usr/bin/jq --version` | no (shown in the runner report) |
| `bytecode_policy` | `fresh_compile_prefix` | yes |
| `launch`, `generated`, `status`, `excludes` | descriptive | no |

The freeze record never contains the owner configuration or its hash. The owner configuration anchors the freeze record, never the reverse.

`REQUIRED_FILES`:
- scripts 40 (historical driver, loaded by V1-inner), 41, 42, 43 and both 44 files;
- `src/__init__.py`, `src/paths.py`, `src/reproducibility.py`;
- `src/phase9_execution.py`, `src/pathway_scoring.py`, `src/score_standardisation.py`.

A missing, unexpected or modified file refuses execution.

**Blocking environment fields (P3):**
- Python implementation and version;
- resolved executable path and that file's SHA-256;
- `sys.prefix`, `sys.platform`, machine, and the OS family and architecture (`platform.system()`, `platform.architecture()`). The full `platform.platform()` string embeds the macOS point version, so it is recorded only: in the freeze (`environment_record_only.platform_string`) and in the run manifest (`record_only.platform_string`, schema `p9c-manifest-2`). The manifest's blocking platform identity is the same OS family and architecture, and `write_manifest_once` ignores only the `record_only` section, so a macOS point update alone neither blocks launch nor triggers H4 between campaign invocations. A change of OS family or architecture still does;
- installed distributions and versions;
- NumPy BLAS/LAPACK build identity;
- thread variables (`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `MKL_NUM_THREADS`);
- launch flags (`-I`, `-B`, and whether a pycache prefix was set).

**Record-only fields:** the exact macOS build (`/usr/bin/sw_vers`), the jq version and the full `platform.platform()` string.

**Thread variables must match.** `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS` and `MKL_NUM_THREADS` are blocking fields, read from the calling shell by both `scripts/44_phase9_freeze_prepare.sh` and the launcher. Set them identically (or leave them all unset) in the shell used for candidate capture and in every shell used for scientific launch. A difference refuses launch (`environment thread_env differs from the freeze`).

**Not authenticated:**
- shared libraries loaded by the interpreter (libpython, the BLAS/LAPACK binaries, system libraries);
- the kernel and its hardware;
- installed package files, unless the optional `distribution_records` field is added;
- code executing in memory after the checks.

These are source-file, metadata and launch-flag checks. They do not prove complete environment identity.

## 3. Initial-freeze bootstrap (P2)

The scientific launcher refuses (code 11) while no freeze record exists, so the first record is prepared without it:

1. **After implementation, testing and independent review, and only with explicit owner authorisation,** run:
   `scripts/44_phase9_freeze_prepare.sh /absolute/path/outside/the/repository/candidate.json`

   The script:
   - starts `.venv/bin/python -I -B -X pycache_prefix=<fresh, empty, external directory>` on the capture entry point only;
   - writes the candidate exclusively, outside the repository;
   - re-checks every inventory hash with `/usr/bin/shasum`;
   - prints the candidate's SHA-256.

   It never invokes the driver, the launcher, the gate, the runner or any scientific operation, and never writes inside the repository.
2. **The owner independently re-hashes the candidate** with `shasum -a 256 candidate.json`. The owner also re-hashes every file listed in `code`, and confirms the list is complete against `orchestration.REQUIRED_FILES` and the package directory.
3. **Independent review** of the candidate covers the inventory, the environment, the record-only fields and the bytecode policy.
4. **Only after that review, and with separate owner approval:**
   - copy the candidate byte-for-byte to `config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json`;
   - set `D13_driver_freeze.freeze_record_sha256` to the candidate's SHA-256;
   - set `sha256` to the driver's inventory entry;
   - set `approved`, `date` and `reference`.

   Any code change after step 1 invalidates the candidate; repeat from step 1.

The candidate is not a trust anchor. The owner's independent hashes in step 2 are the anchor.

## 4. Manual anchor step (before every execution campaign)

1. `shasum -a 256 config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json` must equal the approved `freeze_record_sha256`.
2. `shasum -a 256 scripts/42_phase9_corrected_launch.sh` must equal the launcher's entry in that freeze record.

The launcher cannot authenticate itself: a modified launcher could skip every check. Step 2 is therefore manual and is never delegated to the launcher.

## 5. Launch

- `scripts/42_phase9_corrected_launch.sh validate|run <owner-config> [args]` starts the driver.
- `scripts/42_phase9_corrected_launch.sh verify <owner-config>` starts the post-run runner.

Both go through the same pre-Python checks and the same `-I -B -X pycache_prefix` flags. The launcher checks `approved` but not `date` or `reference`; the Python `verify_code_identity` checks all three. In-process, `verify_code_identity` repeats the checks and verifies every loaded repository-local module against the freeze, including the historical driver alias `p9_dc_driver_frozen_ro`. A short time-of-check/time-of-use window remains.

**Repository-root shadowing (review D1).** Every entry point puts the repository root first on `sys.path`. A module file or regular package planted at the root could therefore shadow a standard-library or third-party module that is imported lazily, after the identity check (for example `multiprocessing` in `execute_replicates`, or `lifelines` on the first real refit).
- The runner re-checks loaded-module identity, and refuses importable root entries other than `src` and `tests`, immediately before every output it writes. A shadowed process therefore cannot write verification evidence.
- This detects a shadow only after it was imported: it does not authenticate code before that code runs.
- The scientific driver's own lazy imports are not re-checked, because `execute_replicates`, the driver and the launcher are outside the authorised scope. A pre-launch root-shadow refusal there needs separate authorisation.

## 6. Runner verdicts (P4)

- Each check reports `PASS`, `FAIL`, `NOT_APPLICABLE` or `NOT_COMPARABLE`, with counts and reasons.
- Missing, malformed or empty evidence is `FAIL`.
- `NOT_APPLICABLE` is accepted only for V4 when the approved configuration has `as_run_variant = false` under the independent convention.
- Any `NOT_COMPARABLE` makes the verdict `OWNER_REVIEW_REQUIRED`. The only owner-approved `NOT_COMPARABLE` classification is the V4 skipped-record policy (D2, below); it still requires owner review.
- Only `PASS` permits reporting.
- V5 `FAIL` is hard failure H5.

Outputs (window results and reports) are created exclusively under `output_dir` and never replaced. Executing the runner performs real-data refits (Stage A, V1-inner, V5) and requires Decision 14.

**Pre-execution gate evidence (gate-cache remediation).** An adversarial test showed that the original shared cache (`gate/pre_execution_gate.json`, reused whenever its key matched and `all_pass` was true) let a failed gate be edited into a PASS and authorise dispatch. That cache is no longer used.
- **Fresh execution, every attempt.** Before any data is loaded, any authorisation is consumed or any start record is written, `run_corrected` executes the complete six-check gate into a new file, `<output_dir>/gate/run<r>.inv<mmmm>.attempt<aa>.json`. Because the file does not exist beforehand, `gate.py` never finds or reuses a cache.
- **PASS established from the written evidence.** Dispatch requires `orchestration.gate_evidence_problems` to find no problem in the file:
  - exactly the six approved checks, all PASS;
  - the key equals this configuration, manifest and NumPy version;
  - the approved V2-domain domain and state count;
  - V1-outer-historical covering 500 replicates for every run.

  A returned `all_pass` flag alone never authorises dispatch.
- **Failure.** A failed or malformed attempt keeps its evidence file, logs `gate_failed` and raises `SystemExit` with the reasons. It writes no start record, consumes no authorisation and dispatches nothing. The next attempt uses a new file.
- **Data loading.** Data is loaded only after a genuine gate PASS and before any authorisation is consumed. If loading fails, a pending runtime-halt or unterminated-invocation authorisation remains unused and the same reference can be supplied again.
- **Binding.** The start record holds `gate_evidence = {path, sha256, all_pass, attempt}`.
- **Re-verification.** `assess_run` re-verifies the evidence of every invocation, so dispatch, Stage B and sizing all depend on it:
  - a missing reference or file is `blocking` (`gate_evidence_unverifiable`);
  - a SHA-256 mismatch or content that is not a full PASS is a `disqualifier` (`gate_evidence_invalid`).

  The runner also reports one `gate-evidence` result per invocation.
- **Bound to its own invocation.** The evidence path must be `gate/run<r>.inv<mmmm>.attempt<aa>.json` for the start record's own run, invocation and attempt. A reference to another invocation's (genuine) evidence is a `disqualifier`.
- **Malformed evidence.** Validation never raises. Wrong types, excessive nesting, unreadable or invalid paths and unparseable JSON are reported as reasons (fail closed).
- **If `gate.py` itself raises** (an error inside a check, not a FAIL result), the exception propagates. Nothing is dispatched and no start record or authorisation is written, but no evidence file and no `gate_failed` event are left. The owner diagnoses the error from the process output; a later attempt runs the gate afresh.
- **Interrupting a gate run.** The gate runs while the run's checkpoint lock (`runs/run<r>.lock`) is held. Killing the process during the long V2-domain enumeration (about 4.5M states) leaves that lock behind. Further invocations then refuse ("another writer is active or a previous ... "), and the post-run runner refuses on the incomplete artifact. After confirming that no process is running, the owner removes the stale lock manually; no other record needs repair.
- **Limitation.** The evidence files are hashed and bound to start records, but not cryptographically authenticated. A coordinated rewrite of an evidence file *and* its start record, by someone with write access to `output_dir`, cannot be detected by hashing alone; the same holds for the other create-only records. Skipping the gate itself would require changing frozen code, which the freeze detects. The `gate.py` docstring still describes the former single-file cache; `gate.py` is protected and unchanged.

**V4 skipped records (review D2): owner-approved NOT_COMPARABLE.**

| V4 condition (k < 500, the historical comparison domain) | Classification |
|---|---|
| historical row missing | NOT_COMPARABLE |
| historical status not `ok` | NOT_COMPARABLE |
| corrected as-run fit failed | NOT_COMPARABLE |
| corrected as-run comparison field missing | NOT_COMPARABLE |
| valid historical and corrected pair | compared; a mismatch is FAIL |
| k ≥ 500 | outside the comparison domain; counted separately (`n_beyond_historical_domain`) |

- Every skip is counted, with the run, the k positions and the reason, in the result `V4:not_comparable`. A skipped record is never a successful comparison and never a scientific FAIL.
- If every in-domain record is skipped, V4 itself is NOT_COMPARABLE.
- The runner cross-checks V4's comparison count against its own skip accounting; a disagreement is FAIL.
- An error in the accounting makes V4 NOT_COMPARABLE, with the exception recorded. It is never PASS.

**V5 (reviews D4–D6).**
- The approved `!=` comparison is kept. Under it:
  - a top-level NaN field always differs (`NaN != NaN`), so it FAILs and is H5;
  - NaN nested in a list or dict compares equal, because json returns one shared NaN object and container comparison checks identity first;
  - ±Infinity compares equal.

  The report classifies every V5 field's non-finite content (`top_level_nan`, `nested_nan_in_list`, `nested_nan_in_dict`, `other_nonfinite`). It also lists the differing fields whose canonical JSON is identical, i.e. those that differ only because of a top-level NaN. The result is unchanged.
- V5 uses run 0 (plan line 383). If run 0 is not selected, V5 is NOT_COMPARABLE and the verdict is OWNER_REVIEW_REQUIRED; no other run is substituted.
- An exception during V5 recomputation is FAIL, with its type, message, run and stage recorded. It is not H5.

**Campaign artifacts (review D5).** Every file in `runs/` and `halts/` is classified:
- complete: `run<r>.jsonl`, and `run<r>.inv<mmmm>.<start|end|auth>.json`;
- incomplete: `run<r>.lock` (an active writer, or a stale lock left by a killed one), `run<r>.jsonl.tmp`, and `*.json.tmp`;
- unrecognised: any other name beginning with `run`;
- unrelated: anything else, which is disclosed and never counted as a run.

Verification refuses if any incomplete or unrecognised artifact is present, if a selected run has no artifacts, or if any artifact belongs to an unselected run.

**Root-shadow detection (review D1).** The root scan reports top-level module files, extension modules and regular packages, including packages whose `__init__` is a versioned (`.cpython-*`) or `abi3` extension. A namespace directory (no `__init__`) is not reported, because it cannot shadow a regular module found later on `sys.path`.

## 7. Path C: the missing sizing-decision step (owner decision outstanding)

A top-up beyond B_out_initial (Path C) is dispatched only when `run_corrected` verifies `<output_dir>/sizing/level<B_cur, 5 digits>.json` (schema `p9c-sizing-1`, `verify_sizing_decision`). **No tool creates that file.** `build_sizing_decision` computes its content, but nothing writes it. The runner verifies sizing files and must never create them (Decision 13 B3).

A future, separately authorised sizing-decision workflow must:
1. Start through the launcher, under the same freeze, flags and identity checks.
2. Take as inputs only: the approved configuration; the manifest; the verified freeze; and the five checkpoints at B_cur, each run complete at B_cur, CLEAN, ELIGIBLE, and with every window a bound Stage A PASS. Run the post-run runner first at that level to produce the window results.
3. Compute the content with `orchestration.build_sizing_decision`, the single source of truth (Decisions 5A, 5B, 6A, 6B, 7A and 7B), and write it once with `durable_create`. It must never overwrite an existing file.
4. Write only when the decision is `top_up`. Terminal or owner-review decisions are reported, not written as authorisations.
5. Stay separate from verification. It must not dispatch a top-up (each run's top-up still needs the owner's explicit `--topup-authorization` reference through `run_corrected`), and the runner must keep re-deriving every sizing file independently.

## 8. Open items

- **Sizing files:** see section 7 (owner decision).
- **Platform:** Linux paths are untested. The procedure assumes macOS (`/usr/bin/jq`, `/usr/bin/shasum`, `/usr/bin/sw_vers`, `F_FULLFSYNC`).
