# Does a pathway signature add prognostic value in colon cancer? A fully auditable nested-bootstrap analysis

> **Status (2026-11-06):** the corrected Phase 9 pipeline is **implemented, tested and committed**. A **candidate code freeze** has been generated and checked; final freeze approval is pending. The previously approved owner decisions (1–11) have been entered into the run configuration. **The corrected scientific bootstrap has not been run.** No corrected confidence interval or scientific conclusion exists yet. All results below are labelled **historical** or **pending**.

## 1. The question

Gene-expression signatures are often reported to improve recurrence prediction in colorectal cancer. Many such claims do not survive careful validation: the apparent gain is inflated by **optimism**, because a model with more parameters fits its own training data better than new data.

This project asks a narrow, pre-specified question on a public cohort: **does adding four frozen pathway scores to a clinical Cox model improve discrimination after correcting for optimism, and how precisely can that improvement be estimated?**

| Item | Specification |
|---|---|
| Data | GEO **GSE39582** (Affymetrix GPL570), stage II–III colon cancer, cohort **S3 = 453 samples**, 0 exclusions |
| Endpoint | cause-specific recurrence, as deposited |
| M1 | clinical-only Cox model (5 parameters) |
| M2 (primary) | M1 + four pre-registered pathway scores, frozen panel (9 parameters) |
| Estimand | **ΔC_corrected = C_corrected(M2) − C_corrected(M1)**: the optimism-corrected difference in Harrell's C |
| Deliverable | a 95 % percentile interval for ΔC_corrected, with Monte-Carlo error reported |

Everything was fixed before any outcome analysis: the cohort, endpoint, panel membership, scoring, representation (median, primary) and model forms. Each decision is recorded in a dated, hash-pinned protocol (`docs/PHASE9_MODELING_PROTOCOL_DRAFT.md`, `docs/PHASE9_DELTAC_INTERVAL_PLAN.md`).

## 2. Why a nested bootstrap

A single Harrell optimism correction gives a **point estimate** of ΔC_corrected but no valid interval. An interval needs the sampling distribution of the *whole correction procedure*:

```
for each outer replicate k = 1..B_out          (B_out = 500, top-up ladder 750 … 3000)
    D*_k  = stratified bootstrap of the cohort
    fit M1, M2 on D*_k                          -> apparent C on D*_k
    for each inner replicate b = 1..B_in        (B_in = 100)
        J*_k,b = bootstrap of D*_k
        fit M1, M2 on J*_k,b; evaluate on D*_k   -> optimism(M, k, b)
    ΔC_corrected(k) = ΔC_apparent(k) − mean_b[opt(M2) − opt(M1)]
interval = percentiles of {ΔC_corrected(k)}
```

M2's ridge penalty λ is re-tuned by 5-fold CV inside **every** fit (about 10⁵ penalised Cox fits per run). Design choices that are argued in the plan rather than assumed:
- **M1 and M2 are paired** on the same resamples. Pairing removes about 44 % of the inner Monte-Carlo noise (ρ ≈ 0.68).
- **B_in = 100** inflates the outer SD by about 1.2 %. The inflation is conservative and is disclosed, not deconvolved.
- Inner replicates where either model fails are **excluded as pairs**, with a bounding analysis.
- **Five independent full-size runs** measure Monte-Carlo reproducibility directly.
- **Endpoint Monte-Carlo SE** is estimated by kernel quantile methods at three bandwidths. EP-1 uses the most conservative; reliability screening uses the least.
- **Adaptive sizing:** B_out grows along a fixed ladder until the endpoint MC-SE is below 5 % of the interval width, or a W5 stability breach forces a top-up.

## 3. Why there is a "corrected" pipeline

The original Phase 9 run (historical, 2026-11-03/04) produced a point estimate, **ΔC_corrected = +0.0043**, but **no interval for it**. A later audit found a protocol deviation: every λ-tuning cross-validation consumed **the same PCG64 starting state as the bootstrap resample it was tuning on**, instead of the independent child stream the plan specified. This was verified empirically, 9/9 cases (`CV_SEED_DEVIATION_DRAFT_v2.md`).

The historical outputs are **preserved byte-for-byte and hash-pinned**. They are used only as verification references. The corrected pipeline (`src/phase9_corrected/`) re-implements the procedure with the specified stream separation and reuses the frozen fitting code (`src/phase9_execution.py`) unchanged. It reproduces the historical *outer* resamples exactly, which is checked against stored hashes.

## 4. Engineering for reproducibility and auditability

### RNG stream design
- Every random draw comes from a **named `SeedSequence` child**: outer, apparent-CV (with its own entropy word `T`), inner, inner-CV and inner-reserve, per run, k and b. No stream is ever reused.
- A **pre-execution gate** enumerates every initial PCG64 state the run could use. For the approved domain that is about 4.5 M states (5 runs × k < 3,000 × b < 100 × 5 families, plus the primary seeds), and the gate proves they are all distinct.
- Every replicate record stores SHA-256 fingerprints of the stream states and resample multiplicities it consumed, so stream assignment is **verifiable after the fact**.

### Execution control
- **Append-only checkpoints** with a header binding config and manifest hashes. Each invocation's durable `start`/`end`/`auth` records seal byte boundaries, and history is re-derived from them before any dispatch.
- **Crash consistency:** create-only writes with `F_FULLFSYNC`. An interrupted invocation is detected as UNTERMINATED and needs an explicit owner acknowledgement before it resumes.
- **Runtime monitor:** projects completion time after warm-up and halts at twice the declared budget. Resuming needs a single-use authorisation reference.
- **Non-estimable stop rule:** a strict literal 5 % predicate (`20·b > n`) applied to every prefix of the record stream.
- **Eligibility engine:** classifies every run as INELIGIBLE, UNVERIFIABLE, INCOMPLETE or ELIGIBLE, with all reasons retained. The same assessment drives dispatch, sizing and reporting.
- **Top-ups (Path C):** allowed only from a sizing decision that the frozen code re-derives from all five checkpoints.

### Code identity: "the code that runs is the code that was reviewed"
- A **freeze record** holds SHA-256 hashes of 27 source files and a blocking environment fingerprint: interpreter path and hash, packages, the NumPy BLAS/LAPACK build, thread variables and launch flags.
- The owner anchors it by recording its SHA-256 in the configuration.
- A **bash launcher** (`scripts/42`) verifies the anchor and every file with `/usr/bin/jq` and `/usr/bin/shasum` **before Python starts**. It then runs `python -I -B -X pycache_prefix=<fresh empty dir>`, so no stale bytecode can be loaded.
- In-process checks confirm that every loaded repository module matches the freeze.
- The run manifest and the gate evidence bind every invocation to the same configuration, code and data hashes.

### Post-run verification runner (`scripts/43`)
The runner works through these steps in order:
1. Identity binding.
2. **Stage A:** recomputes every interrupted window and compares it byte-for-byte.
3. **Stage B:** eligibility.
4. Sizing recomputation at every ladder level.
5. **V1–V5:**
   - V1: outer and inner resamples regenerated and compared with stored hashes and with the frozen historical code path;
   - V2: stream assignment;
   - V3: invariance of the M1 quantities to the CV correction;
   - V4: reproduction of the as-run variant against history;
   - V5: full recomputation of 20 replicates.

Verdicts are PASS, FAIL, NOT_APPLICABLE or NOT_COMPARABLE, each with counts and reasons. Missing evidence is never PASS, and only an overall PASS permits reporting.

## 5. Verification of the software (not the science)

- **1,075 automated tests** pass on synthetic data; 549 of them cover the corrected pipeline.
- **Independent oracles:**
  - sizing and MC-SE rules are re-derived in a pure-Python reference implementation, at 1e-10 tolerance;
  - the gate's state-count formula is checked against the real enumerator;
  - non-finite JSON semantics are checked empirically.
- **Adversarial tests**, for example:
  - a failed pre-execution gate edited to PASS (it now cannot authorise dispatch, because the gate runs freshly per invocation);
  - tampered or forged evidence files;
  - stale `.pyc` shadowing;
  - repository-root module shadowing;
  - launcher refusal paths.
- **Process:** owner-approved decisions, each implementation independently code-reviewed, findings triaged and fixed, protected files verified by hash after every cycle.

Passing tests show the software implements the approved protocol. They are **not** evidence about colon-cancer biology or about the size of ΔC.

## 6. Status

| Component | Status |
|---|---|
| Protocol, estimand, cohort, models | Frozen (historical phases) |
| Historical point estimate ΔC_corrected = +0.0043 | Historical; **no interval**; affected by the CV-seed deviation |
| Corrected pipeline, launcher, runner, freeze tooling | **Implemented, tested, committed** (`bba7f6d`) |
| Candidate freeze record | **Generated and checked against the committed code; final approval (Decision 13) pending** (27 files, SHA-256 `618f20eb…a790`) |
| Owner decisions 1–11 in the run configuration | **Populated** from the recorded approvals, with dates and references; the configuration passes the pipeline's static validation except for the execution authorisation |
| Decision 12 (λ at failure, Option C) | **Implemented and tested**; its entry in the configuration's decision record is **unresolved**, because the template's allowed values do not include Option C (the pipeline does not read this entry) |
| Decision 13 final freeze approval | **Pending** |
| Decision 14 scientific execution authorisation | **Pending** |
| Corrected bootstrap (5 runs × 500 × 100) | **Not run** |
| Corrected interval for ΔC_corrected, and any corrected conclusion | **Pending**: no corrected result exists |

## 7. Remaining steps to a scientific result

The run configuration (`config/analysis/PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json`) already holds the approved scientific decisions. What remains:
1. **Decision 13:** the owner independently re-hashes and reviews the candidate freeze, installs it as the freeze record, and records its SHA-256 and the driver hash.
2. **Manual anchor check:** re-hash the freeze record and the launcher; confirm the repository root has no importable entries other than `src` and `tests`.
3. **Decision 14:** the owner authorises execution in the configuration.
4. Run `scripts/42_phase9_corrected_launch.sh run …` for runs 0–4 (with the four BLAS/OpenMP thread variables unset, as at freeze capture), then `… verify`.

## 8. Repository map

| Path | Contents |
|---|---|
| `src/phase9_corrected/` | engine, RNG streams, resampling, fitting adapter, checkpointing, orchestration, statistics (EP-1, W5, sizing), reliability, gate, verification, post-run runner |
| `src/phase9_execution.py`, `src/pathway_scoring.py`, `src/score_standardisation.py` | frozen historical fitting and scoring code (unchanged) |
| `scripts/40` | historical driver (frozen reference) |
| `scripts/41`–`44` | corrected driver, launcher, verification entry point, freeze preparation |
| `tests/test_phase9_corrected_*.py` | synthetic, adversarial and oracle tests |
| `docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md` | freeze, trust anchor and launch procedure |

```bash
.venv/bin/python -m pytest -q tests/test_phase9_corrected_*.py   # synthetic tests only; no real data
```

## 9. Skills demonstrated

Survival analysis (Cox models, Harrell's C, optimism correction) · resampling inference (nested bootstrap, Monte-Carlo error, kernel quantile SEs) · reproducible-research engineering (RNG stream design, content-addressed code freezes, crash-consistent checkpoints) · verification design (independent oracles, adversarial tests) · protocol discipline (pre-registration, deviation records, honest scoping of claims).
