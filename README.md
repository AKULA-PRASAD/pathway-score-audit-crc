# Optimism-corrected prognostic value of a pathway signature in colon cancer

A pre-specified analysis of whether four pathway scores improve recurrence prediction beyond a clinical Cox model, with a nested bootstrap to put a confidence interval on the optimism-corrected gain.

**Status (6 November 2026).** The corrected Phase 9 pipeline described below is implemented, tested and committed. A candidate code freeze has been generated and checked against the committed code but has not yet been approved. The previously approved analysis decisions are entered in the run configuration. **The corrected bootstrap has not been run.** No corrected interval or conclusion exists yet; any number quoted here is either historical or explicitly marked as pending.

## The question

Expression signatures are frequently reported to improve prognosis in colorectal cancer, but part of the apparent improvement is usually optimism: a model with more parameters fits its own training data better than it fits new patients. I wanted to ask a narrow question on a public cohort and answer it carefully. After correcting for optimism, does adding a small, frozen set of pathway scores to a clinical model improve discrimination, and how precisely can that improvement be estimated?

| | |
|---|---|
| Data | GEO GSE39582 (Affymetrix GPL570), stage II–III colon cancer; analysis cohort of 453 samples, no exclusions |
| Endpoint | cause-specific recurrence, as deposited |
| M1 | clinical-only Cox model (5 parameters) |
| M2 (primary) | M1 plus four pre-registered pathway scores (9 parameters) |
| Estimand | ΔC_corrected = C_corrected(M2) − C_corrected(M1), the optimism-corrected difference in Harrell's C |
| Target output | a 95 % percentile interval for ΔC_corrected, with its Monte Carlo error |

The cohort, endpoint, panel membership, scoring method, gene-level representation and model forms were all fixed before any outcome was examined. They are recorded in dated protocol documents that the pipeline pins by hash: `docs/PHASE9_MODELING_PROTOCOL_DRAFT.md` and `docs/PHASE9_DELTAC_INTERVAL_PLAN.md`.

## Why a nested bootstrap

Harrell's bootstrap optimism correction gives a point estimate of ΔC_corrected but no valid interval for it. To get an interval, the whole correction procedure has to be repeated on resampled data:

```
for each outer replicate k = 1..B_out          (B_out = 500; may grow to 750 ... 3000)
    D*_k = stratified bootstrap of the cohort
    fit M1 and M2 on D*_k                       -> apparent C on D*_k
    for each inner replicate b = 1..B_in        (B_in = 100)
        J*_k,b = bootstrap of D*_k
        fit M1 and M2 on J*_k,b, evaluate on D*_k -> optimism(M, k, b)
    ΔC_corrected(k) = ΔC_apparent(k) − mean_b[opt(M2) − opt(M1)]
interval = percentiles of {ΔC_corrected(k)}
```

M2's ridge penalty is re-tuned by five-fold cross-validation inside every fit, so a single run involves on the order of 10⁵ penalised Cox fits. The main design choices, each argued in the plan:

- **Pairing.** M1 and M2 are fitted on the same resamples. In the earlier runs this removed about 44 % of the inner Monte Carlo variance of the optimism difference (correlation ≈ 0.68).
- **Inner bootstrap size.** B_in = 100 inflates the outer standard deviation by about 1.2 %. That makes the interval slightly conservative; I report the inflation rather than trying to remove it.
- **Fit failures.** Inner replicates in which either model fails are dropped as a pair, and a bounding analysis checks how much this could matter.
- **Reproducibility.** Five independent full-size runs measure Monte Carlo reproducibility directly.
- **Endpoint Monte Carlo error.** It is estimated with kernel quantile methods at three bandwidths. The acceptance rule (EP-1) uses the largest of the three estimates; the reliability screens use the smallest.
- **Adaptive sizing.** B_out is increased along a fixed ladder until the endpoint Monte Carlo error is below 5 % of the interval width. A failed stability check (W5) also forces a top-up.

## Why there is a corrected pipeline

The original Phase 9 run (3–4 November 2026) produced a point estimate of **ΔC_corrected = +0.0043**, with no interval. A later audit of that run found a deviation from the plan. The cross-validation used to tune λ started from the same random-number-generator state as the bootstrap resample it was tuning on, instead of drawing from a separate stream as the plan specifies. I confirmed this empirically in every case examined (9 of 9). The deviation is documented in a protocol-deviation note that is kept with the historical outputs; those outputs are not included in this repository.

The historical outputs are kept unchanged and are pinned by hash; the new pipeline uses them only as reference data for verification. The corrected pipeline (`src/phase9_corrected/`) implements the specified stream separation and calls the original fitting code (`src/phase9_execution.py`) without modification. Its outer bootstrap resamples are checked against the hashes stored by the historical run, so the outer resamples match the historical run exactly.

## Reproducibility and audit design

**Random-number streams**
- Every random draw comes from its own `SeedSequence` child, indexed by run, outer replicate k and inner replicate b. The streams are: outer resample, apparent-fit CV (seeded from a separate entropy word), inner resample, inner CV and an inner reserve. No stream is used twice.
- Before any replicate is computed, a gate enumerates every initial generator state the run could use and checks they are all distinct. For the configured design that is about 4.5 million states.
- Each replicate record stores fingerprints of the stream states and resample multiplicities it actually consumed, so stream assignment can be checked afterwards.

**Long-running execution**
- Checkpoints are append-only. Each run's checkpoint carries a header that binds it to the configuration and input manifest.
- Each invocation writes durable start, end and authorisation records. Before anything new is dispatched, the run history is reconstructed from these records and checked against the checkpoint bytes.
- Records are created exclusively and synced with `F_FULLFSYNC`. An invocation that died without an end record is classed as unterminated, and it cannot be resumed until it is explicitly acknowledged.
- A runtime monitor projects completion time after a warm-up period. It halts the run if the projection exceeds twice the declared budget; resuming then requires a single-use authorisation reference.
- The run stops if more than 5 % of completed outer replicates are non-estimable (strict literal rule, `20·b > n`), checked on every prefix of the record stream.
- Each run is classified as ineligible, unverifiable, incomplete or eligible, with every reason kept. The same assessment controls dispatch, sizing and reporting.
- Increasing B_out (a top-up) is allowed only on the basis of a sizing decision that the frozen code recomputes from all five checkpoints.

**Code identity**
- A freeze record lists SHA-256 hashes of the 27 source files the analysis depends on. It also records an environment fingerprint: interpreter path and hash, installed package versions, the NumPy BLAS/LAPACK build, thread settings and interpreter launch flags.
- The run configuration stores the freeze record's own hash.
- A shell launcher (`scripts/42_phase9_corrected_launch.sh`) checks that hash and every listed file, using system `jq` and `shasum`, before Python starts. It then starts Python in isolated mode with a fresh, empty bytecode directory, so no stale `.pyc` file can be loaded.
- Once running, the pipeline re-checks every repository module it has loaded against the freeze.
- Every invocation records the freshly executed pre-execution gate result it relied on.

**Post-run verification** (`scripts/43_phase9_corrected_verify.py`)

The runner works through these checks in order:
1. Every invocation is bound to the same configuration, manifest and freeze.
2. Interrupted stretches of work are recomputed and compared byte-for-byte with the checkpoint.
3. Run eligibility is re-assessed.
4. The sizing decision is recomputed at every ladder level.
5. The record-level checks:
   - V1 regenerates outer and inner resamples and compares them with stored hashes and with the historical code path;
   - V2 checks stream assignment;
   - V3 checks that the clinical-only model is unaffected by the seed correction;
   - V4 checks the as-run variant against the historical results;
   - V5 recomputes 20 replicates in full.

Each check reports PASS, FAIL, NOT_APPLICABLE or NOT_COMPARABLE, with counts and reasons. Missing evidence never counts as a pass, and results may be reported only if the overall verdict is PASS.

## Software testing

The corrected pipeline has 549 automated tests; the full local suite, including earlier phases of the project, has 1,075, all passing. Most of the corrected-pipeline tests use synthetic data and stand-in model fits. A few check hashes or loading of the locally held study inputs and historical outputs. Those files are not distributed here, so those tests will not pass on a fresh clone.

Besides ordinary unit and integration tests, the suite includes:
- **independent reference calculations**, for example a separate pure-Python implementation of the sizing and Monte Carlo error rules, compared to a relative tolerance of 10⁻¹⁰;
- **checks against the real enumerator**, for the gate's state count;
- **tamper tests**:
  - a failed gate result edited to read PASS;
  - altered or substituted evidence files;
  - stale bytecode;
  - a module placed at the repository root to shadow a library;
  - every refusal path in the launcher.

These tests show that the software does what the protocol says. They say nothing about colon-cancer biology or about the size of ΔC.

## Current status

| Component | Status |
|---|---|
| Protocol, estimand, cohort, models | fixed in earlier phases |
| Historical point estimate, ΔC_corrected = +0.0043 | historical only: no interval, and affected by the seed deviation |
| Corrected pipeline, launcher, verification runner, freeze tooling | implemented, tested, committed |
| Candidate freeze record (27 files, SHA-256 `618f20eb…a790`) | generated and checked against the committed code; approval pending |
| Analysis decisions 1–11 in the run configuration | entered, with dates and references; the configuration passes static validation except for the execution authorisation |
| Decision 12 (recording λ when a fit fails) | implemented and tested; its entry in the configuration's decision log is unresolved, because the template does not list the chosen option (the pipeline does not read this entry) |
| Freeze approval (Decision 13) | pending |
| Authorisation to run (Decision 14) | pending |
| Corrected bootstrap (5 runs × 500 outer × 100 inner replicates) | not run |
| Corrected interval for ΔC_corrected, and any conclusion | pending: none exists yet |

## Remaining steps

1. Re-hash and review the candidate freeze, install it as `config/analysis/PHASE9_CORRECTED_FREEZE_RECORD.json`, and record its hash and the driver hash in the configuration.
2. Before running, re-hash the freeze record and the launcher, and check that the repository root contains no importable Python files or packages other than `src` and `tests`.
3. Record the execution authorisation in the configuration.
4. Run `scripts/42_phase9_corrected_launch.sh run <config> <run> 500` for runs 0–4, with the BLAS/OpenMP thread variables unset as they were at freeze capture. Then run `scripts/42_phase9_corrected_launch.sh verify <config>`.

## Repository layout

| Path | Contents |
|---|---|
| `src/phase9_corrected/` | replicate engine, random streams, resampling, fitting adapter, checkpointing, orchestration, statistics (EP-1, W5, sizing), reliability checks, gate, verification and post-run runner |
| `src/phase9_execution.py`, `src/pathway_scoring.py`, `src/score_standardisation.py` | original fitting and scoring code, used unchanged |
| `scripts/40_phase9_dc_interval.py` | original Phase 9 driver, kept as the reference code path |
| `scripts/41`–`44` | corrected driver, launcher, verification entry point, freeze preparation |
| `config/analysis/PHASE9_CORRECTED_RUN_CONFIG_OWNER_BLANK.json` | run configuration with the recorded analysis decisions |
| `docs/PHASE9_DECISION13_FREEZE_PROCEDURE.md` | freeze, hash-anchoring and launch procedure |
| `tests/test_phase9_corrected_*.py` | tests for the corrected pipeline |

```bash
.venv/bin/python -m pytest -q tests/test_phase9_corrected_*.py
```

## Methods and tools

Cox proportional-hazards models with ridge penalties, Harrell's C and bootstrap optimism correction, nested bootstrap inference with Monte Carlo error estimation, NumPy `SeedSequence` stream design, hash-based code and data provenance, crash-consistent checkpointing, and test design with independent reference implementations. Python (NumPy, pandas, lifelines, pyarrow, pytest) and Bash.
