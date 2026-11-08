# DRAFT — NOT YET ADJUDICATED

## Protocol deviation record: inner-CV seed derivation in the Phase 9 nested bootstrap

| | |
|---|---|
| Status | **DRAFT — NOT YET ADJUDICATED.** Prepared for review on request (2026-11-04). It is not a decision and does not change any result. |
| Analysis | Phase 9 nested bootstrap CI for the optimism-corrected ΔC (`reports/phase9/dc_interval/`, runs 0–4, all complete) |
| Plan | `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` v2.1 (sha256/16 `d27c52a29930a338`) |
| Driver | `scripts/40_phase9_dc_interval.py` (sha256/16 `4b1738a4d0ee0e39`, frozen; unchanged) |
| Shared fitting code | `src/phase9_execution.py` (sha256/16 `571d0c233f367be6`, frozen; unchanged) |
| Found during | Pre-execution review for the §8.4 / H5 reproducibility spot-check, 2026-11-04 |
| Historical results | **Not modified.** All five checkpoints keep their recorded hashes. |

---

### 1. File references

**Plan, §8.2 "Deterministic stream hierarchy"** (`docs/PHASE9_DELTAC_INTERVAL_PLAN.md`):

- line 361: `root_k = numpy.random.SeedSequence(entropy=[BASE, run_id])`
- line 362: `outer  = root_k.spawn(B_out)                 # outer[k]  -> outer replicate k`
- line 363: `inner  = outer[k].spawn(B_in)                # inner[b]  -> inner replicate b of outer k`
- line 364: `cvseed = inner[b].spawn(2)                   # [0] = inner 5-fold CV, [1] = reserve`
- lines 367–372: the stated consequences, including "Independence between layers — an inner stream descends from its own outer stream; no sharing."

**Driver** (`scripts/40_phase9_dc_interval.py`):

- line 3, docstring: "Implements docs/PHASE9_DELTAC_INTERVAL_PLAN.md v2.1 exactly."
- lines 96–98: `ss_root = SeedSequence(entropy=[BASE, run_id])`; `ss_outer = ss_root.spawn(k + 1)[k]`; `rng_o = default_rng(ss_outer)`
- line 100: outer resample `I = strat_boot(e, rng_o)`
- lines 116–117: apparent fit, `pipeline_fit(..., LAM_GRID, np.random.default_rng(ss_outer))`
- line 130: `ss_inner = ss_outer.spawn(B_in)`
- lines 134–135: `rng_i = default_rng(ss_inner[b])`; inner resample `J = strat_boot(eO, rng_i)`
- lines 141–142: inner fit, `pipeline_fit(..., LAM_GRID, np.random.default_rng(ss_inner[b]))`

**Shared fitting code** (`src/phase9_execution.py`):

- lines 131–138: `_event_stratified_folds(e, k, rng)` uses `rng.permutation` for the events, then for the censored cases
- lines 141–144: `tune_lambda(...)` builds the 5 folds from the `rng` it receives
- lines 196–197: `pipeline_fit` calls `tune_lambda(..., rng)` **only if** `p > n_clin`, i.e. only when a molecular (penalised) block is present
- `rng` is not used anywhere else in `pipeline_fit`

### 2. What the plan specifies

The inner 5-fold CV that tunes λ draws from its own child stream, `cvseed = inner[b].spawn(2)[0]`. It is separate from the stream that produces the inner bootstrap resample. The plan doesn't give the apparent-level CV seed on a separate line. By the same hierarchy it would be a child of `outer[k]`, but the plan does not state this explicitly (see §10, Q1).

### 3. What was actually executed in all five runs

- **Inner level** (100 × 500 × 5 fits per model): the inner resample `J` is drawn from `default_rng(ss_inner[b])` (line 134–135). The CV `rng` handed to `pipeline_fit` is a **new** `default_rng(ss_inner[b])` (line 142), built from the same SeedSequence. Its starting state is therefore identical to the starting state of the resampling generator. The CV fold permutations consume the same initial bit stream that the resample consumed.
- **Apparent level** (500 × 5 fits per model): the outer resample `I` comes from `default_rng(ss_outer)` (lines 98, 100). The apparent-fit CV `rng` is a new `default_rng(ss_outer)` (line 117), so the same relationship holds one level up.
- **No `spawn(2)` child stream is created anywhere in the driver.**

The outer-stream derivation `spawn(k + 1)[k]` (line 97) **is equivalent** to the plan's `spawn(B_out)[k]`; see §4. It is not part of this deviation.

### 4. Numerical verification already performed (2026-11-04; read-only, nothing written)

These checks were run with the project interpreter `.venv/bin/python -I`, using numpy from the frozen environment:

| Check | Result |
|---|---|
| `SeedSequence([BASE,0]).spawn(k+1)[k]` vs `.spawn(500)[k]`, k = 0, 25, 475: same spawn key and same `generate_state(8)` | **True** (equivalent) |
| `default_rng(ss).bit_generator.state` (CV, as executed) == `default_rng(ss).bit_generator.state` (resample) | **True** (identical starting state) |
| `default_rng(ss.spawn(2)[0])` state (CV, as planned) == resample starting state | **False** |

### 5. Affected component

**M2 λ tuning by inner 5-fold CV**, at both the apparent level and the inner level, in every outer replicate of every run. In total that is 5 × 500 apparent M2 fits plus 5 × 500 × 100 inner M2 fits.

Through λ, the deviation reaches the fitted M2 model, then C_boot and C_outer for M2, the M2 mean optimism, Δopt, and `ΔC_corr`. **Whether and how much these quantities would differ under plan-conformant seeding has not been measured** (see §10).

### 6. Component identified as unaffected so far

- **M1** has no molecular block, so `p == n_clin` and `pipeline_fit` never calls `tune_lambda`; it uses λ = 0 (`src/phase9_execution.py:196–201`). Supporting evidence: `lambda_apparent.M1` is `0.0` in **all 2,500** stored outer records. Inner-level λ is not stored, so the inner-level M1 claim rests on the code path, not on stored values.
- **The resamples themselves** (outer `I`, inner `J`, and thus `multiplicity_sha`) are drawn first from generators seeded as the plan intends. The deviation concerns only the CV generator.
- **Between-run and between-replicate independence:** the plan's spawn hierarchy for (BASE, run_id, k, b) is intact. The overlap is within a single fit, between that fit's resample and that fit's CV.

### 7. Potential methodological implication (stated, not evaluated)

The CV fold assignment used to select M2's λ is **not seeded independently of the bootstrap resample it is applied to**. Both come from the same initial bits of one PCG64 stream. The algorithms that consume those bits differ: `Generator.choice(..., replace=True)` for the resample, and `Generator.permutation` for the folds. So the fold assignment is not a copy of the resample. But it is a deterministic function of the same random input, rather than an independent draw as §8.2 intends.

**No evidence has been gathered about whether this changes λ selection, the M2 optimism, or `ΔC_corr`, or the interval derived from them.** This record does not describe the deviation as harmless or as consequential.

### 8. Phase 9 primary convention — verified, and **different**

An earlier verbal summary (2026-11-04, before this draft) said Phase 9 primary "used the same convention". **That statement was imprecise and is corrected here.** Source: `scripts/38_phase9_execute.py` (sha256/16 `82afbcf05f3b7514`):

- **Bootstrap** (lines 117–129): one generator per seed, `rng = default_rng(seed)`. `boot_indices(e, rng)` draws the resample, and then the **same, advancing** generator object is passed to `pipeline_fit` (lines 128–129). The CV draws therefore come from the bits **after** the resample draws, continuing sequentially across models and replicates.
- **Apparent fits** (lines 90–92): a fresh `default_rng(SEEDS[0])` for each model, with no resample at that stage.

**Shared with the dc_interval driver:** no dedicated CV child stream; the CV uses the replicate's resampling generator.

**Different:** Phase 9 primary does **not** re-create the generator, so its CV does not replay the resample's starting state. The dc_interval driver does.

Whether Phase 9 primary conforms to its **own** protocol is outside this record and has not been checked (see §10, Q5).

### 9. Deterministic reproducibility vs. conformance to the written protocol

These are separate properties:

- **Deterministic reproducibility** asks whether the executed procedure regenerates its stored values bit-identically. That is what §8.4 / H5 tests, and it was tested with the frozen driver (Option A, approved 2026-11-04). The result is recorded in `h5_comparison.json` in this directory. An H5 PASS establishes reproducibility of the procedure as executed. **It does not establish conformance to §8.2.**
- **Conformance** asks whether the executed procedure is the one the plan specifies. For the CV seed it is not: §8.2 line 364 versus driver lines 117 and 142. The driver's docstring (line 3) says it implements v2.1 "exactly"; for this element, it does not.

The seeding as executed was covered by the driver hash frozen in `RESUME_AUTHORIZATION.json` (`driver_sha256_16 4b1738a4d0ee0e39`, "no change to … seeds"). That record froze the code as it was. It does not bear on whether the code matched §8.2.

### 10. Unresolved questions for review

1. **Which text is binding for the CV seed:** the plan §8.2 line, or the frozen driver as approved and hash-locked? Relatedly, the plan does not state the apparent-level CV seed explicitly.
2. **Statistical impact:** does seeding the CV fold assignment from the same bits as the resample change λ selection, the M2 optimism, or `ΔC_corr`, and if so by how much compared with the Monte Carlo error? This has not been measured. One way to measure it, listed for adjudication only and **not** a recommendation to replace any historical result, would be a separately labelled sensitivity recomputation of a subset of replicates under plan-conformant `spawn(2)[0]` seeding, compared with the stored values.
3. **Direction of any effect:** if there is an effect, could it bias the M2 optimism systematically (affecting the interval's centre), or only add noise (affecting its width)? Unknown.
4. **Reporting:** how should the deviation be disclosed alongside any interval computed from the five runs, whatever the adjudication?
5. **Scope beyond dc_interval:** should the Phase 9 primary and S-7 runs (`scripts/38_*`, and the S-7 script) be checked against their own protocol's seeding text? Their convention differs from both §8.2 and the dc_interval driver (§8).
6. **Process:** the deviation was not caught before execution, despite the "implements v2.1 exactly" claim. Should a seed-hierarchy conformance check be added to future pre-execution audits?

---

*This record changes nothing: no plan, driver, frozen artifact, checkpoint or existing audit record (`HALT_RECORD.json`, `RESUME_AUTHORIZATION.json`) has been modified. It stays a draft until adjudicated.*
