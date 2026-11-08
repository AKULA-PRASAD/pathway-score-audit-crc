# DRAFT — NOT YET ADJUDICATED (version 2)

## Protocol deviation record: inner-CV seed derivation in the Phase 9 nested bootstrap

| | |
|---|---|
| Status | **DRAFT — NOT YET ADJUDICATED.** This is a review version. It is not a decision and changes no result. |
| Supersedes | nothing. Version 1 (`reports/phase9/h5_spotcheck/CV_SEED_DEVIATION_DRAFT.md`, sha256/16 `2866571cecf471af`) is **retained unchanged**. |
| What changed from v1 | Every source reference was re-verified line by line. Seed-state evidence was rerun as a saved, reproducible artifact (§4). Stored and inferred evidence are separated (§6). Stored M2 λ is described (§6.3). Phase 9 primary's sequential stream dependence is added (§8). Possible future investigations (§11) and disclosure considerations (§12) are added. |
| Prepared | 2026-11-04, under the master work order, Stage A |
| Analysis | Phase 9 nested bootstrap CI for the optimism-corrected ΔC, runs 0–4 (`reports/phase9/dc_interval/`) |
| Plan | `docs/PHASE9_DELTAC_INTERVAL_PLAN.md` v2.1, sha256/16 `d27c52a29930a338` |
| Driver | `scripts/40_phase9_dc_interval.py`, sha256/16 `4b1738a4d0ee0e39` (frozen, unchanged) |
| Shared fitting code | `src/phase9_execution.py`, sha256/16 `571d0c233f367be6` (frozen, unchanged) |
| Phase 9 primary script | `scripts/38_phase9_execute.py`, sha256/16 `82afbcf05f3b7514` |
| Historical results | **Not modified.** All five checkpoints keep their recorded hashes (re-verified 2026-11-04). |

---

### 1. Protocol requirement (plan v2.1 §8.2)

`docs/PHASE9_DELTAC_INTERVAL_PLAN.md`, re-verified 2026-11-04:

| line | text |
|---|---|
| 361 | `root_k = numpy.random.SeedSequence(entropy=[BASE, run_id])` |
| 362 | `outer  = root_k.spawn(B_out)                 # outer[k]  -> outer replicate k` |
| 363 | `inner  = outer[k].spawn(B_in)                # inner[b]  -> inner replicate b of outer k` |
| 364 | `cvseed = inner[b].spawn(2)                   # [0] = inner 5-fold CV, [1] = reserve` |
| 367–372 | "Each replicate's stream is a pure function of `(BASE, run_id, k, b)`", followed by four consequences: independence between runs, independence between layers, scheduling invariance, and resumability |

The plan requires the inner 5-fold CV that tunes λ to draw from a **dedicated child stream**, `inner[b].spawn(2)[0]`, separate from the stream that produces the inner bootstrap resample. The plan has no explicit line for the apparent-level CV seed, the one used when fitting on `D*ₖ` itself (see §10, Q1).

### 2. Actual implementation (frozen driver)

`scripts/40_phase9_dc_interval.py`, re-verified 2026-11-04:

| line(s) | code | role |
|---|---|---|
| 3 | "Implements docs/PHASE9_DELTAC_INTERVAL_PLAN.md v2.1 exactly." | docstring claim |
| 96–98 | `ss_root = SeedSequence(entropy=[BASE, run_id])`; `ss_outer = ss_root.spawn(k + 1)[k]`; `rng_o = default_rng(ss_outer)` | outer stream |
| 100 | `I = strat_boot(e, rng_o)` | outer resample, drawn from `ss_outer` |
| 116–117 | `pipeline_fit(..., LAM_GRID, np.random.default_rng(ss_outer))` | **apparent-fit CV rng: fresh generator on `ss_outer`** |
| 130 | `ss_inner = ss_outer.spawn(B_in)` | inner streams |
| 134–135 | `rng_i = default_rng(ss_inner[b])`; `J = strat_boot(eO, rng_i)` | inner resample, drawn from `ss_inner[b]` |
| 141–142 | `pipeline_fit(..., LAM_GRID, np.random.default_rng(ss_inner[b]))` | **inner-fit CV rng: fresh generator on `ss_inner[b]`** |

`src/phase9_execution.py`, re-verified 2026-11-04:

| line(s) | role |
|---|---|
| 131–138 | `_event_stratified_folds(e, k, rng)`: `rng.permutation` of event indices, then of censored indices, dealt into 5 folds |
| 141–144 | `tune_lambda(...)` builds the folds from the `rng` it is given |
| 196–201 | `pipeline_fit` calls `tune_lambda(..., rng)` **only if `p > n_clin`**; otherwise it sets `lam = 0.0` |

`rng` is not used anywhere else in `pipeline_fit`. No `spawn(2)` call exists anywhere in the driver.

### 3. Random-seed derivation: planned vs executed

| stream | planned (§8.2) | executed |
|---|---|---|
| outer resample `I*ₖ` | `root.spawn(B_out)[k]` | `root.spawn(k+1)[k]`, **equivalent** (§4) |
| apparent-fit CV (on `D*ₖ`) | not stated explicitly | fresh `default_rng(ss_outer)`, **the same starting state as the outer resample** |
| inner resample `J*ₖ,b` | `outer[k].spawn(B_in)[b]` | `ss_outer.spawn(B_in)[b]`, same |
| inner-fit CV | `inner[b].spawn(2)[0]`, a **separate child stream** | fresh `default_rng(ss_inner[b])`, **the same starting state as the inner resample** |

Executed behaviour: in every fit, the CV fold permutation consumes **the same initial random bits** that produced that fit's bootstrap resample. The algorithms reading those bits differ (`Generator.choice(replace=True)` for the resample, `Generator.permutation` for the folds), so the fold assignment is not a copy of the resample. But it is a deterministic function of the same random input. It is not the independent draw that §8.2 specifies.

### 4. Numerical evidence (stored, reproducible)

Script: `reports/phase9/final_project/validation/cv_seed_state_evidence.py`. Output: `validation/cv_seed_state_evidence.json`. It covers 9 cases from run 0, for k ∈ {0, 25, 475} at the apparent level and at inner b ∈ {0, 99}:

| check | result |
|---|---|
| `spawn(k+1)[k]` vs `spawn(500)[k]`: same spawn key and `generate_state(8)` | **True** in 9/9 |
| CV generator as executed: same PCG64 state as the resampling generator | **True** in 9/9 |
| CV generator as planned (`spawn(2)[0]`): same state as the resampling generator | **False** in 9/9 |
| first raw 64-bit outputs, e.g. k = 0, inner b = 0 | resample `12275464072138857609, 7070011973245522162 …`; CV executed **identical**; CV planned `5070398370141261833, 4941475203673471867 …` |

This repeats and extends the console checks reported on 2026-11-04 (k = 0, 25, 475; outer equivalence; state equality), which were not saved to a file.

### 5. Affected components

The **M2 λ selection** by inner 5-fold CV, at both levels:
- 5 × 500 = 2,500 apparent M2 fits;
- 5 × 500 × 100 = 250,000 inner M2 fits.

Through λ, the deviation reaches the fitted M2 models, M2's apparent C on `D*ₖ`, M2's `C_boot` and `C_outer`, `opt_M2_mean`, `dC_app`, `dOpt` and **`ΔC_corr(k)`**. It therefore reaches every interval, MC-SE, warning and EP-1 result computed from them.

**The size of any difference under plan-conformant seeding has not been measured.**

### 6. Unaffected components and the evidence for each

#### 6.1 Stored evidence (directly observed in the checkpoints)

- `lambda_apparent.M1 = 0.0` in **all 2,500** outer records (`results/s1_interval_results.json`, `distribution.lambda_apparent_M1_distinct`).
- The outer-resample multiplicities (`multiplicity_sha`) are reproduced exactly by the frozen driver for the 20 H5 replicates (`reports/phase9/h5_spotcheck/h5_comparison.json`).
- Event/censoring composition is 136/316 in every outer record.

#### 6.2 Inferred from the code path (not directly measured)

- **Inner-level M1 λ = 0.** The outer records don't store inner-level λ. The claim rests on `src/phase9_execution.py:196–201`: with no molecular block `p == n_clin`, so `tune_lambda` is never called. The stored apparent-level values in §6.1 agree with that code path but don't measure the inner level.
- **The resamples are drawn as planned.** `I*ₖ` and `J*ₖ,b` are drawn first, from generators seeded as §8.2 intends (lines 98–100 and 134–135). Creating a fresh CV generator afterwards cannot change draws that have already been made.
- **Independence between runs and between replicates holds.** The spawn hierarchy over `(BASE, run_id, k, b)` is intact. The overlap is *within* a single fit, between its resample and its CV.

#### 6.3 Stored M2 λ (descriptive only)

`lambda_apparent.M2` is stored for every outer record. In each run its median is ≈ 0.152 and its maximum is 10, the top of the λ grid (`results/s1_interval_results.json`, `distribution.lambda_apparent_M2`). These are the values **as executed**. On their own they say nothing about what λ would have been under the planned seeding.

### 7. Potential methodological implication (stated, not evaluated)

The CV fold assignment used to select M2's λ is **not seeded independently of the bootstrap resample it is applied to**.
- Whether this changes λ selection, the M2 optimism, `ΔC_corr`, or the interval is **unknown**.
- Whether any effect would be systematic or only add noise is **unknown**.
- This record does not call the deviation harmless or consequential.
- The large observed shift between the outer `ΔC_corr(k)` distribution and θ̂ (z₀ ≈ −0.43 to −0.58; final report §22) is **not attributed** to this deviation. Nothing links them.

### 8. Historical implementation comparison: Phase 9 primary

`scripts/38_phase9_execute.py`, re-verified 2026-11-04:

| line(s) | behaviour |
|---|---|
| 90–92 | apparent fits: a fresh `default_rng(SEEDS[0])` per model, passed to `pipeline_fit`. No resample at this stage |
| 117–118 | bootstrap: **one** generator per seed, `rng = default_rng(seed)` |
| 122 | `bi = boot_indices(e, rng)`: the resample comes from that generator |
| 126–129 | the **same, advancing** generator is passed to `pipeline_fit` for each model in turn |

**Correction carried forward from v1:** a verbal summary given before v1 said Phase 9 primary "used the same convention". That was imprecise.
- **Shared with the dc_interval driver:** no dedicated CV child stream; the CV uses the replicate's resampling generator.
- **Different:** Phase 9 primary never re-creates the generator, so its CV draws come from the bits **after** the resample. It does not replay the resample's starting state.
- **A further property, observed from the code:** since one generator advances through resampling, then M2's CV, then the next replicate's resampling, replicate *b+1*'s resample depends on how many draws replicate *b*'s CV consumed. Phase 9 primary replicates are therefore **not** a pure function of `(seed, b)`.

Whether Phase 9 primary conforms to its own protocol text has **not** been examined.

θ̂ = +0.0043, which the dc_interval basic interval and z₀ use, comes from Phase 9 primary. It was therefore produced under that **different** seeding convention.

### 9. Reproducibility versus conformance

- **Reproducibility.** H5 (plan §8.4) recomputed run-0 k = 0, 25, …, 475 with the frozen driver: 20/20 bit-identical `ΔC_corr`, 20/20 identical full records. That establishes that the **as-executed** procedure is deterministic and resumable.
- **Conformance.** For the CV seed, the as-executed procedure is not the one the plan specifies: plan line 364 versus driver lines 117 and 142. The driver's line-3 claim ("implements v2.1 exactly") is not true for this element.
- **H5 does not establish conformance.** Conformance doesn't establish correctness either: a conformant implementation could still be wrong in other respects.
- **The resume record froze the code as it was.** `RESUME_AUTHORIZATION.json` ("no change to … seeds", driver hash `4b1738a4d0ee0e39`) locked the code as run. It does not bear on whether that code matched §8.2.

### 10. Unresolved questions

1. **Which text is binding for the CV seed:** the plan's §8.2 line 364, or the frozen, hash-locked driver approved at resume? The plan also leaves the apparent-level CV seed unstated.
2. **Effect size:** does seeding the folds from the same bits as the resample change λ selection, the M2 optimism or `ΔC_corr`, and by how much compared with the Monte-Carlo error?
3. **Direction:** if there is an effect, does it bias (moving the interval's centre) or only add noise (changing its width)?
4. **Interaction with θ̂:** θ̂ comes from Phase 9 primary, which used a third convention (§8). Does that affect whether θ̂ and the outer `ΔC_corr(k)` are comparable in the basic interval and in z₀?
5. **Scope:** should Phase 9 primary and S-7 be checked against their own protocol's seeding text?
6. **Process:** should future pre-execution audits add a check that the seed hierarchy matches the plan?

### 11. Possible future investigations (for adjudication; none is performed or recommended here)

- A **separately labelled** sensitivity recomputation of a fixed subset of outer replicates (for example the H5 indices) under plan-conformant `spawn(2)[0]` seeding, compared replicate by replicate with the stored values. It would need explicit authorisation, would be reported beside the historical results and would never replace them. **It is not authorised by the current work order.**
- A code-level conformance audit of the Phase 9 primary and S-7 seeding against their protocols.
- A plan-level clarification of the apparent-level CV seed.

### 12. Disclosure considerations

Whatever the adjudication, the following would need to accompany any reported interval:
- **(a)** the deviation itself: plan line 364 against driver lines 117 and 142;
- **(b)** that every `ΔC_corr`-derived result depends on M2 λ tuning carried out under the as-executed seeding;
- **(c)** that H5 establishes reproducibility, not conformance;
- **(d)** that θ̂ comes from a third seeding convention (§8);
- **(e)** the adjudication outcome, once one exists.

Until then, the final-project outputs label every M2-dependent result: *"PROVISIONAL — calculated using the as-run implementation; CV-seed protocol deviation not yet adjudicated."*

---

*This record changes nothing: no plan, driver, frozen artifact, checkpoint or existing audit record has been modified, and version 1 is retained. It stays a draft until adjudicated.*
