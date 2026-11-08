# Phase 9 corrected campaign — independent checkpoint-level audit

**Date:** 8 November 2026. **Classification: INDEPENDENTLY REPRODUCED.**
**Scope:** read-only. No checkpoint, configuration, frozen code or official result was changed.
**Tool:** `scripts/48_phase9_independent_audit.py` (SHA-256 `1227fce9…eccbd`). Python standard library only: no import of `src.phase9_corrected`, `scripts/46` or NumPy.
**Output:** `independent_audit_results.json` in this directory (SHA-256 `b67560a8…e58d`), created exclusively.

## 1. Independence of the calculation path

| Quantity | Official path (`scripts/46`) | Independent path (`scripts/48`) |
|---|---|---|
| ΔC_corrected(k) | stored aggregate `dC_corrected` | rebuilt from the 100 per-inner-replicate `C_fit` − `C_eval` values and the apparent fits, with pairing re-derived from each fit's `ok` flag |
| Quantiles | `numpy.quantile(method="linear")` | hand-written Hyndman–Fan type 7 |
| Normal CDF and inverse | `scipy.stats.norm` | `statistics.NormalDist` |
| Kernel SE and bandwidth | frozen `statistics.py` with NumPy | pure-Python sums; Silverman (1986) eq. 3.31; plan v2.1 section 2.2 |
| Checkpoint integrity | frozen `checkpoint.py` | raw bytes, `hashlib`; canonical JSON re-serialised with the standard library |
| θ̂ | `theta_hat_from_execution` | recomputed from `monte_carlo` in the execution record, and separately as ΔC_app(D) − ΔOpt(D) |

## 2. Completeness and integrity

| Check | Result |
|---|---|
| Runs present | 5 (0–4) |
| Outer replicates | 750 per run; **3,750 in total**; k = 0–749 in order, each exactly once |
| Inner bootstrap pairs | 100 attempted and 100 paired-complete for every outer replicate; **375,000 in total**; 0 fit failures |
| Record status | all `ok`, no failure entries |
| Per-record SHA-256 (canonical JSON) | 3,750 of 3,750 valid |
| Checkpoint headers | all bound to configuration `bf939a4e…6a56` and manifest `1ff24387…3f27` |
| Checkpoint SHA-256 = final end record | 5 of 5 |
| Original 500-record prefix | byte-identical to the 500-level end record in 5 of 5 runs |
| Top-up provenance | every second invocation started from the 500-level checkpoint, with path `topup`, n_start 500 and n_new 250 |
| Rebuilt ΔC_corrected, ΔC_app and ΔOpt vs stored | maximum absolute difference 2.1×10⁻¹⁷ |

| Run | Bytes | Checkpoint SHA-256 | 500-record prefix SHA-256 |
|---|---:|---|---|
| 0 | 112,951,341 | `0e27494ff875dcd7512355d213ee22f6bda2ce9a58595790b3a8b3b56b5dab15` | `e55a71e06d0ba603…` |
| 1 | 112,957,953 | `9dfe640b32793b5e21d40fa4723c65f8455b15282f0ad443f7a26f645a9b7f43` | `ef360ca0ab6f9aa4…` |
| 2 | 112,960,361 | `0204da98f37ce4af5b5c7cdc61c114c602178f190bb25ca82dcba47b360e91a3` | `9dbf6359e08569a7…` |
| 3 | 112,949,339 | `05274f0e5309454f9340e8a427304cd0d94a1c3538059784d85ba539fffc1e70` | `9c47c29e38c1cdeb…` |
| 4 | 112,959,127 | `af7530077f5ce68ea874034d796a129327056afe20951451699d30559e58e326` | `5d007755ec6f5814…` |

## 3. Reconstructed results

θ̂ = **0.004316342804030704**, matching the official value exactly. The identity θ̂ = ΔC_app(D) − ΔOpt(D) = 0.013963503027202173 − 0.009647160223171487 holds to 1.7×10⁻¹⁷.

| Run | Percentile 95 % (primary) | Width | Basic | BC (z₀ only) | z₀ (# < θ̂) | A7 | EP-1 ratio lower / upper |
|---|---|---:|---|---|---|---:|---|
| 0 | [−0.008618, +0.051887] | 0.06050 | [−0.04325, +0.01725] | [−0.01118, +0.02993] | −0.475 (238) | 0.572 | 0.0073 / 0.0308 |
| 1 | [−0.008962, +0.046893] | 0.05585 | [−0.03826, +0.01759] | [−0.01163, +0.02486] | −0.548 (219) | 0.525 | 0.0092 / 0.0424 |
| 2 | [−0.009292, +0.048372] | 0.05766 | [−0.03974, +0.01792] | [−0.01182, +0.02687] | −0.544 (220) | 0.528 | 0.0089 / 0.0315 |
| 3 | [−0.009150, +0.051784] | 0.06093 | [−0.04315, +0.01778] | [−0.01171, +0.03186] | −0.431 (250) | 0.558 | 0.0073 / 0.0379 |
| 4 | [−0.009078, +0.049436] | 0.05851 | [−0.04080, +0.01771] | [−0.01174, +0.02649] | −0.591 (208) | 0.542 | 0.0093 / 0.0257 |

- **Between-run SD (ddof = 1):** lower endpoint 0.000255, upper endpoint 0.002170.
- **EP-1:** ρ_max = 0.04245 ≤ 0.05, so EP-1 is met.
- **Pooled draws (n = 3,750):** percentile range [−0.009062, +0.050232]; mean +0.013936; median +0.012125; 20.24 % of draws below zero.

**Pooling rule.** Plan v2.1 designates **no** pooled interval. The primary result is the per-run percentile interval (sections 2.3, 7 and 9.1 step 5). The pooled range is reproduced here only because the official output reports it as descriptive. It is not a primary or separately established interval.

## 4. Reconciliation

- **Official results file** (`corrected_interval_results.json`, `1c504b13…eb60`). 138 quantities were compared:
  - per run: percentile and basic endpoints, BC endpoints, width, z₀ and its count, A7, A8, Silverman bandwidth, all 6 kernel SEs, distribution mean, median, SD and share below zero, mean ΔC_app, mean ΔOpt, and mean corrected C for each model;
  - between-run SD and range;
  - pooled endpoints;
  - ρ_max;
  - θ̂.

  The **maximum absolute difference is 2.2×10⁻¹⁶** (run 1 A7), which is floating-point summation order.
- **Final report** (`docs/PHASE9_CORRECTED_FINAL_REPORT.md`). Every interval endpoint, width, SD, z₀, A7 and EP-1 value quoted there matches the reconstruction at its printed precision.
- **Frozen verification** (`L00750.report.json`, `9558aaa5…ff25`): PASS, 60/0/0/0, reporting permitted. The independent counts agree with its V1 counts: 750 outer and 75,000 inner per run.
- **Plan v2.1.** The estimator rebuilt here is the section 3.1 procedure: paired exclusion, optimism evaluated on D*ₖ, ΔC_corrected = ΔC_app − ΔOpt. The interval is the type-7 percentile interval at α = 0.05.
- **θ̂ alignment audit** (`docs/PHASE9_THETA_ALIGNMENT_AUDIT.md`). The θ̂ components used there are reproduced exactly.

## 5. Limitations of this audit

- It does not refit any model. It verifies that the reported statistics follow exactly from the recorded per-fit C values, and that those records are complete, untampered and correctly linked across invocations. That the recorded C values are themselves correct model outputs rests on the frozen verification's V5: a bit-identical full recomputation of 20 replicates in Run 0 (k = 0, 25, …, 475).
- It shares the record format, but not the calculation code, with the official path.

**Outcome: INDEPENDENTLY REPRODUCED.** No discrepancy was found, and no official result was changed.
