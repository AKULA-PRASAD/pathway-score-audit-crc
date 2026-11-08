# Phase 9 — Estimand S-7: Representative-IQR Representation Sensitivity

**Version:** 1.0 · **Executed:** 2026-11-03T05:24:20Z · **Runtime:** 4520.3 s · **Exit code 0**
**Script:** `scripts/39_phase9_execute_s7.py` · **Engine:** `src/phase9_execution.py` (unchanged)
**Results:** `reports/phase9/execution_s7/phase9_s7_results.json`

> **SENSITIVITY-ONLY.** Median probe aggregation **remains the primary representation** (D8D-3,
> Phase 2B A-04), regardless of anything below. M3 was not fitted; the batch investigation was
> **not** repeated (S-7's frozen definition prescribes an aggregation-choice re-fit only).

---

## 1. Pre-execution verification

| Check | Result |
|---|---|
| Frozen S-7 definition | `representative_iqr`, `discovery_443`, 21 815 genes, 0 ties, `validation_samples_used: false` |
| **Selection population** | 443 samples, **443/443 discovery**, **0 validation, 0 non-tumoural**; equals *all* full-series discovery arrays |
| Frozen mappings | `frozen_gene_to_probe_discovery443.csv` `c386d5f1…` (21 815 rows, all `iqr`/`discovery_443`); matrix `e636e259…` |
| Gene index / samples | **identical** to the median matrix |
| A-19 population | **452 / 136 events / 316 censored** — identical to the primary run |
| Primary results | `05d13b8b339fa7d4`, untouched; S-7 wrote to a **separate** directory |
| Config | seeds ×5, B=500, patient-level, event-stratified, inner CV k=5, MC-1…MC-4, 5 % threshold |

**One check initially failed and was run down.** The freeze record's `selection_samples_sha256`
did not reproduce. Cause: `sha_list()` (`scripts/22_freeze_representation.py` L46) hashes
`"\n".join(items)` **without** a trailing newline; appending one yields the *file* hash. Recomputed
per the actual definition: **MATCH**. **Not a discrepancy** — the verification was wrong, not the
record.

## 2. Population and completion

452 patients, **136 events**, 316 censored — byte-identical to the primary run. B = 500 × 5 seeds,
`B_effective = 500` in every seed for every model (**7 500 fits**), exit 0, no traceback.

## 3. Results

| model | df | λ (S-7) | λ (median) | C apparent | C corrected | calib. slope corrected |
|---|---:|---|---|---|---|---|
| M1 | 6 | 0 | 0 | 0.6242 | **0.6016** | 0.8327 |
| M2 four-set | 10 | 0.1520 | 0.3511 | 0.6387 | **0.6059** | 0.8324 |
| M2 five-set | 11 | 0.1520 | 0.1520 | 0.6388 | **0.6041** | 0.8193 |

**M1 is bit-identical across representations** — max |Δβ| = **0.000e+00** over all six clinical
terms, identical C, logLik and knots. This is the built-in control: M1 contains no molecular term,
so the representation cannot reach it. **Control passes.**

### Molecular coefficients, M2 four-set (per 1 training-fold SD)

| pathway | β median | β S-7 | Δ | HR S-7 | 95 % CI (S-7) | sign |
|---|---|---|---|---|---|---|
| EMT | +0.0672 | +0.0804 | +0.0132 | 1.084 | 0.937–1.254 | same |
| E2F_TARGETS | −0.0991 | −0.1370 | −0.0379 | 0.872 | 0.758–1.003 | same |
| WNT_BETA_CATENIN | +0.0154 | +0.0105 | −0.0049 | 1.011 | 0.879–1.162 | same |
| P53_PATHWAY | −0.0047 | −0.0067 | −0.0019 | 0.993 | 0.852–1.158 | same |

**No sign flips. Every CI includes 1 — in both representations.**

## 4. Median vs representative-IQR

| model | C corr median | C corr S-7 | absolute Δ | relative | MC-SE | exceeds MC-SE? |
|---|---|---|---|---|---|---|
| M1 | 0.6016 | 0.6016 | **+0.0000** | 0.00 % | 0.00058 | no |
| **M2 four-set** | 0.6059 | 0.6059 | **−0.0000** | −0.00 % | 0.00047 | **no** |
| **M2 five-set** | 0.6066 | 0.6041 | **−0.0025** | −0.41 % | 0.00041 | **YES** |

Corrected calibration slopes: M1 0.8327 → 0.8327, four-set 0.8620 → 0.8324 (−0.0296), five-set
0.8364 → 0.8193 (−0.0171). All well below 1 under **both** representations — overfitting is a
property of the analysis, not of the representation.

### ΔC (M2 four-set − M1)

| | paired out-of-original ΔC | 95 % pct | corrected ΔC |
|---|---|---|---|
| median | +0.0121 | **[+0.0002, +0.0240]** | **+0.0043** |
| **S-7** | +0.0105 | **[−0.0001, +0.0224]** | **+0.0043** |

> **The one marginal result in the primary analysis does not replicate.** Under median the paired
> interval *barely* excluded zero (lower bound +0.0002); under representative-IQR it **includes
> zero** (−0.0001). The two bounds differ by 0.0003 — this is a knife-edge, not a reversal of
> evidence, and it shows the median-run exclusion was never robust.
>
> **The optimism-corrected ΔC is +0.0043 under both** — identical to four decimals. That is the
> estimate that matters, and it is unchanged by representation.

## 5. MC-1 … MC-4 (S-7)

**MC-1** final-40 % trace spread: M1 **0.00033**, four-set **0.00021**, five-set **0.00022** — stable.
**MC-2** five-seed range: 0.00157 / 0.00115 / 0.00058.
**MC-3** — `MC-SE ≤ 0.10 × sampling SE`, reported separately:

| model | MC-SE | sampling SE | 0.10 × SE | ratio | result |
|---|---|---|---|---|---|
| M1 | 0.00058 | 0.00944 | 0.00094 | 0.062 | **PASS** |
| M2 four-set | 0.00047 | 0.00990 | 0.00099 | 0.048 | **PASS** |
| M2 five-set | 0.00025 | 0.01030 | 0.00103 | 0.024 | **PASS** |

**All pass at B = 500** — no B-increase required; B was not lowered; seeds, thresholds and stopping
rules unchanged. **MC-4:** MC-SE accompanies every corrected estimate and is never summed with the
sampling SE. *Caveat:* sampling SE is a normal-approximation proxy, (p97.5 − p2.5)/3.92.

## 6. Non-estimable replicates

Attempted **7 500** (2 500 × 3 models); successful **7 500**; **non-estimable 0 (0.0 %)**; failure
reasons none; `redrawn: false`; 5 % threshold not approached. No estimate is flagged UNRELIABLE.

## 7. WNT-specific robustness

| pathway | score r (med vs IQR) | β median | β S-7 | relative change |
|---|---|---|---|---|
| EMT | 0.9919 | +0.0672 | +0.0804 | 19.6 % |
| **E2F_TARGETS** | **0.9975** | −0.0991 | −0.1370 | **38.2 %** |
| **WNT_BETA_CATENIN** | **0.8573** | +0.0154 | +0.0105 | 31.6 % |
| P53_PATHWAY | 0.9718 | −0.0047 | −0.0067 | 41.0 % |

> **Score-level correlation did not predict coefficient stability.** WNT — the weakest-agreeing
> score (0.857) — showed a *smaller* relative coefficient change than E2F, whose scores agree at
> **0.9975**. P53, agreeing at 0.972, moved most in relative terms.
>
> The explanation is arithmetic, not biological: these coefficients are all near zero, so a tiny
> absolute shift is a large *relative* one. **None of these movements is interpretable as
> stability or instability of a real effect, because none of the effects is distinguishable from
> zero in either representation.** The Phase 8D warning — *do not assume score correlation
> guarantees coefficient stability* — is borne out, though not in the direction anticipated.

## 8. P53–PI3K collinearity

| | P53↔PI3K r | condition | min eigenvalue | pairs \|r\| > 0.5 | VIF P53 | VIF PI3K |
|---|---|---|---|---|---|---|
| median | +0.6360 | 7.084 | 0.2762 | **1** | 1.96 | 2.19 |
| **rep-IQR** | **+0.6883** | **9.550** | 0.2262 | **2** | 2.26 | 2.67 |

**Collinearity is worse under representative-IQR**, and a **second** pair crosses 0.5
(E2F ↔ PI3K, r = +0.5448). Consistent with this, the five-set is the only model whose corrected C
moves beyond Monte-Carlo noise (−0.0025), and PI3K's coefficient falls from +0.0417 to +0.0103.
**This strengthens, not weakens, the case for PI3K remaining sensitivity-only** — and it was
decided outcome-blind, before any of these fits.

## 9. Integrity

**526 tests pass.** Bundle VERIFICATION PASSED; `27_freeze_m2_specification` and
`29_freeze_a19_missing_data` exit 0; Phase 9 audit **36/36**; raw checksums OK.
**Primary median results unchanged** (`05d13b8b339fa7d4`). **14/14 frozen artifacts unchanged**,
including both expression matrices and the frozen IQR mapping. Amendment chain `M2-A1 → M2-A2`
resolves at every revision. All 11 S-7 input hashes match disk. Python 3.13.13, numpy 2.5.4,
pandas 2.3.3, lifelines 0.30.3; seeds `[20261102…20261106]`, B = 500. **Zero commits.**

## 10. Interpretation and limitations

**Are the main findings robust to gene-representation choice?**

**For the primary four-set model: yes, to the precision available.** Corrected C is
**0.6059 under both** representations — agreement to four decimals, far inside the MC-SE. The
corrected ΔC is **+0.0043 under both**. No molecular coefficient changes sign, and every CI
includes 1 in both. M1 is identical by construction.

**But robustness here is agreement about a null-consistent result.** Showing that two
summarisations of the same matrix produce the same near-zero estimate is a statement about
pipeline stability, **not** evidence of a real effect — and it is **not independent replication**:
same patients, same genes, same deposited matrix, summarised differently.

**Three qualifications:**

1. **The marginal ΔC interval does not replicate.** Median excluded zero by +0.0002; S-7 includes zero by −0.0001. Any reading that leaned on that exclusion is unsupported.
2. **The five-set is representation-sensitive.** It is the only model whose corrected C moves beyond MC-SE, in the direction predicted by its worsened collinearity.
3. **Overfitting persists in both** — corrected calibration slopes 0.82–0.86, so predictions are too extreme by roughly 15 % regardless of representation.

**Not addressed here, by instruction:** the corrected-ΔC confidence interval and the
batch-perturbation coefficient instability. Both remain open and require separate authorisation.
No claim of clinical utility, biological causality, or external validity is made or supported.
