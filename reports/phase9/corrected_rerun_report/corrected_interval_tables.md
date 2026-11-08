# Corrected nested-bootstrap interval tables (B_out = 750, B_in = 100)

Generated 2026-11-08T05:03:27.285952+00:00 by `scripts/46_phase9_corrected_report.py`. Conditional on the observed 136 recurrence events and 316 censored observations in the 452-patient A-19 population.

## T1. Per-run 95 % percentile intervals (primary) and diagnostics

| Run | n usable | Percentile (primary) | Width | Basic (diagnostic) | BC z0 (diagnostic) | z0 | A7 | A8 |
|---|---:|---|---:|---|---|---:|---:|---:|
| 0 | 750 | [-0.00862, +0.05189] | 0.06050 | [-0.04325, +0.01725] | [-0.01118, +0.02993] | -0.475 | 0.572 | 0.475 |
| 1 | 750 | [-0.00896, +0.04689] | 0.05585 | [-0.03826, +0.01759] | [-0.01163, +0.02486] | -0.548 | 0.525 | 0.548 |
| 2 | 750 | [-0.00929, +0.04837] | 0.05766 | [-0.03974, +0.01792] | [-0.01182, +0.02687] | -0.544 | 0.528 | 0.544 |
| 3 | 750 | [-0.00915, +0.05178] | 0.06093 | [-0.04315, +0.01778] | [-0.01171, +0.03186] | -0.431 | 0.558 | 0.431 |
| 4 | 750 | [-0.00908, +0.04944] | 0.05851 | [-0.04080, +0.01771] | [-0.01174, +0.02649] | -0.591 | 0.542 | 0.591 |

## T2. Between-run endpoint variability (SD ddof = 1, df = 4)

| Endpoint | Five values | SD | Range |
|---|---|---:|---:|
| lower | -0.00862, -0.00896, -0.00929, -0.00915, -0.00908 | 0.00025 | 0.00067 |
| upper | +0.05189, +0.04689, +0.04837, +0.05178, +0.04944 | 0.00217 | 0.00499 |

## T3. Kernel endpoint MC-SE (bandwidth multipliers 0.5, 1, 2) and EP-1 ratio

| Run | Endpoint | SE ×0.5 | SE ×1 | SE ×2 | EP-1 ratio (Option B) |
|---|---|---:|---:|---:|---:|
| 0 | lower | 0.00032 | 0.00038 | 0.00044 | 0.0073 |
| 0 | upper | 0.00168 | 0.00186 | 0.00182 | 0.0308 |
| 1 | lower | 0.00041 | 0.00048 | 0.00051 | 0.0092 |
| 1 | upper | 0.00219 | 0.00237 | 0.00188 | 0.0424 |
| 2 | lower | 0.00042 | 0.00047 | 0.00051 | 0.0089 |
| 2 | upper | 0.00155 | 0.00182 | 0.00177 | 0.0315 |
| 3 | lower | 0.00033 | 0.00039 | 0.00044 | 0.0073 |
| 3 | upper | 0.00206 | 0.00231 | 0.00205 | 0.0379 |
| 4 | lower | 0.00047 | 0.00051 | 0.00054 | 0.0093 |
| 4 | upper | 0.00117 | 0.00144 | 0.00150 | 0.0257 |

EP-1 rho_max = 0.0424; sizing action at this level: `stop_ep1_pass`.

## T4. Pooled draws (descriptive only, not a primary interval)

n = 3750; percentile [-0.00906, +0.05023]; mean +0.01394, median +0.01213, share below zero 0.202.

