# Statistical conventions adopted where plan v2.1 is silent

Every convention below was fixed **before** any result was examined. Each is either named in project precedent or is the textbook default. None of them changes a definition the plan does state. Where a convention could plausibly matter, the alternative is computed alongside as a cross-check, and the decision basis is stated.

| # | Item | Plan text | Convention adopted | Basis | Alternative reported? |
|---|---|---|---|---|---|
| 1 | Confidence level | not stated explicitly | **95 %**, two-sided (α = 0.05; quantiles 0.025 / 0.975) | Phase 9 primary reports "bootstrap 95 % pct"; the predecessor audit uses the 2.5–97.5 percentile spread; plan §2.2 arithmetic (expected SE ≈ 3.0 % of width with proxy SD 0.0061) implies width ≈ 2 × 1.96 SD | no |
| 2 | Quantile estimator | "percentile" | Hyndman–Fan **type 7** (numpy `linear`, the default) | plan §3.4 uses type 7 for knots; Phase 9 primary used numpy's default | no |
| 3 | θ̂ for basic and z₀ | θ̂ = +0.0043 (plan header, §3.2) | `monte_carlo.M2_primary.corrected_C_mean − monte_carlo.M1.corrected_C_mean` from the frozen `phase9_execution_results.json` = 0.004316342804030704 | the plan's +0.0043 is this difference (Phase 9 execution report) | ±2 MC-SE sensitivity of z₀ (report §18) |
| 4 | z₀ ties | `#{ΔC_corr(k) < θ̂}` | strict "<" as written; ties counted and reported, never split | plan text | n/a (0 ties) |
| 5 | z₀ boundary | not stated | if the proportion is 0 or 1: z₀ = ±∞, flagged, and no BC interval | standard | n/a (never occurred) |
| 6 | BC interval | §7 lists z₀ as a diagnostic; the audit (line 98) names a "bias-corrected (z₀-only) variant" | percentiles at Φ(2z₀ + z_{α/2}) and Φ(2z₀ + z_{1−α/2}), with **a = 0** (BCa excluded by decision 4) | Efron (1987) BC interval | — |
| 7 | Bandwidth "Silverman's rule" | named, variant unstated | **Silverman (1986) eq. 3.31**: h = 0.9 · min(s, IQR/1.34) · n^(−1/5), s with ddof = 1, IQR by type-7 quantiles | the rule of thumb Silverman recommends | yes: scipy's `silverman` (≈ 1.06 s n^(−1/5)) as a cross-check column (T2) |
| 8 | Kernel density | "Gaussian kernel density evaluated at the endpoint" | f̂(t) = (1/(nh)) Σ φ((t − xᵢ)/h), evaluated at the type-7 quantile estimate | plan text | — |
| 9 | Between-run SD | "SD of each endpoint across the five runs", df = 4 | **sample SD, ddof = 1** | the plan's df = 4 implies ddof = 1 | — |
| 10 | Which MC-SE governs EP-1 | not stated | **all prescribed estimators evaluated** (kernel × {0.5, 1, 2} and between-run SD); the outcome is reported per estimator; no single estimator is chosen | avoids choosing the estimator after seeing results | yes: the full 40-row grid |
| 11 | EP-1 interval | "(upper − lower)" | the **percentile (primary)** interval of each run | plan §7: percentile is primary | — |
| 12 | MC-1 grid for W5 | "final 40 % of B_out" | B ∈ {300, 400, 500} from the Phase 6B MC-1 grid {50, 100, 200, 300, 400, 500}, first B replicates **in k order** | Phase 6B MC-1; k is the canonical, order-independent index | fine trace (step 10) shown in F4 for context only |
| 13 | Bounding imputation (§6.4) | "missing paired optimism values imputed at the min and max observed within that outer replicate" | imputes the **paired difference** d_b = opt_M2,b − opt_M1,b at its within-replicate min / max | preserves pairing (§6.2) | n/a (vacuous: 0 missing) |
| 14 | W3 inner-noise inflation | "inner-noise SD inflation > 5 %" | noise var = mean_k var(d_b, ddof = 1) / n_b; SD inflation = √(Var_obs / (Var_obs − noise)) − 1 | plan §1.2 decomposition | variance inflation also reported |
| 15 | Pooled 2 500-draw interval | not designated | computed **descriptively only**, and labelled as such | — | — |
| 16 | Stored inner C precision | the driver rounds `inner_C` to 8 decimals | ΔOpt recomputed from `inner_C` agrees with the stored `dOpt` to ≤ 2.7 × 10⁻⁹. All intervals use the **stored** `dC_corrected`; W3 and bounding use `inner_C` | driver line 167 | — |
