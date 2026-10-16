# Generated v0.1.1 pilot results

Regenerate with `Rscript scripts/run_all.R`. Source configuration: `R/config.R`.

## Run and interpretation

Pilot: 2000 designs; calibration: 300 designs per case; fixed-design oracle check: 20000 noise draws.

n=1500; p_G=0.12; alpha=(-0.3, 0.4, 0.9, 1); beta=(0, 0.7, 0.5, 0.3); gamma=0.10; tau0=0.25; Delta=0.25; sigma=1.00; theta=0.50.

Known-propensity trimming: [0.10, 0.90]. RNG: Mersenne-Twister / Inversion / Rejection. Calibration seed = 100000+r; pilot seed = 200000+r. The same calibration seeds are used in every case.

The outcome includes gamma*G*X, while the five candidate regressions are unchanged. Analysis 2 is correctly specified within G=1; analysis 5 gains precision by pooling but imposes an imperfect common X slope. The treatment effect remains constant within G.

## Calibration grid

Overall shares below use all five candidates. Pairwise shares compare only analyses 2 and 5.

| Effect case | gamma | tau0 | Delta | theta_5 | L2 mean | L4 mean | L5 mean | Overall 1/2/3/4/5 (%) | Pairwise 2 wins (%) |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|
| v01 | 0.00 | 0.33 | 0.17 | 0.500000 | 0.032138 | 0.025857 | 0.024421 | 0.00 / 0.00 / 0.00 / 14.33 / 85.67 | 0.00 |
| v01 | 0.05 | 0.33 | 0.17 | 0.546545 | 0.032138 | 0.024844 | 0.026604 | 0.00 / 0.00 / 0.00 / 87.67 / 12.33 | 0.00 |
| v01 | 0.10 | 0.33 | 0.17 | 0.593089 | 0.032138 | 0.023856 | 0.033153 | 0.00 / 0.00 / 0.00 / 100.00 / 0.00 | 85.33 |
| v01 | 0.15 | 0.33 | 0.17 | 0.639634 | 0.032138 | 0.022891 | 0.044069 | 0.00 / 0.00 / 0.00 / 100.00 / 0.00 | 100.00 |
| v01 | 0.20 | 0.33 | 0.17 | 0.686179 | 0.032138 | 0.021950 | 0.059351 | 0.00 / 0.00 / 0.00 / 100.00 / 0.00 | 100.00 |
| restored | 0.00 | 0.25 | 0.25 | 0.500000 | 0.032138 | 0.052560 | 0.024421 | 0.00 / 0.00 / 0.00 / 0.00 / 100.00 | 0.00 |
| restored | 0.05 | 0.25 | 0.25 | 0.546545 | 0.032138 | 0.051066 | 0.026604 | 0.00 / 0.00 / 0.00 / 0.00 / 100.00 | 0.00 |
| restored | 0.10 | 0.25 | 0.25 | 0.593089 | 0.032138 | 0.049595 | 0.033153 | 0.00 / 85.33 / 0.00 / 0.00 / 14.67 | 85.33 |
| restored | 0.15 | 0.25 | 0.25 | 0.639634 | 0.032138 | 0.048148 | 0.044069 | 0.00 / 100.00 / 0.00 / 0.00 / 0.00 | 100.00 |
| restored | 0.20 | 0.25 | 0.25 | 0.686179 | 0.032138 | 0.046725 | 0.059351 | 0.00 / 100.00 / 0.00 / 0.00 / 0.00 | 100.00 |

See `simulation_notes.md` for the calibration rationale; no target winner frequency is imposed.

## Final pilot

| Analysis | theta_a | b_a | Mean T | Mean V | Mean loss (MCSE) | Oracle optimal (MCSE, pp) |
|---|---:|---:|---:|---:|---:|---:|
| 1: Subgroup linear | 0.370005 | -0.129995 | 0.375355 | 0.031746 | 0.058810 (0.000755) | 0.00% (0.00) |
| 2: Subgroup flexible | 0.500000 | -0.000000 | 0.504775 | 0.032063 | 0.032063 (0.000079) | 87.35% (0.74) |
| 3: Subgroup flexible trimmed | 0.500000 | 0.000000 | 0.505376 | 0.034181 | 0.034181 (0.000087) | 0.00% (0.00) |
| 4: Full constant effect | 0.284201 | -0.215799 | 0.284138 | 0.002889 | 0.049532 (0.000027) | 0.00% (0.00) |
| 5: Full interaction | 0.593089 | 0.093089 | 0.595972 | 0.024360 | 0.033133 (0.000072) | 12.65% (0.74) |

## Pairwise design loss: Delta L25 = L2 - L5

Valid comparisons: 2000/2000. Negative differences favor 2; positive differences favor 5. Effective ties satisfy |Delta L25| <= 1e-10 (absolute loss units).

Mean difference: -0.001070 (MCSE 0.000023); median: -0.001119.

Strict L2<L5: 87.35%; strict L5<L2: 12.65%. With the tie tolerance: 2 wins 87.35%, 5 wins 12.65%, ties 0.0000%. Winner-share MCSE: 0.74 percentage points.

Absolute gap: mean 0.001255, median 0.001167. Relative gap |L2-L5|/min(L2,L5): mean 4.02%, median 3.68%.

Fractions with absolute gaps below 1e-6 / 1e-4 / 1e-3: 0.15% / 3.00% / 41.40%. These are descriptive scales, not selection rules.

Pairwise oracle mean loss: 0.031971. Relative to always choosing 2, the mean gain is 0.000092 (0.288% of mean L2); relative to always choosing 5 it is 0.001162 (3.508% of mean L5). This bounds the gain available from learning the realized pairwise ranking in this calibration.

| Percentile | Delta L25 | Absolute gap | Relative gap (%) |
|---:|---:|---:|---:|
| 0 | -0.00454944 | 3.41416e-07 | 0.000921 |
| 1 | -0.00340078 | 3.00907e-05 | 0.0823 |
| 5 | -0.00267608 | 0.000155105 | 0.48 |
| 25 | -0.00175109 | 0.000619363 | 1.89 |
| 50 | -0.00111932 | 0.00116665 | 3.68 |
| 75 | -0.000457655 | 0.00178081 | 5.66 |
| 95 | 0.0007 | 0.00271607 | 8.79 |
| 99 | 0.00173969 | 0.00340171 | 11.2 |
| 100 | 0.00490841 | 0.00490841 | 14.5 |

## Realized design comparison

Regions use the pairwise tolerance above, independently of which candidate wins overall. Entries are region means; CSV output also includes medians. These are descriptive associations, not causal effects.

Region sizes: 2-wins = 1747; 5-wins = 253; ties = 0.

| Feature | 2-wins | 5-wins |
|---|---:|---:|
| n_subgroup | 181.618775 | 169.525692 |
| n_subgroup_treated | 110.035489 | 104.588933 |
| n_subgroup_control | 71.583286 | 64.936759 |
| treatment_share_subgroup | 0.605820 | 0.616791 |
| ps_subgroup_q05 | 0.159150 | 0.181506 |
| ps_subgroup_q95 | 0.946195 | 0.939312 |
| ps_subgroup_spread | 0.787045 | 0.757806 |
| overlap_information_subgroup | 0.174592 | 0.178644 |
| subgroup_share_extreme_ps | 0.152654 | 0.131916 |
| ps_range_overlap_subgroup | 0.783811 | 0.723573 |
| V_2 | 0.031517 | 0.035834 |
| V_5 | 0.024113 | 0.026063 |
| conditional_bias_2 | 0.000000 | 0.000000 |
| conditional_bias_5 | 0.092687 | 0.094123 |
| variance_gain | 0.007404 | 0.009771 |
| bias_penalty | 0.008735 | 0.009040 |
| max_leverage_2 | 0.223801 | 0.227124 |
| max_leverage_5 | 0.057338 | 0.057362 |
| scaled_condition_number_2 | 3.940662 | 4.155019 |
| scaled_condition_number_5 | 3.936021 | 3.934353 |
| delta_L25 | -0.001331 | 0.000730 |
| absolute_gap | 0.001331 | 0.000730 |
| relative_gap | 0.043069 | 0.020135 |

## Checks and scope

Candidate fit failures: 0; SE failures: 0; complete oracle designs: 2000. Maximum population refinement difference: 8.33e-16.

Maximum fixed-design variance discrepancy: 0.83%; maximum absolute MSE discrepancy in Monte Carlo SE units: 1.76. All fixed-design coherence gates passed.

Analysis 3: 2000/2000 complete designs satisfy L2 <= L3 + 1e-12. Trimming remains weakly dominated under the constant-within-group effect and correct subgroup mean model.

Generated local outputs include `output/results/pairwise_designs.csv`, `output/tables/pairwise_summary.csv`, `pairwise_quantiles.csv`, `pairwise_regions.csv`, `calibration_pairwise.csv`, and `output/figures/pairwise_tradeoff.png`, alongside the existing pilot outputs. They are excluded from Git; this compact report is committed.

This is baseline DGP calibration only. No synthetic-information regimes, proposition experiments, validation budgets, or synthesis distortion are implemented.
