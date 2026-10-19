# v0.2 researcher-observable design information: generated results

Regenerate with `Rscript scripts/04_researcher_information.R`. See `researcher_information_model.md` for the maintained assumptions, information boundary, algebra, and interpretation.

10000 replications per case. Canonical n=1500, p_G=0.12; sensitivity n=2000, p_G=0.08. All remaining parameters are unchanged. Common draws use seed 700000+r, maximum length 2000. Independent-noise evaluation uses seed 800000+r. RNG: Mersenne-Twister / Inversion / Rejection.

The researcher sees exactly Y, X, D, G. Baseline variance proxies use the augmented model's residual variance times each candidate's squared weight norm. Estimated restriction penalty is [gamma_hat*h5'(GX)]^2. Analysis 2 has no assigned restriction penalty by a maintained reference-model assumption. Ties within 1e-10 favor 2. HC3 is a diagnostic comparison, not a replacement rule. True propensity, gamma, sigma, m, theta_a, biases and losses enter evaluation only.

## Assessment for review

canonical: feasible selection changes assigned loss relative to the best fixed rule by +0.000457 (MCSE 0.000008); it captures -528.9% of the attainable gain. The paired pure-noise shift relative to always 2 is -0.030654 (MCSE 0.000465).

sensitivity: feasible selection changes assigned loss relative to the best fixed rule by +0.000155 (MCSE 0.000010); it captures -45.4% of the attainable gain. The paired pure-noise shift relative to always 2 is -0.032602 (MCSE 0.000500).

The benchmark separates noisy structural evidence from actual design quality. Interpret the raw loss changes before the normalized fractions. Pure-noise selection is assessed using paired Monte Carlo errors; a small net treatment bias can conceal offsetting specification and noise distortions. HC3 is reported without selecting a preferred rule after seeing results. The conceptual document records the decisions needed before extending this layer.

## Design-loss performance

D is the mean assigned oracle design loss E[L_choice(Z)]. Because the feasible choice depends on Y, D is NOT generally the actual MSE of its same-sample treatment estimate. The latter is reported separately below. The oracle comparator is exact pairwise min loss. Best fixed means the smaller of the two unconditional mean losses, not a choice made separately in every replication.

| case | policy | choice_5 | agreement | D_feasible | D_feasible_mcse | D_fixed | D_oracle | attainable_gain | attainable_gain_mcse | excess_over_oracle | improvement_over_fixed | improvement_mcse | fraction_captured | fraction_mcse | denominator_fragile |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | always_2 | 0.000000 | 0.870500 | 0.032103 | 0.000034 | 0.032103 | 0.032017 | 0.000086 | 0.000003 | 0.000086 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | FALSE |
| canonical | always_5 | 1.000000 | 0.129500 | 0.033174 | 0.000031 | 0.032103 | 0.032017 | 0.000086 | 0.000003 | 0.001158 | -0.001071 | 0.000010 | -12.395041 | 0.548109 | FALSE |
| canonical | oracle_design | 0.129500 | 1.000000 | 0.032017 | 0.000033 | 0.032103 | 0.032017 | 0.000086 | 0.000003 | 0.000000 | 0.000086 | 0.000003 | 1.000000 | 0.000000 | FALSE |
| canonical | feasible | 0.454400 | 0.543900 | 0.032560 | 0.000034 | 0.032103 | 0.032017 | 0.000086 | 0.000003 | 0.000543 | -0.000457 | 0.000008 | -5.288508 | 0.251489 | FALSE |
| canonical | HC3_diagnostic | 0.459000 | 0.540700 | 0.032566 | 0.000034 | 0.032103 | 0.032017 | 0.000086 | 0.000003 | 0.000549 | -0.000463 | 0.000009 | -5.353918 | 0.254464 | FALSE |
| sensitivity | always_2 | 0.000000 | 0.664000 | 0.036190 | 0.000042 | 0.036190 | 0.035850 | 0.000340 | 0.000008 | 0.000340 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | FALSE |
| sensitivity | always_5 | 1.000000 | 0.336000 | 0.036601 | 0.000036 | 0.036190 | 0.035850 | 0.000340 | 0.000008 | 0.000751 | -0.000411 | 0.000013 | -1.207401 | 0.063170 | FALSE |
| sensitivity | oracle_design | 0.336000 | 1.000000 | 0.035850 | 0.000038 | 0.036190 | 0.035850 | 0.000340 | 0.000008 | 0.000000 | 0.000340 | 0.000008 | 1.000000 | 0.000000 | FALSE |
| sensitivity | feasible | 0.472500 | 0.527500 | 0.036345 | 0.000039 | 0.036190 | 0.035850 | 0.000340 | 0.000008 | 0.000495 | -0.000155 | 0.000010 | -0.454216 | 0.035619 | FALSE |
| sensitivity | HC3_diagnostic | 0.470100 | 0.532700 | 0.036342 | 0.000039 | 0.036190 | 0.035850 | 0.000340 | 0.000008 | 0.000492 | -0.000152 | 0.000010 | -0.446946 | 0.034974 | FALSE |

All shares/fractions above are on a 0-1 scale, not percentages. Fraction captured = (D_fixed-D_feasible)/(D_fixed-D_oracle), with no clipping; zero denominator returns NA. MCSE for this ratio uses paired influence values and is branch-local. The small denominator amplifies uncertainty; raw losses/gains take priority. For a coherent same-draw comparison, assigned loss is pointwise at least the oracle minimum, so a fraction above one would indicate inconsistent evaluation, not ordinary Monte Carlo variation. Negative fractions are allowed and mean the rule is worse than the best fixed choice.

## Confusion matrix

| case | oracle | feasible | count | fraction |
| --- | --- | --- | --- | --- |
| canonical | 2 | 2 | 4800 | 0.480000 |
| canonical | 5 | 2 | 656 | 0.065600 |
| canonical | 2 | 5 | 3905 | 0.390500 |
| canonical | 5 | 5 | 639 | 0.063900 |
| sensitivity | 2 | 2 | 3595 | 0.359500 |
| sensitivity | 5 | 2 | 1680 | 0.168000 |
| sensitivity | 2 | 5 | 3045 | 0.304500 |
| sensitivity | 5 | 5 | 1680 | 0.168000 |

Counts are oracle rows versus feasible columns; fractions divide by all replications in the case. Agreement alone does not measure loss performance.

## Systematic and realization-specific consequences

epsilon_a=T_a-theta_a includes both a realized-design offset E[T_a|Z]-theta_a and pure outcome noise eta_a=T_a-E[T_a|Z]. Report all three rather than interpreting mean epsilon alone as result-noise selection. Systematic + epsilon distortion equals total bias by construction. Noise selection versus 2 is the paired mean eta_choice-eta_2, which removes the baseline's Monte Carlo mean-noise error.

| case | policy | systematic_distortion | systematic_mcse | epsilon_distortion | epsilon_mcse | design_offset | outcome_noise_distortion | outcome_noise_mcse | noise_selection_vs_2 | noise_selection_mcse | total_bias |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | always_2 | -8.88e-16 | 0.000000 | -0.000833 | 0.001815 | 8.77e-16 | -0.000833 | 0.001815 | 0.000000 | 0.000000 | -0.000833 |
| canonical | always_5 | 0.093089 | 0.000000 | -0.001012 | 0.001589 | 0.000052 | -0.001064 | 0.001586 | -0.000231 | 0.000887 | 0.092077 |
| canonical | oracle_design | 0.012055 | 0.000313 | -0.000730 | 0.001780 | 0.000174 | -0.000904 | 0.001780 | -0.000071 | 0.000357 | 0.011325 |
| canonical | feasible | 0.042300 | 0.000464 | -0.031456 | 0.001701 | 0.000031 | -0.031486 | 0.001699 | -0.030654 | 0.000465 | 0.010844 |
| canonical | HC3_diagnostic | 0.042728 | 0.000464 | -0.030523 | 0.001704 | 0.000028 | -0.030551 | 0.001703 | -0.029718 | 0.000469 | 0.012205 |
| sensitivity | always_2 | 1.67e-15 | 0.000000 | -0.000920 | 0.001897 | -1.62e-15 | -0.000920 | 0.001897 | 0.000000 | 0.000000 | -0.000920 |
| sensitivity | always_5 | 0.096414 | 0.000000 | -0.000770 | 0.001655 | 0.000018 | -0.000788 | 0.001650 | 0.000132 | 0.000950 | 0.095644 |
| sensitivity | oracle_design | 0.032395 | 0.000455 | -0.001258 | 0.001813 | 0.000153 | -0.001411 | 0.001811 | -0.000491 | 0.000596 | 0.031137 |
| sensitivity | feasible | 0.045556 | 0.000481 | -0.033451 | 0.001777 | 0.000071 | -0.033522 | 0.001774 | -0.032602 | 0.000500 | 0.012105 |
| sensitivity | HC3_diagnostic | 0.045324 | 0.000481 | -0.031993 | 0.001781 | 0.000041 | -0.032034 | 0.001779 | -0.031114 | 0.000500 | 0.013331 |

## Same-sample MSE versus assigned design loss

Independent-noise evaluation keeps the realized design and the already-made choice fixed, then evaluates the treatment estimate using a fresh independent error vector. It is a simulation diagnostic, not a feasible sample-splitting procedure or a new information regime. Its expected MSE equals the assigned oracle loss. Same-sample MSE need not equal that loss because choice and treatment errors share outcome noise.

| case | policy | D_feasible | actual_mse | actual_mse_mcse | independent_noise_mse | independent_noise_mse_mcse | actual_mse_minus_design_loss | mse_minus_loss_mcse | fresh_noise_distortion | fresh_noise_mcse |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | always_2 | 0.032103 | 0.032922 | 0.000467 | 0.032671 | 0.000464 | 0.000819 | 0.000465 | -0.000845 | 0.001808 |
| canonical | always_5 | 0.033174 | 0.033738 | 0.000464 | 0.033315 | 0.000454 | 0.000564 | 0.000463 | -0.000158 | 0.001565 |
| canonical | oracle_design | 0.032017 | 0.032923 | 0.000465 | 0.032419 | 0.000454 | 0.000906 | 0.000464 | -0.000819 | 0.001768 |
| canonical | feasible | 0.032560 | 0.033701 | 0.000473 | 0.033060 | 0.000462 | 0.001141 | 0.000471 | -0.000679 | 0.001704 |
| canonical | HC3_diagnostic | 0.032566 | 0.033739 | 0.000474 | 0.032882 | 0.000457 | 0.001173 | 0.000473 | -0.000486 | 0.001700 |
| sensitivity | always_2 | 0.036190 | 0.035975 | 0.000525 | 0.037003 | 0.000521 | -0.000215 | 0.000524 | -0.000119 | 0.001924 |
| sensitivity | always_5 | 0.036601 | 0.036550 | 0.000505 | 0.036875 | 0.000508 | -0.000051 | 0.000503 | 0.000182 | 0.001654 |
| sensitivity | oracle_design | 0.035850 | 0.036076 | 0.000510 | 0.036443 | 0.000499 | 0.000226 | 0.000508 | -0.000096 | 0.001820 |
| sensitivity | feasible | 0.036345 | 0.036905 | 0.000526 | 0.037027 | 0.000514 | 0.000561 | 0.000525 | -0.000330 | 0.001802 |
| sensitivity | HC3_diagnostic | 0.036342 | 0.036929 | 0.000528 | 0.036795 | 0.000512 | 0.000586 | 0.000526 | -0.000114 | 0.001798 |

## Noisy structural evidence and precision

| case | B | gamma_mean | gamma_mean_mcse | gamma_sd | gamma_rmse | gamma_q025 | gamma_median | gamma_q975 | gamma_mean_se | gamma_negative | gamma_coverage | sigma2_mean |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | 10000 | 0.099447 | 0.000923 | 0.092308 | 0.092305 | -0.079912 | 0.099593 | 0.280354 | 0.091741 | 0.140000 | 0.949000 | 1.000664 |
| sensitivity | 10000 | 0.099678 | 0.000957 | 0.095746 | 0.095742 | -0.087481 | 0.100492 | 0.286218 | 0.095644 | 0.151100 | 0.952500 | 1.000589 |

True gamma=0.10 appears here only as an evaluation reference. Coverage uses a descriptive normal 1.96 interval. The squared plug-in penalty includes estimation noise: under the maintained augmented model, its conditional expectation is b5(Z)^2 + c5(Z)^2 Var(gamma_hat|Z). No de-biasing or shrinkage is silently applied.

| case | bias_hat_mean | oracle_bias_mean | bias_tracking_correlation | bias_rmse | loss_gap_correlation | mean_estimated_penalty | mean_oracle_penalty | mean_expected_penalty_inflation | mean_variance_proxy_2 | mean_variance_proxy_5 | mean_oracle_variance_2 | mean_oracle_variance_5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | 0.092651 | 0.093142 | 0.137506 | 0.086945 | 0.041524 | 0.016289 | 0.008820 | 0.007458 | 0.032126 | 0.024371 | 0.032103 | 0.024354 |
| sensitivity | 0.096235 | 0.096432 | 0.149032 | 0.093204 | 0.027646 | 0.018145 | 0.009473 | 0.008707 | 0.036212 | 0.027144 | 0.036190 | 0.027128 |

## Result-selection correlations

Choice indicator is 1 for analysis 5. Correlations with eta isolate pure outcome noise from realized-design offsets. The pooled treatment estimate and gamma_hat have zero conditional noise covariance under the augmented Gaussian model; the subgroup treatment estimate generally does not. The squared-gamma choice can therefore select treatment noise without directly comparing treatment coefficients.

| case | choice_cor_epsilon_2 | choice_cor_epsilon_5 | choice_cor_noise_2 | choice_cor_noise_5 | choice_cor_noise_difference | gamma_error_cor_noise_2 | gamma_error_cor_noise_5 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | 0.334628 | -0.003922 | 0.334628 | -0.004018 | -0.691692 | -0.472741 | 0.003594 |
| sensitivity | 0.337393 | -0.007916 | 0.337393 | -0.008705 | -0.689063 | -0.476529 | 0.011041 |

## Design and algebra checks

| case | subgroup_mean | subgroup_q05 | subgroup_q95 | treated_mean | control_mean | max_ovb_error | max_bias_identity_error | feasible_ties |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| canonical | 180.296600 | 159 | 201 | 109.734200 | 70.562400 | 5.22e-15 | 1.15e-14 | 0 |
| sensitivity | 160.314800 | 140 | 180 | 97.561500 | 62.753300 | 4.16e-15 | 1.69e-14 | 0 |

All reported cases completed without coefficient-fit or feasible-precision failure; such failures stop this benchmark rather than silently change its denominator. The full five-analysis estimates, HC3 SEs, population decompositions, oracle moments, overlap and conditioning diagnostics remain in the local menu draws and menu-summary tables. True propensity diagnostics are evaluation-only and never cross the observed-data interface.

Detailed generated files are `output/results/researcher_*_draws.csv`, `researcher_*_menu.csv`, configuration/session records, and `output/tables/researcher_*.csv`. The figure is `output/figures/researcher_information.png`. Only this compact generated report is tracked; detailed outputs remain ignored. No DGP recalibration, synthetic datasets, fidelity/leakage/sensitivity parameters, regimes, validation, or synthesis distortion are implemented.
