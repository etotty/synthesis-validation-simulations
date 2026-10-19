write_researcher_report <- function(policy,diagnostic,confusion) {
  p <- policy; d <- diagnostic
  assessment <- vapply(config$researcher$cases$case,function(case) {
    x <- p[p$case==case & p$policy=="feasible",]
    sprintf("%s: feasible selection changes assigned loss relative to the best fixed rule by %+.6f (MCSE %.6f); it captures %.1f%% of the attainable gain. The paired pure-noise shift relative to always 2 is %+.6f (MCSE %.6f).",
      case,-x$improvement_over_fixed,x$improvement_mcse,100*x$fraction_captured,
      x$noise_selection_vs_2,x$noise_selection_mcse)
  },character(1))
  lines <- c("# v0.2 researcher-observable design information: generated results", "",
    "Regenerate with `Rscript scripts/04_researcher_information.R`. See `researcher_information_model.md` for the maintained assumptions, information boundary, algebra, and interpretation.", "",
    sprintf("%d replications per case. Canonical n=1500, p_G=0.12; sensitivity n=2000, p_G=0.08. All remaining parameters are unchanged. Common draws use seed %d+r, maximum length 2000. Independent-noise evaluation uses seed %d+r. RNG: %s.",
      config$researcher$reps,config$seeds$researcher,config$seeds$researcher_fresh_noise,paste(config$rng,collapse=" / ")), "",
    "The researcher sees exactly Y, X, D, G. Baseline variance proxies use the augmented model's residual variance times each candidate's squared weight norm. Estimated restriction penalty is [gamma_hat*h5'(GX)]^2. Analysis 2 has no assigned restriction penalty by a maintained reference-model assumption. Ties within 1e-10 favor 2. HC3 is a diagnostic comparison, not a replacement rule. True propensity, gamma, sigma, m, theta_a, biases and losses enter evaluation only.", "",
    "## Assessment for review", "", paste(assessment,collapse="\n\n"), "",
    "The benchmark separates noisy structural evidence from actual design quality. Interpret the raw loss changes before the normalized fractions. Pure-noise selection is assessed using paired Monte Carlo errors; a small net treatment bias can conceal offsetting specification and noise distortions. HC3 is reported without selecting a preferred rule after seeing results. The conceptual document records the decisions needed before extending this layer.", "",
    "## Design-loss performance", "",
    "D is the mean assigned oracle design loss E[L_choice(Z)]. Because the feasible choice depends on Y, D is NOT generally the actual MSE of its same-sample treatment estimate. The latter is reported separately below. The oracle comparator is exact pairwise min loss. Best fixed means the smaller of the two unconditional mean losses, not a choice made separately in every replication.", "",
    design_value_table(p[,c("case","policy","choice_5","agreement","D_feasible","D_feasible_mcse",
      "D_fixed","D_oracle","attainable_gain","attainable_gain_mcse","excess_over_oracle","improvement_over_fixed",
      "improvement_mcse","fraction_captured","fraction_mcse","denominator_fragile")]), "",
    "All shares/fractions above are on a 0-1 scale, not percentages. Fraction captured = (D_fixed-D_feasible)/(D_fixed-D_oracle), with no clipping; zero denominator returns NA. MCSE for this ratio uses paired influence values and is branch-local. The small denominator amplifies uncertainty; raw losses/gains take priority. For a coherent same-draw comparison, assigned loss is pointwise at least the oracle minimum, so a fraction above one would indicate inconsistent evaluation, not ordinary Monte Carlo variation. Negative fractions are allowed and mean the rule is worse than the best fixed choice.", "",
    "## Confusion matrix", "", design_value_table(confusion), "",
    "Counts are oracle rows versus feasible columns; fractions divide by all replications in the case. Agreement alone does not measure loss performance.", "",
    "## Systematic and realization-specific consequences", "",
    "epsilon_a=T_a-theta_a includes both a realized-design offset E[T_a|Z]-theta_a and pure outcome noise eta_a=T_a-E[T_a|Z]. Report all three rather than interpreting mean epsilon alone as result-noise selection. Systematic + epsilon distortion equals total bias by construction. Noise selection versus 2 is the paired mean eta_choice-eta_2, which removes the baseline's Monte Carlo mean-noise error.", "",
    design_value_table(p[,c("case","policy","systematic_distortion","systematic_mcse",
      "epsilon_distortion","epsilon_mcse","design_offset","outcome_noise_distortion",
      "outcome_noise_mcse","noise_selection_vs_2","noise_selection_mcse","total_bias")]), "",
    "## Same-sample MSE versus assigned design loss", "",
    "Independent-noise evaluation keeps the realized design and the already-made choice fixed, then evaluates the treatment estimate using a fresh independent error vector. It is a simulation diagnostic, not a feasible sample-splitting procedure or a new information regime. Its expected MSE equals the assigned oracle loss. Same-sample MSE need not equal that loss because choice and treatment errors share outcome noise.", "",
    design_value_table(p[,c("case","policy","D_feasible","actual_mse","actual_mse_mcse",
      "independent_noise_mse","independent_noise_mse_mcse","actual_mse_minus_design_loss",
      "mse_minus_loss_mcse","fresh_noise_distortion","fresh_noise_mcse")]), "",
    "## Noisy structural evidence and precision", "",
    design_value_table(d[,c("case","B","gamma_mean","gamma_mean_mcse","gamma_sd","gamma_rmse",
      "gamma_q025","gamma_median","gamma_q975","gamma_mean_se","gamma_negative","gamma_coverage",
      "sigma2_mean")]), "",
    "True gamma=0.10 appears here only as an evaluation reference. Coverage uses a descriptive normal 1.96 interval. The squared plug-in penalty includes estimation noise: under the maintained augmented model, its conditional expectation is b5(Z)^2 + c5(Z)^2 Var(gamma_hat|Z). No de-biasing or shrinkage is silently applied.", "",
    design_value_table(d[,c("case","bias_hat_mean","oracle_bias_mean","bias_tracking_correlation",
      "bias_rmse","loss_gap_correlation","mean_estimated_penalty","mean_oracle_penalty",
      "mean_expected_penalty_inflation","mean_variance_proxy_2","mean_variance_proxy_5",
      "mean_oracle_variance_2","mean_oracle_variance_5")]), "",
    "## Result-selection correlations", "",
    "Choice indicator is 1 for analysis 5. Correlations with eta isolate pure outcome noise from realized-design offsets. The pooled treatment estimate and gamma_hat have zero conditional noise covariance under the augmented Gaussian model; the subgroup treatment estimate generally does not. The squared-gamma choice can therefore select treatment noise without directly comparing treatment coefficients.", "",
    design_value_table(d[,c("case","choice_cor_epsilon_2","choice_cor_epsilon_5","choice_cor_noise_2",
      "choice_cor_noise_5","choice_cor_noise_difference","gamma_error_cor_noise_2",
      "gamma_error_cor_noise_5")]), "",
    "## Design and algebra checks", "",
    design_value_table(d[,c("case","subgroup_mean","subgroup_q05","subgroup_q95","treated_mean",
      "control_mean","max_ovb_error","max_bias_identity_error","feasible_ties")]), "",
    "All reported cases completed without coefficient-fit or feasible-precision failure; such failures stop this benchmark rather than silently change its denominator. The full five-analysis estimates, HC3 SEs, population decompositions, oracle moments, overlap and conditioning diagnostics remain in the local menu draws and menu-summary tables. True propensity diagnostics are evaluation-only and never cross the observed-data interface.", "",
    "Detailed generated files are `output/results/researcher_*_draws.csv`, `researcher_*_menu.csv`, configuration/session records, and `output/tables/researcher_*.csv`. The figure is `output/figures/researcher_information.png`. Only this compact generated report is tracked; detailed outputs remain ignored. No DGP recalibration, synthetic datasets, fidelity/leakage/sensitivity parameters, regimes, validation, or synthesis distortion are implemented.")
  writeLines(lines,"paper/researcher_information_results.md")
}

plot_researcher_results <- function(policy) {
  png("output/figures/researcher_information.png",width=1200,height=650,res=125)
  on.exit(dev.off())
  par(mfrow=c(1,2),mar=c(6,4,3,1))
  for (case in config$researcher$cases$case) {
    p <- policy[policy$case==case,]
    p <- p[match(c("always_2","always_5","oracle_design","feasible","HC3_diagnostic"),p$policy),]
    values <- rbind(p$D_feasible,p$actual_mse,p$independent_noise_mse)
    colors <- c("#2166ac","#b2182b","#4d9221")
    barplot(values,beside=TRUE,col=colors,names.arg=c("Always 2","Always 5","Oracle","Feasible","HC3"),
            las=2,ylab="Mean squared-error loss",main=case,ylim=c(0,max(values)*1.35))
    legend("top",legend=c("Assigned design loss","Same-sample MSE","Independent-noise MSE"),
           fill=colors,bty="n",cex=.8)
  }
}
