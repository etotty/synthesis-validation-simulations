# All numerical report content is generated from the current run and calibration.
write_pilot_report <- function(summary, pair, population, calibration, calibration_pair,
                               coherence, results) {
  p <- config$dgp
  ps <- summarize_pairwise(pair)
  quantiles <- pairwise_quantiles(pair)
  regions <- pairwise_regions(pair)
  lines <- c("# Generated v0.1.1 pilot results", "",
    "Regenerate with `Rscript scripts/run_all.R`. Source configuration: `R/config.R`.", "",
    "## Run and interpretation", "",
    sprintf("Pilot: %d designs; calibration: %d designs per case; fixed-design oracle check: %d noise draws.",
            config$pilot_reps, config$calibration_reps, config$oracle_check_reps), "",
    sprintf("n=%d; p_G=%.2f; alpha=(%s); beta=(%s); gamma=%.2f; tau0=%.2f; Delta=%.2f; sigma=%.2f; theta=%.2f.",
      p$n, p$p_G, paste(p$alpha, collapse=", "), paste(p$beta, collapse=", "),
      p$gamma, p$tau0, p$delta, p$sigma, target_effect(p)), "",
    sprintf("Known-propensity trimming: [%.2f, %.2f]. RNG: %s. Calibration seed = %d+r; pilot seed = %d+r. The same calibration seeds are used in every case.",
      config$trim[1], config$trim[2], paste(config$rng, collapse=" / "),
      config$seeds$calibration, config$seeds$pilot), "",
    "The outcome includes gamma*G*X, while the five candidate regressions are unchanged. Analysis 2 is correctly specified within G=1; analysis 5 gains precision by pooling but imposes an imperfect common X slope. The treatment effect remains constant within G.", "",
    "## Calibration grid", "",
    "Overall shares below use all five candidates. Pairwise shares compare only analyses 2 and 5.", "",
    "| Effect case | gamma | tau0 | Delta | theta_5 | L2 mean | L4 mean | L5 mean | Overall 1/2/3/4/5 (%) | Pairwise 2 wins (%) |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---|---:|")
  for (case in names(config$calibration_effects)) {
    for (gamma in config$gamma_grid) {
      d <- calibration[calibration$effect_case == case & calibration$gamma == gamma, ]
      d <- d[order(d$analysis), ]
      cp <- calibration_pair[calibration_pair$effect_case == case & calibration_pair$gamma == gamma, ]
      stopifnot(nrow(d) == 5L, nrow(cp) == 1L)
      lines <- c(lines, sprintf("| %s | %.2f | %.2f | %.2f | %.6f | %.6f | %.6f | %.6f | %s | %.2f |",
        case, gamma, d$tau0[1], d$delta[1], d$theta_a[5], d$mean_loss[2], d$mean_loss[4], d$mean_loss[5],
        paste(sprintf("%.2f", 100*d$oracle_frequency), collapse=" / "), 100*cp$fraction_2_wins))
    }
  }
  lines <- c(lines, "", "See `simulation_notes.md` for the calibration rationale; no target winner frequency is imposed.", "",
    "## Final pilot", "",
    "| Analysis | theta_a | b_a | Mean T | Mean V | Mean loss (MCSE) | Oracle optimal (MCSE, pp) |",
    "|---|---:|---:|---:|---:|---:|---:|")
  for (i in 1:5) lines <- c(lines, sprintf("| %d: %s | %.6f | %.6f | %.6f | %.6f | %.6f (%.6f) | %.2f%% (%.2f) |",
    i, analysis_labels[i], summary$theta_a[i], summary$b_a[i], summary$mean_T[i], summary$mean_V[i],
    summary$mean_loss[i], summary$mean_loss_mcse[i], 100*summary$oracle_frequency[i],
    100*summary$oracle_frequency_mcse[i]))
  lines <- c(lines, "", "## Pairwise design loss: Delta L25 = L2 - L5", "",
    sprintf("Valid comparisons: %d/%d. Negative differences favor 2; positive differences favor 5. Effective ties satisfy |Delta L25| <= %.0e (absolute loss units).", ps$n_valid, ps$n_designs, ps$tie_tolerance), "",
    sprintf("Mean difference: %.6f (MCSE %.6f); median: %.6f.", ps$mean_delta, ps$mean_delta_mcse, ps$median_delta), "",
    sprintf("Strict L2<L5: %.2f%%; strict L5<L2: %.2f%%. With the tie tolerance: 2 wins %.2f%%, 5 wins %.2f%%, ties %.4f%%. Winner-share MCSE: %.2f percentage points.",
      100*ps$fraction_2_strict, 100*ps$fraction_5_strict, 100*ps$fraction_2_wins,
      100*ps$fraction_5_wins, 100*ps$fraction_tied, 100*ps$fraction_2_mcse), "",
    sprintf("Absolute gap: mean %.6f, median %.6f. Relative gap |L2-L5|/min(L2,L5): mean %.2f%%, median %.2f%%.",
      ps$mean_absolute_gap, ps$median_absolute_gap, 100*ps$mean_relative_gap, 100*ps$median_relative_gap), "",
    sprintf("Fractions with absolute gaps below 1e-6 / 1e-4 / 1e-3: %.2f%% / %.2f%% / %.2f%%. These are descriptive scales, not selection rules.",
      100*ps$fraction_gap_below_1e6, 100*ps$fraction_gap_below_1e4, 100*ps$fraction_gap_below_1e3), "",
    sprintf("Pairwise oracle mean loss: %.6f. Relative to always choosing 2, the mean gain is %.6f (%.3f%% of mean L2); relative to always choosing 5 it is %.6f (%.3f%% of mean L5). This bounds the gain available from learning the realized pairwise ranking in this calibration.",
      mean(pmin(pair$L_2[pair$valid], pair$L_5[pair$valid])),
      mean(pmax(pair$delta_L25[pair$valid], 0)),
      100*mean(pmax(pair$delta_L25[pair$valid], 0))/mean(pair$L_2[pair$valid]),
      mean(pmax(-pair$delta_L25[pair$valid], 0)),
      100*mean(pmax(-pair$delta_L25[pair$valid], 0))/mean(pair$L_5[pair$valid])), "",
    "| Percentile | Delta L25 | Absolute gap | Relative gap (%) |",
    "|---:|---:|---:|---:|")
  for (i in seq_len(nrow(quantiles))) lines <- c(lines, sprintf("| %.0f | %.6g | %.6g | %.3g |",
    100*quantiles$probability[i], quantiles$delta_L25[i], quantiles$absolute_gap[i], 100*quantiles$relative_gap[i]))
  lines <- c(lines, "", "## Realized design comparison", "",
    "Regions use the pairwise tolerance above, independently of which candidate wins overall. Entries are region means; CSV output also includes medians. These are descriptive associations, not causal effects.", "",
    sprintf("Region sizes: 2-wins = %d; 5-wins = %d; ties = %d.", sum(pair$region=="2_wins"), sum(pair$region=="5_wins"), sum(pair$region=="tie")), "",
    "| Feature | 2-wins | 5-wins |", "|---|---:|---:|")
  features <- unique(regions$feature)
  for (feature in features) {
    d <- regions[regions$feature==feature, ]
    lines <- c(lines, sprintf("| %s | %.6f | %.6f |", feature,
      d$mean[d$region=="2_wins"], d$mean[d$region=="5_wins"]))
  }
  L <- matrix(results$L, ncol=5, byrow=TRUE)
  lines <- c(lines, "", "## Checks and scope", "",
    sprintf("Candidate fit failures: %d; SE failures: %d; complete oracle designs: %d. Maximum population refinement difference: %.2e.",
      sum(summary$failures), sum(summary$se_failures), summary$valid_oracle_designs[1], max(abs(population$refinement_difference))), "",
    sprintf("Maximum fixed-design variance discrepancy: %.2f%%; maximum absolute MSE discrepancy in Monte Carlo SE units: %.2f. All fixed-design coherence gates passed.",
      100*max(abs(coherence$simulated_variance/coherence$analytic_variance-1)),
      max(abs(coherence$simulated_mse-coherence$analytic_loss)/coherence$mse_mcse)), "",
    sprintf("Analysis 3: %d/%d complete designs satisfy L2 <= L3 + 1e-12. Trimming remains weakly dominated under the constant-within-group effect and correct subgroup mean model.",
      sum(L[,2] <= L[,3]+1e-12, na.rm=TRUE), sum(complete.cases(L))), "",
    "Generated local outputs include `output/results/pairwise_designs.csv`, `output/tables/pairwise_summary.csv`, `pairwise_quantiles.csv`, `pairwise_regions.csv`, `calibration_pairwise.csv`, and `output/figures/pairwise_tradeoff.png`, alongside the existing pilot outputs. They are excluded from Git; this compact report is committed.", "",
    "This is baseline DGP calibration only. No synthetic-information regimes, proposition experiments, validation budgets, or synthesis distortion are implemented.")
  writeLines(lines, "paper/pilot_results.md")
}
