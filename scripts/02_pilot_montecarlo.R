source("R/simulation_helpers.R")
load_project()
ensure_output()
# Coherence gate MUST pass before the pilot runs.
check <- check_oracle()
stopifnot(all(abs(check$simulated_variance - check$analytic_variance) < 6 * check$variance_mcse),
          all(abs(check$simulated_mean - check$expected) < 6 * check$mean_mcse),
          all(abs(check$simulated_mse - check$analytic_loss) < 6 * check$mse_mcse))
write_csv(check, "output/tables/oracle_coherence.csv")
source("scripts/01_population_estimands.R")
results <- run_montecarlo(config$pilot_reps, config$dgp, population$theta_a,
                         config$seeds$pilot, progress = TRUE)
write_csv(results, "output/results/pilot_draws.csv")
summary <- summarize_pilot(results)
write_csv(summary, "output/tables/pilot_summary.csv")
# Entire joint rows only: no pairwise deletion that could obscure failures.
epsilon <- matrix(results$epsilon, ncol = 5, byrow = TRUE)
colnames(epsilon) <- paste0("analysis_", 1:5)
covariance <- cov(epsilon, use = "complete.obs")
write_csv(data.frame(analysis = 1:5, covariance), "output/tables/epsilon_covariance.csv")
metrics <- c("n_subgroup", "n_subgroup_treated", "n_subgroup_control", "n_used", "n_trimmed",
             "treatment_share_subgroup", "ps_subgroup_q05", "ps_subgroup_q95",
             "overlap_information_subgroup", "overlap_information_other",
             "subgroup_share_extreme_ps", "ps_range_overlap_subgroup",
             "max_leverage", "scaled_condition_number")
diagnostic_table <- do.call(rbind, lapply(1:5, function(a) {
  d <- results[results$analysis == a, ]
  do.call(rbind, lapply(metrics, function(metric) {
    x <- d[[metric]]
    q <- quantile(x, c(0, .05, .5, .95, 1), na.rm = TRUE)
    data.frame(analysis = a, metric = metric, mean = mean(x, na.rm = TRUE),
               min = q[1], q05 = q[2], median = q[3], q95 = q[4], max = q[5], row.names = NULL)
  }))
}))
write_csv(diagnostic_table, "output/tables/design_diagnostics.csv")
# PNGs have no variable creation-date metadata and can be hash-checked on rerun.
png("output/figures/estimates_and_loss.png", width = 1400, height = 650, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
boxplot(T ~ analysis, data = results, outline = FALSE, xlab = "Analysis", ylab = "Treatment estimate",
        main = sprintf("Estimates across %s samples", format(config$pilot_reps, big.mark = ",")), col = "lightblue")
abline(h = target_effect(config$dgp), col = "firebrick", lty = 2)
points(1:5, population$theta_a, pch = 18, col = "navy", cex = 1.4)
legend("topright", c("Common target", "Population counterpart"), col = c("firebrick", "navy"),
       lty = c(2, NA), pch = c(NA, 18), bty = "n", cex = .8)
boxplot(L ~ analysis, data = results, outline = FALSE, xlab = "Analysis", ylab = "Conditional MSE",
        main = "Oracle design loss", col = "wheat")
dev.off()
png("output/figures/design_diagnostics.png", width = 1400, height = 1000, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
one <- results[results$analysis == 1, ]
trimmed <- results[results$analysis == 3, ]
hist(one$n_subgroup, main = "Subgroup size", xlab = "Number of G = 1 observations", col = "lightblue")
hist(one$treatment_share_subgroup, main = "Subgroup treatment share", xlab = "Treated fraction", col = "lightblue")
hist(trimmed$n_trimmed, main = "Overlap trimming", xlab = "Subgroup observations removed", col = "wheat")
plot(one$ps_subgroup_q05, one$ps_subgroup_q95, pch = 16, cex = .35,
     col = adjustcolor("navy", .2), xlab = "Subgroup PS: 5th percentile",
     ylab = "Subgroup PS: 95th percentile", main = "Realized propensity spread")
dev.off()
# Generated review summary: values never copied manually into tracked prose.
lines <- c("# Generated v0.1 pilot results", "", "Regenerate with `Rscript scripts/run_all.R`.", "",
           "| Analysis | theta_a | b_a | Mean loss | Oracle optimal |",
           "|---|---:|---:|---:|---:|")
for (i in 1:5) lines <- c(lines, sprintf("| %d: %s | %.6f | %.6f | %.6f | %.2f%% |", i,
  analysis_labels[i], summary$theta_a[i], summary$b_a[i], summary$mean_loss[i], 100 * summary$oracle_frequency[i]))
lines <- c(lines, "", sprintf("Failures: %d candidate fits. Complete oracle designs: %d.",
                              sum(summary$failures), summary$valid_oracle_designs[1]), "",
           "The raw draws, covariance matrix, Monte Carlo standard errors, and diagnostic quantiles are in output/.")
writeLines(lines, "paper/pilot_results.md")
capture.output(sessionInfo(), file = "output/results/session_info.txt")
capture.output(dput(config), file = "output/results/run_config.R")
print(summary, digits = 5, row.names = FALSE)
