source("R/simulation_helpers.R")
load_project()
ensure_output()
# The same replication seeds are reused for every gamma/effect case. Only the
# explicitly listed outcome coefficients change; X,D,G and u stay identical.
calibration_summaries <- list()
calibration_pairs <- list()
for (effect_case in names(config$calibration_effects)) {
  for (gamma in config$gamma_grid) {
    parameters <- config$dgp
    parameters$gamma <- gamma
    parameters$tau0 <- unname(config$calibration_effects[[effect_case]]["tau0"])
    parameters$delta <- unname(config$calibration_effects[[effect_case]]["delta"])
    population <- population_estimands(parameters)
    calibration <- run_montecarlo(config$calibration_reps, parameters, population$theta_a,
                                 config$seeds$calibration)
    summary <- summarize_pilot(calibration)
    summary$effect_case <- effect_case
    summary$gamma <- gamma
    summary$tau0 <- parameters$tau0
    summary$delta <- parameters$delta
    pair_summary <- summarize_pairwise(pairwise_designs(calibration))
    pair_summary$effect_case <- effect_case
    pair_summary$gamma <- gamma
    key <- paste(effect_case, gamma, sep = "_")
    calibration_summaries[[key]] <- summary
    calibration_pairs[[key]] <- pair_summary
    message("Calibration ", effect_case, ", gamma=", gamma, ": pairwise 2 wins ",
            round(100 * pair_summary$fraction_2_wins, 1), "%; overall 2/5 wins ",
            paste(round(100 * summary$oracle_frequency[c(2, 5)], 1), collapse = "/"), "%")
  }
}
write_csv(do.call(rbind, calibration_summaries), "output/tables/calibration_summary.csv")
write_csv(do.call(rbind, calibration_pairs), "output/tables/calibration_pairwise.csv")

capture.output(dput(config), file = "output/results/calibration_config.R")
