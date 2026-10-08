source("R/simulation_helpers.R")
load_project()
ensure_output()
# Three documented calibration cases; shared seeds isolate the effect of
# reducing Delta while retaining the common target theta = 0.50.
cases <- list(initial = config$dgp, intermediate = config$dgp, final = config$dgp)
cases$initial$tau0 <- 0.25
cases$initial$delta <- 0.25
cases$intermediate$tau0 <- 0.32
cases$intermediate$delta <- 0.18
calibration_summaries <- list()
for (case in names(cases)) {
  parameters <- cases[[case]]
  population <- population_estimands(parameters)
  calibration <- run_montecarlo(config$calibration_reps, parameters, population$theta_a,
                               config$seeds$calibration)
  summary <- summarize_pilot(calibration)
  summary$calibration_case <- case
  summary$tau0 <- parameters$tau0
  summary$delta <- parameters$delta
  calibration_summaries[[case]] <- summary
}
summary <- do.call(rbind, calibration_summaries)
write_csv(summary, "output/tables/calibration_summary.csv")
print(summary[, c("calibration_case", "analysis", "theta_a", "mean_loss", "oracle_frequency")],
      row.names = FALSE)
