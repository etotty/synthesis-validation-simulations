# Run from the repository root. Optional argument: grid, confirmation, or report;
# default all runs both stages and regenerates the combined tracked report.
source("R/simulation_helpers.R")
load_project()
ensure_output()
args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args)) args[1] else "all"
stopifnot(mode %in% c("all", "grid", "confirmation", "report"))
original_dgp <- config$dgp
for (stage in c("grid", "confirmation")) {
  if (mode != "all" && mode != stage) next
  cells <- if (stage == "grid") expand.grid(n=config$design_value$n, p_G=config$design_value$p_G) else config$design_value$confirmation_cells
  if (!nrow(cells)) next
  B <- if (stage == "grid") config$design_value$grid_reps else config$design_value$confirmation_reps
  seed <- if (stage == "grid") config$seeds$design_grid else config$seeds$design_confirmation
  runs <- lapply(seq_len(nrow(cells)), function(i)
    run_design_value_cell(cells$n[i], cells$p_G[i], B, seed, stage))
  dput(design_value_stage_config(stage), file=sprintf("output/results/design_value_%s_config.R", stage))
  for (name in c("summary", "quantiles", "regions", "failures")) {
    values <- do.call(rbind, lapply(runs, `[[`, name))
    write_csv(values, sprintf("output/tables/design_value_%s_%s.csv", stage, name))
  }
}
for (stage in c("grid", "confirmation")) {
  path <- sprintf("output/results/design_value_%s_config.R", stage)
  if (file.exists(path) && !identical(dget(path), design_value_stage_config(stage)))
    stop("Stale ", stage, " outputs: rerun with mode all")
}
# Persist exact diagnostic configuration and environment, separately from pilot.
capture.output(dput(list(dgp=config$dgp, diagnostic=config$design_value,
  seeds=config$seeds[c("design_grid", "design_confirmation")], rng=config$rng,
  trim=config$trim, pairwise_tolerance=config$pairwise_tolerance,
  oracle_tolerance=config$oracle_tolerance)), file="output/results/design_value_config.R")
capture.output(sessionInfo(), file="output/results/design_value_session_info.txt")
for (name in c("summary", "quantiles", "regions", "failures")) {
  paths <- sprintf("output/tables/design_value_%s_%s.csv", c("grid", "confirmation"), name)
  combined <- do.call(rbind, lapply(paths[file.exists(paths)], read.csv))
  if (name=="summary") {
    combined$zero_observed_gain <- combined$oracle_gain == 0
    combined$rare_win_probability_upper <- ifelse(combined$zero_observed_gain,
      1-.05^(1/combined$n_valid), NA_real_)
    # A zero observed excess-loss distribution gives a zero local MCSE but
    # cannot supply an informative upper bound on rare unobserved losses.
    combined$gain_upper[combined$zero_observed_gain] <- NA_real_
    combined$proportional_upper[combined$zero_observed_gain] <- NA_real_
  }
  assign(name, combined)
  write_csv(combined, sprintf("output/tables/design_value_%s.csv", if (name=="summary") "grid" else name))
}
stopifnot(identical(config$dgp, original_dgp))
write_design_value_report(summary, failures)
plot_design_value(summary)
message("Design-value diagnostic generated; canonical DGP unchanged.")
