# Source from the repository root. No packages needed beyond R and renv.
load_project <- function() {
  for (file in c("config", "dgp", "analysis_menu", "oracle_loss", "diagnostics", "population_estimands"))
    source(file.path("R", paste0(file, ".R")), local = .GlobalEnv)
}

ensure_output <- function() {
  for (directory in c("figures", "tables", "results"))
    dir.create(file.path("output", directory), recursive = TRUE, showWarnings = FALSE)
}

run_replication <- function(r, p, theta_a, seed_base) {
  dat <- simulate_sample(p, seed_base + r)
  rows <- lapply(1:5, function(a) {
    fit <- fit_analysis(dat, a)
    loss <- oracle_loss(dat, fit, p)
    diagnostic <- design_diagnostics(dat, fit)
    data.frame(replication = r, seed = seed_base + r, analysis = a,
               status = fit$status, se_status = fit$se_status, T = fit$T, SE = fit$SE, t = fit$t,
               SE_ols = fit$SE_ols, theta = target_effect(p), theta_a = theta_a[a],
               b_a = theta_a[a] - target_effect(p), epsilon = fit$T - theta_a[a],
               as.list(loss), as.list(diagnostic), check.names = FALSE)
  })
  result <- do.call(rbind, rows)
  # Never silently drop failed candidates or pick from an incomplete menu.
  result$oracle_best <- NA_integer_
  if (all(result$status == "ok") && all(is.finite(result$L))) {
    # Ties within roundoff: smallest analysis index. Ties are recorded below.
    eligible <- which(result$L <= min(result$L) + 1e-12)
    result$oracle_best <- as.integer(result$analysis == eligible[1])
    result$oracle_tie_count <- length(eligible)
  } else result$oracle_tie_count <- NA_integer_
  result
}

run_montecarlo <- function(B, p, theta_a, seed_base, progress = FALSE) {
  records <- vector("list", B)
  for (r in seq_len(B)) {
    records[[r]] <- run_replication(r, p, theta_a, seed_base)
    if (progress && r %% 250L == 0L) message("Completed ", r, "/", B)
  }
  do.call(rbind, records)
}

summarize_pilot <- function(results) {
  do.call(rbind, lapply(split(results, results$analysis), function(d) {
    valid <- d$status == "ok"
    e2 <- (d$T[valid] - d$theta[valid])^2
    best <- d$oracle_best[!is.na(d$oracle_best)]
    frequency <- if (length(best)) mean(best) else NA_real_
    data.frame(analysis = d$analysis[1], label = analysis_labels[d$analysis[1]],
               theta_a = d$theta_a[1], b_a = d$b_a[1],
               mean_T = mean(d$T[valid]), sd_T = sd(d$T[valid]),
               mean_T_mcse = sd(d$T[valid]) / sqrt(sum(valid)),
               mean_epsilon = mean(d$epsilon[valid]), mean_SE = mean(d$SE[valid], na.rm = TRUE),
               mean_conditional_bias = mean(d$conditional_bias[valid]),
               mean_V = mean(d$V[valid]), mean_loss = mean(d$L[valid]),
               mean_loss_mcse = sd(d$L[valid]) / sqrt(sum(valid)),
               realized_mse = mean(e2), realized_mse_mcse = sd(e2) / sqrt(length(e2)),
               oracle_frequency = frequency,
               oracle_frequency_mcse = sqrt(frequency * (1 - frequency) / length(best)),
               failures = sum(!valid), se_failures = sum(d$se_status != "ok"), valid_oracle_designs = length(best))
  }))
}

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "NA")
