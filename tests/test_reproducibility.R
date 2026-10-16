population <- population_estimands()
x <- run_montecarlo(5L, config$dgp, population$theta_a, config$seeds$test)
y <- run_montecarlo(5L, config$dgp, population$theta_a, config$seeds$test)
stopifnot(identical(x, y), all(x$epsilon == x$T - x$theta_a),
          all(x$L >= 0), all(x$status == "ok"))
# Replications are invariant to execution order.
for (r in c(5L, 2L, 1L)) {
  independent <- run_replication(r, config$dgp, population$theta_a, config$seeds$test)
  original <- x[x$replication == r, ]
  rownames(original) <- NULL
  stopifnot(identical(original, independent))
}
# Fresh R process without workspace/startup state or attached packages.
script <- tempfile(fileext = ".R")
result <- tempfile(fileext = ".rds")
writeLines(c('source("R/simulation_helpers.R"); load_project()',
             'p <- population_estimands()',
             sprintf('saveRDS(run_montecarlo(5L, config$dgp, p$theta_a, config$seeds$test), %s)',
                     encodeString(result, quote = '"'))), script)
status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)))
stopifnot(status == 0L, identical(x, readRDS(result)))
unlink(c(script, result))
# Only trimming dominance applies at nonzero gamma.
L <- matrix(x$L, ncol = 5, byrow = TRUE)
stopifnot(all(L[, 2] <= L[, 3] + 1e-12))

# Pairwise results retain exact replication alignment and decomposition.
pair <- pairwise_designs(x)
stopifnot(identical(pair, pairwise_designs(y)))
expect_close(pair$delta_L25, L[, 2] - L[, 5])
expect_close(pair$delta_L25, pair$variance_gain - pair$bias_penalty)
summary <- summarize_pairwise(pair)
expect_close(summary$fraction_2_wins + summary$fraction_5_wins + summary$fraction_tied, 1)
# Explicit tolerance-edge cases: disjoint, exhaustive classification including
# near-ties of both signs and an invalid comparison.
fixture <- x
five <- which(fixture$analysis == 5)
two <- which(fixture$analysis == 2)
fixture$L[five] <- fixture$L[two] - c(-1e-3, 1e-3, -5e-11, 5e-11, 0)
fixture$status[five[5]] <- "rank_deficient"
regions <- pairwise_designs(fixture)$region
stopifnot(identical(regions, c("2_wins", "5_wins", "tie", "tie", "invalid")))
