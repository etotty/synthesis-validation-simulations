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
# Dominance is a mathematical implication, not a calibration target.
L <- matrix(x$L, ncol = 5, byrow = TRUE)
stopifnot(all(L[, 5] <= L[, 2] + 1e-12), all(L[, 2] <= L[, 3] + 1e-12))
