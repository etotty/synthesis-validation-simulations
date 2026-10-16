# Hand-computable crossing, equal-mean, dominance, and zero-gain examples.
g <- oracle_adaptation_gain(c(1, 3), c(2, 2))
expect_close(c(g$mean_min_loss, g$oracle_gain, g$proportional_gain), c(1.5, .5, .25))
stopifnot(!g$fixed_choice_separated, g$gain_lower <= .5, g$gain_upper >= .5)
for (loss in list(list(c(1, 1), c(2, 3)), list(c(2, 3), c(1, 1)),
                 list(c(2, 2), c(2, 2)), list(c(1, 3, 2), c(4, 1, 2)))) {
  g <- oracle_adaptation_gain(loss[[1]], loss[[2]])
  expect_close(g$oracle_gain, min(mean(loss[[1]]), mean(loss[[2]])) - mean(pmin(loss[[1]], loss[[2]])))
  expect_close(g$proportional_gain, g$oracle_gain / min(g$mean_L2, g$mean_L5))
  stopifnot(g$oracle_gain >= 0, g$mean_min_loss <= min(g$mean_L2, g$mean_L5)+1e-12,
            g$gain_lower <= g$oracle_gain, g$gain_upper >= g$oracle_gain)
  swapped <- oracle_adaptation_gain(loss[[2]], loss[[1]])
  expect_close(g$oracle_gain, swapped$oracle_gain)
  expect_close(g$proportional_lower, swapped$proportional_lower)
}
stopifnot(inherits(try(oracle_adaptation_gain(c(1, NA), c(2, 3)), silent=TRUE), "try-error"))
# Diagnostic CRNs are reproducible and nested in n and p_G without modifying
# the pilot sampler. Thresholding independent uniforms implements Bernoulli p.
p <- config$dgp
original <- p
z <- design_value_draws(config$seeds$test+91L)
stopifnot(identical(z, design_value_draws(config$seeds$test+91L)))
a <- design_value_sample(p, z)
p$n <- 1000L
b <- design_value_sample(p, z)
stopifnot(identical(a[1:1000, ], b))
p$p_G <- .08
c <- design_value_sample(p, z)
stopifnot(all(c$G <= b$G), identical(c$X, b$X), identical(p$gamma, original$gamma))
expect_close(c$m, conditional_mean(c$X, c$D, c$G, p))
expect_close(c$Y-c$m, p$sigma*z$u[1:1000])
stopifnot(identical(config$dgp, original))
# Saved diagnostic summaries must equal the existing oracle on supplied data.
pop <- population_estimands(p)
x <- lapply(1:3, function(r) run_replication(r, p, pop$theta_a, config$seeds$test,
  dat=design_value_sample(p, design_value_draws(config$seeds$test+r))))
pair <- pairwise_designs(do.call(rbind, x))
g <- oracle_adaptation_gain(pair$L_2, pair$L_5)
stopifnot(all(pair$valid), g$oracle_gain >= 0)
# A new R process reproduces the full supplied-data path, independently of
# workspace state, including generated losses, diagnostics and summaries.
script <- tempfile(fileext=".R")
result <- tempfile(fileext=".rds")
writeLines(c('source("R/simulation_helpers.R"); load_project()',
  'p <- config$dgp; p$n <- 1000L; p$p_G <- .08; pop <- population_estimands(p)',
  'x <- lapply(1:3, function(r) run_replication(r,p,pop$theta_a,config$seeds$test,dat=design_value_sample(p,design_value_draws(config$seeds$test+r))))',
  sprintf('saveRDS(x, %s)', encodeString(result, quote='"'))), script)
status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)))
stopifnot(status==0L, identical(x, readRDS(result)))
unlink(c(script, result))
