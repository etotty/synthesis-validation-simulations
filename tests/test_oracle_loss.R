data <- simulate_sample(config$dgp, config$seeds$test)
loss <- matrix(NA_real_, 4, 5)
for (a in 1:5) {
  fit <- fit_analysis(data, a)
  oracle <- oracle_loss(data, fit)
  loss[, a] <- oracle
  stopifnot(oracle["L"] >= 0, oracle["V"] > 0)
  # Independent normal-equation calculation for a well-conditioned test design.
  W <- fit$design$W
  c <- fit$design$contrast
  inverse <- solve(crossprod(W))
  direct <- config$dgp$sigma^2 * drop(t(c) %*% inverse %*% c)
  expect_close(oracle["V"], direct)
  noise_free <- data
  noise_free$Y <- data$m
  expect_close(oracle["conditional_mean_T"], fit_analysis(noise_free, a)$T)
  expect_close(oracle["L"], oracle["conditional_bias"]^2 + oracle["V"])
}
expect_close(loss[2, c(2, 3)], rep(0, 2))
stopifnot(loss[3, 5] <= loss[3, 2] + 1e-12, loss[4, 2] <= loss[4, 3] + 1e-12)
check <- check_oracle()
stopifnot(all(abs(check$simulated_variance - check$analytic_variance) < 6 * check$variance_mcse),
          all(abs(check$simulated_mean - check$expected) < 6 * check$mean_mcse),
          all(abs(check$simulated_mse - check$analytic_loss) < 6 * check$mse_mcse))
population <- population_estimands()
refined <- population_estimands(bound = 12, tolerance = 1e-12)
expect_close(population$theta_a, refined$theta_a, 1e-8)
expect_close(population$theta_a[c(2, 3)], rep(target_effect(config$dgp), 2), 1e-8)
# Remove all three sources of misspecification in a controlled test.
correct <- config$dgp
correct$beta[3] <- 0
correct$delta <- 0
correct$gamma <- 0
expect_close(population_estimands(correct)$theta_a, rep(target_effect(correct), 5), 1e-8)

# Gamma=0 recovers the original structural dominance on multiple fixed designs.
zero <- config$dgp
zero$gamma <- 0
pop_zero <- population_estimands(zero)
expect_close(pop_zero$theta_a[c(2, 3, 5)], rep(target_effect(zero), 3), 1e-8)
for (r in 0:4) {
  dat <- simulate_sample(zero, config$seeds$test + r)
  fits <- lapply(1:5, function(a) fit_analysis(dat, a))
  losses <- vapply(fits, function(f) oracle_loss(dat, f, zero), numeric(4))
  expect_close(losses[2, c(2, 3, 5)], rep(0, 3))
  stopifnot(losses[4, 5] <= losses[4, 2] + 1e-12,
            losses[4, 2] <= losses[4, 3] + 1e-12)
  # With the same design, gamma changes conditional bias but never variance.
  nonzero <- zero
  nonzero$gamma <- 0.10
  alt <- vapply(fits, function(f) oracle_loss(dat, f, nonzero), numeric(4))
  expect_close(alt[3, ], losses[3, ])
  expect_close(alt[2, c(2, 3)], rep(0, 2))
  stopifnot(abs(alt[2, 5]) > 1e-4)
  fit5 <- fits[[5]]
  rows <- fit5$design$rows
  omitted_term_effect <- nonzero$gamma * sum(fit5$design$h * dat$G[rows] * dat$X[rows])
  expect_close(alt[2, 5], omitted_term_effect)
}
nonzero <- zero
nonzero$gamma <- 0.10
pop_nonzero <- population_estimands(nonzero)
expect_close(pop_nonzero$theta_a[1:3], pop_zero$theta_a[1:3], 1e-8)
stopifnot(abs(pop_nonzero$b_a[5]) > 1e-4)
# Two-point integration response confirms the omitted slope's linear effect.
nonzero$gamma <- 0.20
expect_close(population_estimands(nonzero)$theta_a - pop_zero$theta_a,
             2 * (pop_nonzero$theta_a - pop_zero$theta_a), 1e-8)
