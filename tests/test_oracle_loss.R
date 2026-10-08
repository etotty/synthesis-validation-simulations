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
expect_close(loss[2, c(2, 3, 5)], rep(0, 3))
stopifnot(loss[4, 5] <= loss[4, 2] + 1e-12, loss[4, 2] <= loss[4, 3] + 1e-12)
check <- check_oracle()
stopifnot(all(abs(check$simulated_variance - check$analytic_variance) < 6 * check$variance_mcse),
          all(abs(check$simulated_mean - check$expected) < 6 * check$mean_mcse),
          all(abs(check$simulated_mse - check$analytic_loss) < 6 * check$mse_mcse))
population <- population_estimands()
refined <- population_estimands(bound = 12, tolerance = 1e-12)
expect_close(population$theta_a, refined$theta_a, 1e-8)
expect_close(population$theta_a[c(2, 3, 5)], rep(target_effect(config$dgp), 3), 1e-8)
# Remove only the two sources of misspecification in a controlled test.
correct <- config$dgp
correct$beta[3] <- 0
correct$delta <- 0
expect_close(population_estimands(correct)$theta_a, rep(target_effect(correct), 5), 1e-8)
