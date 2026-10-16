first <- simulate_sample(config$dgp, config$seeds$test)
second <- simulate_sample(config$dgp, config$seeds$test)
stopifnot(identical(first, second), nrow(first) == config$dgp$n,
          all(first$G %in% 0:1), all(first$D %in% 0:1),
          all(first$ps > 0 & first$ps < 1),
          !identical(first, simulate_sample(config$dgp, config$seeds$test + 1L)))
expect_close(first$m, conditional_mean(first$X, first$D, first$G, config$dgp))
expect_close(conditional_mean(first$X, 1, first$G, config$dgp) -
             conditional_mean(first$X, 0, first$G, config$dgp),
             config$dgp$tau0 + config$dgp$delta * first$G)

# Direct formula check, independent of conditional_mean().
p <- config$dgp
expected <- p$beta[1] + p$beta[2] * first$X + p$beta[3] * first$X^2 +
  p$beta[4] * first$G + p$gamma * first$G * first$X +
  p$tau0 * first$D + p$delta * first$D * first$G
expect_close(first$m, expected)
zero <- p
zero$gamma <- 0
baseline <- simulate_sample(zero, config$seeds$test)
stopifnot(identical(first[c("X", "D", "G", "ps")], baseline[c("X", "D", "G", "ps")]))
expect_close(first$m - baseline$m, p$gamma * first$G * first$X)
expect_close(first$Y - baseline$Y, p$gamma * first$G * first$X)
# A larger nonzero slope increment gives the same checks independently of the
# final calibration. No D*X treatment heterogeneity is introduced.
p$gamma <- 0.20
alternative <- simulate_sample(p, config$seeds$test)
expect_close(alternative$m - baseline$m, p$gamma * first$G * first$X)
expect_close(conditional_mean(first$X, 1, first$G, p) -
             conditional_mean(first$X, 0, first$G, p),
             p$tau0 + p$delta * first$G)
