oracle_loss <- function(data, fit, p = config$dgp) {
  if (!fit$design$ok) return(c(conditional_mean_T = NA_real_,
                              conditional_bias = NA_real_, V = NA_real_, L = NA_real_))
  design <- fit$design
  # The conditional mean is recalculated from X,D,G and the known DGP.
  rows <- design$rows
  m <- conditional_mean(data$X[rows], data$D[rows], data$G[rows], p)
  expected_T <- sum(design$h * m)
  bias <- expected_T - target_effect(p)
  variance <- p$sigma^2 * sum(design$h^2)
  c(conditional_mean_T = expected_T, conditional_bias = bias,
    V = variance, L = bias^2 + variance)
}

# Controlled repeated-noise experiment, holding X,D,G (and trimming) fixed.
# This checks ALL five estimators, including the misspecified ones.
check_oracle <- function(p = config$dgp, B = config$oracle_check_reps) {
  data <- simulate_sample(p, config$seeds$oracle_check)
  fits <- lapply(1:5, function(a) fit_analysis(data, a))
  stopifnot(all(vapply(fits, function(f) f$design$ok, logical(1))))
  H <- matrix(0, nrow(data), 5)
  for (a in 1:5) H[fits[[a]]$design$rows, a] <- fits[[a]]$design$h
  truth <- vapply(fits, function(f) oracle_loss(data, f, p), numeric(4))
  set.seed(config$seeds$oracle_check + 1L)
  estimates <- matrix(NA_real_, B, 5)
  # Batching avoids an n x B allocation, with an invariant RNG sequence.
  for (start in seq.int(1L, B, by = 250L)) {
    idx <- start:min(start + 249L, B)
    U <- matrix(rnorm(nrow(data) * length(idx), sd = p$sigma), nrow(data))
    estimates[idx, ] <- sweep(crossprod(U, H), 2, truth["conditional_mean_T", ], "+")
  }
  mse <- colMeans((estimates - target_effect(p))^2)
  data.frame(analysis = 1:5, expected = truth[1, ], simulated_mean = colMeans(estimates),
             analytic_variance = truth[3, ], simulated_variance = apply(estimates, 2, var),
             analytic_loss = truth[4, ], simulated_mse = mse,
             variance_mcse = truth[3, ] * sqrt(2 / (B - 1)),
             mean_mcse = sqrt(truth[3, ] / B),
             mse_mcse = sqrt((2 * truth[3, ]^2 + 4 * truth[2, ]^2 * truth[3, ]) / B))
}
