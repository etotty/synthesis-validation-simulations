validate_dgp <- function(p) {
  stopifnot(length(p$n) == 1L, p$n >= 10, p$n == as.integer(p$n),
            p$p_G > 0, p$p_G < 1, length(p$alpha) == 4L,
            length(p$beta) == 4L, p$sigma > 0,
            all(is.finite(unlist(p))))
  invisible(TRUE)
}

target_effect <- function(p) p$tau0 + p$delta

propensity <- function(X, G, p) {
  a <- p$alpha
  plogis(a[1] + a[2] * X + a[3] * G + a[4] * X * G)
}

conditional_mean <- function(X, D, G, p) {
  b <- p$beta
  b[1] + b[2] * X + b[3] * X^2 + b[4] * G +
    p$tau0 * D + p$delta * D * G
}

simulate_sample <- function(p = config$dgp, seed) {
  validate_dgp(p)
  do.call(RNGkind, as.list(config$rng))
  set.seed(seed)
  G <- rbinom(p$n, 1, p$p_G)
  X <- rnorm(p$n)
  ps <- propensity(X, G, p)
  D <- rbinom(p$n, 1, ps)
  m <- conditional_mean(X, D, G, p)
  Y <- m + rnorm(p$n, sd = p$sigma)
  data.frame(X = X, G = G, D = D, ps = ps, m = m, Y = Y)
}
