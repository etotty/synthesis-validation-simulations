# Deterministic numerical integration of E[WW'] and E[W m].
# Sum exactly over binary G,D; integrate over standard-normal X. No outcome
# or population Monte Carlo noise. Integrals use [-bound,bound]; a second
# run with a larger bound and tighter tolerance supplies a numerical check.
population_estimands <- function(p = config$dgp, trim = config$trim,
                                bound = 10, tolerance = 1e-10) {
  validate_dgp(p)
  results <- vector("list", 5)
  for (a in 1:5) {
    template <- regression_matrix(data.frame(X = 0, D = 0, G = 0), a)
    k <- ncol(template)
    A <- matrix(0, k, k)
    rhs <- numeric(k)
    max_error <- 0
    for (g in if (a <= 3) 1 else 0:1) {
      lower <- -bound
      upper <- bound
      if (a == 3) {
        intercept <- p$alpha[1] + p$alpha[3]
        slope <- p$alpha[2] + p$alpha[4]
        if (abs(slope) < 1e-14) {
          if (plogis(intercept) < trim[1] || plogis(intercept) > trim[2])
            stop("Empty trimmed population")
        } else {
          limits <- sort((qlogis(trim) - intercept) / slope)
          lower <- max(lower, limits[1])
          upper <- min(upper, limits[2])
        }
      }
      if (lower >= upper) stop("Empty integration interval")
      for (d in 0:1) {
        integrand <- function(x, j, column_index = NULL) {
          dat <- data.frame(X = x, G = rep(g, length(x)), D = rep(d, length(x)))
          W <- regression_matrix(dat, a)
          ps <- propensity(x, g, p)
          mass <- (if (g == 1) p$p_G else 1 - p$p_G) *
            (if (d == 1) ps else 1 - ps) * dnorm(x)
          other <- if (is.null(column_index)) conditional_mean(x, d, g, p) else W[, column_index]
          mass * W[, j] * other
        }
        for (j in seq_len(k)) {
          value <- integrate(integrand, lower, upper, j = j,
                             rel.tol = tolerance, abs.tol = tolerance, subdivisions = 500L)
          rhs[j] <- rhs[j] + value$value
          max_error <- max(max_error, value$abs.error)
          for (l in seq_len(k)) {
            value <- integrate(integrand, lower, upper, j = j, column_index = l,
                               rel.tol = tolerance, abs.tol = tolerance, subdivisions = 500L)
            A[j, l] <- A[j, l] + value$value
            max_error <- max(max_error, value$abs.error)
          }
        }
      }
    }
    coefficient <- solve(A, rhs)
    theta_a <- sum(effect_contrast(template, a) * coefficient)
    results[[a]] <- data.frame(analysis = a, label = analysis_labels[a],
                              theta = target_effect(p), theta_a = theta_a,
                              b_a = theta_a - target_effect(p),
                              max_integral_error = max_error, moment_condition = kappa(A, exact = TRUE))
  }
  do.call(rbind, results)
}
