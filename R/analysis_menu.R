analysis_labels <- c("Subgroup linear", "Subgroup flexible",
                     "Subgroup flexible trimmed", "Full constant effect",
                     "Full interaction")

# Known-DGP propensity trimming, inclusive endpoints, independent of Y.
# This is an oracle overlap benchmark, not an estimated-propensity procedure.
analysis_rows <- function(data, a, trim = config$trim) {
  stopifnot(a %in% 1:5, length(trim) == 2L,
            0 < trim[1], trim[1] < trim[2], trim[2] < 1)
  if (a == 3L) return(which(data$G == 1 & data$ps >= trim[1] & data$ps <= trim[2]))
  if (a <= 2L) return(which(data$G == 1))
  seq_len(nrow(data))
}

regression_matrix <- function(data, a) {
  W <- cbind(intercept = rep(1, nrow(data)), D = data$D, X = data$X)
  if (a != 1L) W <- cbind(W, X2 = data$X^2)
  if (a >= 4L) W <- cbind(W, G = data$G)
  if (a == 5L) W <- cbind(W, DG = data$D * data$G)
  W
}

effect_contrast <- function(W, a) {
  contrast <- setNames(rep(0, ncol(W)), colnames(W))
  contrast["D"] <- 1
  if (a == 5L) contrast["DG"] <- 1
  contrast
}

# Every estimator is h'Y. QR avoids explicitly inverting W'W.
prepare_analysis <- function(data, a, trim = config$trim) {
  rows <- analysis_rows(data, a, trim)
  W <- regression_matrix(data[rows, , drop = FALSE], a)
  contrast <- effect_contrast(W, a)
  if (nrow(W) <= ncol(W)) return(list(ok = FALSE, reason = "insufficient_rows", rows = rows))
  decomposition <- qr(W, tol = 1e-10)
  if (decomposition$rank < ncol(W)) return(list(ok = FALSE, reason = "rank_deficient", rows = rows))
  Q <- qr.Q(decomposition)
  R <- qr.R(decomposition)
  h <- drop(Q %*% backsolve(R, contrast[decomposition$pivot], transpose = TRUE))
  list(ok = TRUE, reason = "ok", rows = rows, W = W,
       contrast = contrast, qr = decomposition, h = h,
       leverage = rowSums(Q^2))
}

fit_analysis <- function(data, a, trim = config$trim) {
  design <- prepare_analysis(data, a, trim)
  if (!design$ok) return(list(a = a, status = design$reason, T = NA_real_,
                             SE = NA_real_, t = NA_real_, SE_ols = NA_real_,
                             se_status = "fit_failed", design = design))
  y <- data$Y[design$rows]
  coefficients <- qr.coef(design$qr, y)
  residuals <- qr.resid(design$qr, y)
  estimate <- sum(design$h * y)
  # HC3 reported SE: random-design robust inference for misspecified candidates.
  # Distinct from the known-sigma fixed-design oracle variance.
  se_status <- if (any(1 - design$leverage <= 1e-12)) "leverage_one" else "ok"
  se <- if (se_status == "ok") {
    sqrt(sum((design$h * residuals / (1 - design$leverage))^2))
  } else NA_real_
  se_ols <- sqrt(sum(residuals^2) / (length(y) - ncol(design$W)) * sum(design$h^2))
  list(a = a, status = "ok", T = estimate, SE = se,
       t = estimate / se, SE_ols = se_ols, se_status = se_status,
       coefficients = coefficients, design = design)
}
