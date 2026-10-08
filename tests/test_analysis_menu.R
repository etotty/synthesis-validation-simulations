data <- simulate_sample(config$dgp, config$seeds$test)
fits <- lapply(1:5, function(a) fit_analysis(data, a))
stopifnot(length(fits) == 5L)
for (a in 1:5) {
  f <- fits[[a]]
  stopifnot(f$status == "ok", all(c("T", "SE", "t", "coefficients", "design") %in% names(f)),
            is.finite(f$T), f$SE > 0)
  expect_close(f$t, f$T / f$SE)
  expect_close(drop(crossprod(f$design$W, f$design$h)), unname(f$design$contrast))
}
# Independent lm/formula reference for each specified analysis.
sub <- subset(data, G == 1)
trimmed <- subset(sub, ps >= config$trim[1] & ps <= config$trim[2])
references <- list(lm(Y ~ D + X, sub), lm(Y ~ D + X + I(X^2), sub),
                   lm(Y ~ D + X + I(X^2), trimmed),
                   lm(Y ~ D + X + I(X^2) + G, data),
                   lm(Y ~ D + D:G + X + I(X^2) + G, data))
for (a in 1:5) {
  reference <- coef(references[[a]])["D"]
  if (a == 5L) reference <- reference + coef(references[[a]])["D:G"]
  expect_close(fits[[a]]$T, reference)
}
expect_close(fits[[5]]$T, fits[[5]]$coefficients["D"] + fits[[5]]$coefficients["DG"])
stopifnot(length(fits[[3]]$design$rows) == nrow(trimmed))
# Inclusive threshold rule and no dependence on Y.
boundary <- data.frame(G = c(1, 1, 1, 1, 0), ps = c(.1, .9, .099, .901, .5))
stopifnot(identical(analysis_rows(boundary, 3), 1:2))
modified <- data
modified$Y <- -1000 * data$Y
stopifnot(identical(analysis_rows(data, 3), analysis_rows(modified, 3)))
# Explicit failure objects, including an empty subgroup, instead of silent drops.
no_treatment <- data
no_treatment$D <- 0
stopifnot(fit_analysis(no_treatment, 2)$status == "rank_deficient")
no_subgroup <- data
no_subgroup$G <- 0
stopifnot(fit_analysis(no_subgroup, 2)$status == "insufficient_rows")

# One treated observation can yield unit leverage and undefined HC3 even when
# the coefficient and oracle variance remain identified. Flag this explicitly.
one_treated <- data
one_treated$D <- 0
one_treated$D[which(one_treated$G == 1)[1]] <- 1
unit_leverage <- fit_analysis(one_treated, 2)
stopifnot(unit_leverage$status == "ok", unit_leverage$se_status == "leverage_one",
          is.na(unit_leverage$SE), is.finite(oracle_loss(one_treated, unit_leverage)["V"]))
