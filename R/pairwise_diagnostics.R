# One row per realized design. Pairwise validity does not require the other
# three candidates to fit; overall oracle eligibility remains a separate flag.
pairwise_designs <- function(results, tolerance = config$pairwise_tolerance) {
  two <- results[results$analysis == 2L, ]
  five <- results[results$analysis == 5L, ]
  stopifnot(!anyDuplicated(two$replication), !anyDuplicated(five$replication),
            setequal(two$replication, five$replication))
  five <- five[match(two$replication, five$replication), ]
  features <- c("replication", "seed", "n_subgroup", "n_subgroup_treated",
                "n_subgroup_control", "treatment_share_subgroup", "ps_subgroup_q05",
                "ps_subgroup_q95", "overlap_information_subgroup",
                "subgroup_share_extreme_ps", "ps_range_overlap_subgroup")
  out <- two[, features]
  out$ps_subgroup_spread <- two$ps_subgroup_q95 - two$ps_subgroup_q05
  for (name in c("L", "V", "conditional_bias", "max_leverage", "scaled_condition_number")) {
    out[[paste0(name, "_2")]] <- two[[name]]
    out[[paste0(name, "_5")]] <- five[[name]]
  }
  out$valid <- two$status == "ok" & five$status == "ok" & is.finite(two$L) & is.finite(five$L)
  out$delta_L25 <- out$L_2 - out$L_5
  out$absolute_gap <- abs(out$delta_L25)
  out$variance_gain <- out$V_2 - out$V_5
  out$bias_penalty <- out$conditional_bias_5^2 - out$conditional_bias_2^2
  # Excess loss as a fraction of the better of the two candidates' loss.
  out$relative_gap <- out$absolute_gap / pmin(out$L_2, out$L_5)
  out$region <- ifelse(!out$valid, "invalid", ifelse(out$delta_L25 < -tolerance,
    "2_wins", ifelse(out$delta_L25 > tolerance, "5_wins", "tie")))
  out$overall_winner <- NA_integer_
  winner <- results[!is.na(results$oracle_best) & results$oracle_best == 1L, ]
  out$overall_winner <- winner$analysis[match(out$replication, winner$replication)]
  out
}

summarize_pairwise <- function(pair, tolerance = config$pairwise_tolerance) {
  d <- pair[pair$valid, ]
  if (!nrow(d)) stop("No valid 2-versus-5 designs")
  fraction <- function(x) mean(x)
  p2 <- fraction(d$region == "2_wins")
  p5 <- fraction(d$region == "5_wins")
  data.frame(n_designs = nrow(pair), n_valid = nrow(d), n_invalid = sum(!pair$valid),
    tie_tolerance = tolerance, mean_delta = mean(d$delta_L25), median_delta = median(d$delta_L25),
    mean_delta_mcse = sd(d$delta_L25) / sqrt(nrow(d)),
    fraction_2_strict = fraction(d$delta_L25 < 0), fraction_5_strict = fraction(d$delta_L25 > 0),
    fraction_2_wins = p2, fraction_5_wins = p5, fraction_tied = fraction(d$region == "tie"),
    fraction_2_mcse = sqrt(p2 * (1 - p2) / nrow(d)), fraction_5_mcse = sqrt(p5 * (1 - p5) / nrow(d)),
    mean_absolute_gap = mean(d$absolute_gap), median_absolute_gap = median(d$absolute_gap),
    mean_relative_gap = mean(d$relative_gap), median_relative_gap = median(d$relative_gap),
    fraction_gap_below_1e6 = fraction(d$absolute_gap < 1e-6),
    fraction_gap_below_1e4 = fraction(d$absolute_gap < 1e-4),
    fraction_gap_below_1e3 = fraction(d$absolute_gap < 1e-3))
}

pairwise_quantiles <- function(pair) {
  d <- pair[pair$valid, ]
  probabilities <- c(0, .01, .05, .25, .5, .75, .95, .99, 1)
  data.frame(probability = probabilities,
    delta_L25 = unname(quantile(d$delta_L25, probabilities)),
    absolute_gap = unname(quantile(d$absolute_gap, probabilities)),
    relative_gap = unname(quantile(d$relative_gap, probabilities)))
}

pairwise_regions <- function(pair) {
  features <- c("n_subgroup", "n_subgroup_treated", "n_subgroup_control",
    "treatment_share_subgroup", "ps_subgroup_q05", "ps_subgroup_q95", "ps_subgroup_spread",
    "overlap_information_subgroup", "subgroup_share_extreme_ps", "ps_range_overlap_subgroup",
    "V_2", "V_5", "conditional_bias_2", "conditional_bias_5", "variance_gain", "bias_penalty",
    "max_leverage_2", "max_leverage_5", "scaled_condition_number_2", "scaled_condition_number_5",
    "delta_L25", "absolute_gap", "relative_gap")
  do.call(rbind, lapply(c("2_wins", "5_wins", "tie"), function(region) {
    d <- pair[pair$region == region, ]
    do.call(rbind, lapply(features, function(feature) {
      x <- d[[feature]]
      data.frame(region = region, n_designs = nrow(d), feature = feature,
                 mean = if (length(x)) mean(x) else NA_real_,
                 median = if (length(x)) median(x) else NA_real_)
    }))
  }))
}
