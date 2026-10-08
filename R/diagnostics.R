safe_range_overlap <- function(ps, D) {
  if (!all(c(0, 1) %in% D)) return(NA_real_)
  max(0, min(max(ps[D == 0]), max(ps[D == 1])) -
        max(min(ps[D == 0]), min(ps[D == 1])))
}

design_diagnostics <- function(data, fit) {
  sub <- data[data$G == 1, , drop = FALSE]
  used <- data[fit$design$rows, , drop = FALSE]
  used_sub <- used[used$G == 1, , drop = FALSE]
  ps_quantile <- function(d, prob) if (nrow(d)) unname(quantile(d$ps, prob)) else NA_real_
  treatment_share <- function(d) if (nrow(d)) mean(d$D) else NA_real_
  # Conditioning diagnostic: 2-norm condition number of columns scaled to
  # unit Euclidean norm (no centering). Scale dependence remains documented.
  condition <- NA_real_
  leverage <- NA_real_
  if (fit$design$ok) {
    W <- fit$design$W
    scaled <- sweep(W, 2, sqrt(colSums(W^2)), "/")
    condition <- kappa(scaled, exact = TRUE)
    leverage <- max(fit$design$leverage)
  }
  c(n_total = nrow(data), n_subgroup = nrow(sub),
    n_subgroup_treated = sum(sub$D == 1), n_subgroup_control = sum(sub$D == 0),
    treatment_share_total = mean(data$D), treatment_share_subgroup = treatment_share(sub),
    n_used = nrow(used), n_subgroup_used = nrow(used_sub),
    n_subgroup_treated_used = sum(used_sub$D == 1), n_subgroup_control_used = sum(used_sub$D == 0),
    treatment_share_used = treatment_share(used),
    n_trimmed = if (fit$a == 3L) nrow(sub) - nrow(used) else 0,
    n_treated_trimmed = if (fit$a == 3L) sum(sub$D) - sum(used$D) else 0,
    n_control_trimmed = if (fit$a == 3L) sum(sub$D == 0) - sum(used$D == 0) else 0,
    ps_subgroup_min = ps_quantile(sub, 0), ps_subgroup_q05 = ps_quantile(sub, .05),
    ps_subgroup_median = ps_quantile(sub, .5), ps_subgroup_q95 = ps_quantile(sub, .95),
    ps_subgroup_max = ps_quantile(sub, 1),
    ps_used_min = ps_quantile(used, 0), ps_used_max = ps_quantile(used, 1),
    overlap_information_subgroup = if (nrow(sub)) mean(sub$ps * (1 - sub$ps)) else NA_real_,
    overlap_information_other = mean(data$ps[data$G == 0] * (1 - data$ps[data$G == 0])),
    subgroup_share_extreme_ps = if (nrow(sub)) mean(sub$ps < config$trim[1] | sub$ps > config$trim[2]) else NA_real_,
    ps_range_overlap_subgroup = safe_range_overlap(sub$ps, sub$D),
    max_leverage = leverage, scaled_condition_number = condition)
}
