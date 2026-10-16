# v0.1.2 diagnostic: independent replicates, coupled cells within a replicate.
# Prefixes keep X and uniforms identical across n; thresholding the same U_G
# nests subgroup membership across p_G. The original pilot RNG is unchanged.
design_value_draws <- function(seed, max_n = max(config$design_value$n)) {
  do.call(RNGkind, as.list(config$rng))
  set.seed(seed)
  list(U_G = runif(max_n), X = rnorm(max_n), U_D = runif(max_n), u = rnorm(max_n))
}

design_value_sample <- function(p, draws) {
  validate_dgp(p)
  stopifnot(all(lengths(draws) >= p$n))
  i <- seq_len(p$n)
  X <- draws$X[i]
  G <- as.integer(draws$U_G[i] < p$p_G)
  ps <- propensity(X, G, p)
  D <- as.integer(draws$U_D[i] < ps)
  m <- conditional_mean(X, D, G, p)
  data.frame(X, G, D, ps, m, Y = m + p$sigma * draws$u[i])
}

# G = min(E[A2], E[A5]), where Ak = Lk - min(L2,L5).
# Estimate the two paired means before taking their minimum. The same identity
# for g uses ratios E[Ak]/E[Lk]. Simultaneous intervals for both branches avoid
# treating the empirically best fixed choice as known near a mean-loss tie.
oracle_adaptation_gain <- function(L2, L5, level = .95) {
  stopifnot(length(L2) == length(L5), length(L2) >= 2L,
            all(is.finite(c(L2, L5))), all(c(L2, L5) > 0), level > 0, level < 1)
  L <- cbind(L2, L5)
  M <- pmin(L2, L5)
  A <- L - M
  N <- length(M)
  means <- colMeans(L)
  gains <- colMeans(A)
  ratios <- gains / means
  gain_se <- apply(A, 2, sd) / sqrt(N)
  ratio_se <- vapply(1:2, function(k)
    sd((A[, k] - ratios[k] * L[, k]) / means[k]) / sqrt(N), numeric(1))
  # Bonferroni two-sided normal intervals for two branches: joint >= level
  # asymptotically, without differentiating the nonsmooth minimum.
  z <- qnorm(1 - (1 - level) / 4)
  best <- which.min(means)
  data.frame(mean_L2 = means[1], mean_L5 = means[2], mean_min_loss = mean(M),
    best_fixed = c(2L, 5L)[best], oracle_gain = min(gains),
    gain_mcse_local = gain_se[best],
    gain_lower = max(0, min(gains - z * gain_se)),
    gain_upper = min(gains + z * gain_se),
    proportional_gain = min(ratios), proportional_mcse_local = ratio_se[best],
    proportional_lower = max(0, min(ratios - z * ratio_se)),
    proportional_upper = min(1, min(ratios + z * ratio_se)),
    mean_L2_mcse = sd(L2) / sqrt(N), mean_L5_mcse = sd(L5) / sqrt(N),
    mean_min_mcse = sd(M) / sqrt(N),
    fixed_gap_mcse = sd(L2 - L5) / sqrt(N),
    fixed_choice_separated = abs(diff(means)) > z * sd(L2 - L5) / sqrt(N))
}

summarize_design_value <- function(pair, population, failures) {
  # Never silently condition on successful designs for the primary benchmark.
  # A cell with a failed pair requires scientific review; stop rather than
  # assign arbitrary losses or report a selected-sample gain.
  if (!all(pair$valid)) stop("Invalid pairwise design: inspect saved failures before reporting gain")
  gain <- oracle_adaptation_gain(pair$L_2, pair$L_5)
  ps <- summarize_pairwise(pair)
  q <- function(x) setNames(as.list(unname(quantile(x, c(0, .05, .5, .95, 1)))),
                            c("min", "q05", "median", "q95", "max"))
  distribution <- function(x, prefix) {
    values <- c(list(mean = mean(x)), q(x))
    setNames(values, paste0(prefix, "_", names(values)))
  }
  data.frame(gain, ps,
    delta_q05 = unname(quantile(pair$delta_L25, .05)),
    delta_q25 = unname(quantile(pair$delta_L25, .25)),
    delta_q75 = unname(quantile(pair$delta_L25, .75)),
    delta_q95 = unname(quantile(pair$delta_L25, .95)),
    theta_5 = population$theta_a[5], b_5 = population$b_a[5],
    mean_V2 = mean(pair$V_2), mean_V5 = mean(pair$V_5),
    mean_bias5_squared = mean(pair$conditional_bias_5^2),
    distribution(pair$n_subgroup, "subgroup"),
    distribution(pair$n_subgroup_treated, "treated"),
    distribution(pair$n_subgroup_control, "control"),
    fit_failures = sum(failures$fit_failures),
    rank_failures = sum(failures$rank_failures),
    se_failures = sum(failures$se_failures), row.names = NULL)
}

run_design_value_cell <- function(n, p_G, B, seed_base, stage) {
  p <- config$dgp
  p$n <- as.integer(n)
  p$p_G <- p_G
  population <- population_estimands(p)
  pair <- vector("list", B)
  failures <- data.frame(analysis = 1:5, fit_failures = 0L, rank_failures = 0L,
                         se_failures = 0L, oracle_wins = 0L)
  trimming_violations <- 0L
  for (r in seq_len(B)) {
    dat <- design_value_sample(p, design_value_draws(seed_base + r))
    result <- run_replication(r, p, population$theta_a, seed_base, dat = dat)
    pair[[r]] <- pairwise_designs(result)
    failures$fit_failures <- failures$fit_failures + (result$status != "ok")
    failures$rank_failures <- failures$rank_failures + (result$status == "rank_deficient")
    failures$se_failures <- failures$se_failures + (result$se_status != "ok")
    failures$oracle_wins <- failures$oracle_wins + ifelse(is.na(result$oracle_best), 0, result$oracle_best)
    if (all(is.finite(result$L[2:3])))
      trimming_violations <- trimming_violations + (result$L[2] > result$L[3] + 1e-12)
    if (r %% 1000L == 0L) message(stage, ": n=", n, ", p_G=", p_G, ": ", r, "/", B)
  }
  pair <- do.call(rbind, pair)
  id <- data.frame(stage = stage, n = n, p_G = p_G)
  prefix <- sprintf("design_value_%s_n%d_p%02d", stage, n, round(100 * p_G))
  write_csv(pair, file.path("output/results", paste0(prefix, ".csv")))
  write_csv(failures, file.path("output/tables", paste0(prefix, "_failures.csv")))
  out <- summarize_design_value(pair, population, failures)
  out$trimming_violations <- trimming_violations
  list(summary = cbind(id, out),
       quantiles = cbind(id, pairwise_quantiles(pair)),
       regions = cbind(id, pairwise_regions(pair)), failures = cbind(id, failures))
}

# Code-generated Markdown tables keep every reported number tied to a run.
design_value_table <- function(x, digits = 6L) {
  formatted <- lapply(x, function(v) {
    if (!is.numeric(v) || is.integer(v)) return(as.character(v))
    out <- formatC(v, format="f", digits=digits)
    small <- is.finite(v) & v != 0 & abs(v) < 10^(-digits)
    out[small] <- formatC(v[small], format="g", digits=3)
    out
  })
  rows <- do.call(cbind, formatted)
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(rows, 1, function(v) paste0("| ", paste(v, collapse = " | "), " |")))
}

write_design_value_report <- function(summary, failures) {
  s <- summary
  main <- data.frame(stage = s$stage, n = s$n, p_G = s$p_G, B = s$n_designs,
    L2 = s$mean_L2, L5 = s$mean_L5, oracle = s$mean_min_loss, G = s$oracle_gain,
    G_low = s$gain_lower, G_high = s$gain_upper,
    g_percent = 100*s$proportional_gain, g_low = 100*s$proportional_lower,
    g_high = 100*s$proportional_upper, fixed = s$best_fixed)
  gaps <- data.frame(stage = s$stage, n = s$n, p_G = s$p_G,
    win2_percent = 100*s$fraction_2_wins, win5_percent = 100*s$fraction_5_wins,
    tie_percent = 100*s$fraction_tied, mean = s$mean_delta, median = s$median_delta,
    q05 = s$delta_q05, q25 = s$delta_q25, q75 = s$delta_q75, q95 = s$delta_q95,
    mean_absolute = s$mean_absolute_gap, mean_relative_percent = 100*s$mean_relative_gap)
  components <- data.frame(stage = s$stage, n = s$n, p_G = s$p_G,
    theta5 = s$theta_5, b5 = s$b_5, V2 = s$mean_V2, V5 = s$mean_V5,
    bias5_squared = s$mean_bias5_squared)
  counts <- data.frame(stage = s$stage, n = s$n, p_G = s$p_G,
    subgroup_mean = s$subgroup_mean, subgroup_q05 = s$subgroup_q05,
    subgroup_median = s$subgroup_median, subgroup_q95 = s$subgroup_q95,
    treated_mean = s$treated_mean, treated_q05 = s$treated_q05,
    control_mean = s$control_mean, control_q05 = s$control_q05,
    min_treated = s$treated_min, min_control = s$control_min,
    fit_failures = s$fit_failures, rank_failures = s$rank_failures)
  band <- ifelse(s$proportional_gain < .01, "below 1%",
    ifelse(s$proportional_gain < .03, "1-3%", ifelse(s$proportional_gain < .05, "3-5%", "above 5%")))
  overall <- do.call(rbind, lapply(seq_len(nrow(s)), function(i) {
    f <- failures[failures$stage==s$stage[i] & failures$n==s$n[i] & failures$p_G==s$p_G[i], ]
    f <- f[order(f$analysis), ]
    data.frame(stage=s$stage[i], n=s$n[i], p_G=s$p_G[i],
      complete_designs=sum(f$oracle_wins),
      shares_1_to_5=paste(sprintf("%.2f", 100*f$oracle_wins/sum(f$oracle_wins)), collapse=" / "))
  }))
  lines <- c("# v0.1.2: value of realized-design adaptation", "",
    "Generated by `Rscript scripts/03_calibrate_design_value.R`. The canonical v0.1.1 baseline is unchanged: n=1500, p_G=0.12. No later information or validation stages are implemented.", "",
    "## Design and interpretation", "",
    "Grid: n in {1000,1500,2000}, p_G in {0.08,0.12,0.16}. Fixed: alpha=(-0.3,0.4,0.9,1.0), beta=(0,0.7,0.5,0.3), gamma=0.10, tau0=0.25, Delta=0.25, sigma=1, theta=0.50. Five analyses, known-propensity trimming [0.10,0.90], and exact conditional MSE are unchanged.", "",
    "G = min(E[L2],E[L5]) - E[min(L2,L5)]; g = G/min(E[L2],E[L5]). This is the maximum average loss improvement from perfectly observing the realized pairwise ranking, relative to the best fixed choice. Winner shares and average absolute gaps are different quantities: only gains on designs where that fixed choice loses contribute to G.", "",
    sprintf("First pass: %d replications per cell, seed %d+r. Confirmation: %d independent replications per selected cell, seed %d+r. RNG: %s. Within each stage/replication all cells use the same 2000 uniforms for G, normal X values, uniforms for D, and standard-normal errors; n selects a prefix. G and D use inverse-CDF thresholds. This preserves the DGP marginals and provides common random numbers across both dimensions; it intentionally does not reuse the old pilot draw sequence.",
      config$design_value$grid_reps, config$seeds$design_grid,
      config$design_value$confirmation_reps, config$seeds$design_confirmation,
      paste(config$rng, collapse=" / ")), "",
    "## Mean losses and perfect-information gains", "",
    "Intervals are pointwise, asymptotic 95% Monte Carlo intervals, conservative for the minimum: construct Bonferroni simultaneous normal intervals for the two paired excess-loss means and take minima of endpoints. Proportional intervals apply the same construction to delta-method ratios. They remain valid at a tie between mean losses under the joint CLT; no differentiability of the minimum is assumed. They are not simultaneous across cells and do not account for selection of promising cells, hence independent confirmation. If every observed design favors one candidate, the observed gain and local MCSE are zero; the degenerate interval cannot rule out rare, unobserved reversals, so its upper endpoint is shown as NA in the combined table/report. CSV reports the exact one-sided 95% binomial upper bound on such a reversal probability (not its unbounded loss magnitude). The plug-in minimum can have downward Monte Carlo bias near equal mean losses; the simultaneous interval avoids assuming the selected branch is known. CSV also reports branch-local MCSEs and whether the fixed-choice gap is statistically separated. Normal intervals are also fragile when only one or a few reversals are observed. No resampling or additional seed is required.", "",
    design_value_table(main), "",
    "## Pairwise gaps and winner shares", "",
    "Delta=L2-L5; negative favors 2. Effective ties have |Delta|<=1e-10. Relative gap is |Delta|/min(L2,L5), not the proportional oracle gain. Detailed signed/absolute/relative quantiles and small-gap fractions are saved locally.", "",
    design_value_table(gaps), "",
    "## Bias and variance", "", design_value_table(components), "",
    "Population counterparts use the existing deterministic quadrature, depend on p_G but not n, and have no population Monte Carlo uncertainty. Analysis 2 remains conditionally unbiased. V5<=V2 still holds; squared conditional bias makes pooling imperfect.", "",
    "## Sample structure and failures", "", design_value_table(counts, 2L), "",
    sprintf("Across all runs: %d fit failures, %d rank failures, %d SE failures, %d violations of L2<=L3+1e-12. All five candidates are fitted in every design. Analysis 3 remains weakly dominated; it was not tuned to win. Invalid pairwise designs would stop gain reporting rather than be silently excluded.",
      sum(s$fit_failures), sum(s$rank_failures), sum(s$se_failures), sum(s$trimming_violations)), "",
    "## Full-menu context", "",
    "These percentages compare all five candidates, unlike the primary pairwise diagnostic. Analysis 4 can become attractive in sparse subgroups; no candidate is removed to improve the reported gain. Full-menu oracle shares condition on a complete menu, with the denominator shown.", "",
    design_value_table(overall, 2L), "",
    "## Diagnostic assessment", "",
    design_value_table(data.frame(stage=s$stage, n=s$n, p_G=s$p_G, gain_band=band), 2L), "",
    "These ranges describe effect size, not acceptance thresholds or an instruction to select the largest estimate. See `simulation_notes.md` for the scientific assessment and confirmation-cell rationale. No alternative is adopted automatically.", "",
    "All detailed outputs are regenerated locally: `output/tables/design_value_grid.csv`, `design_value_quantiles.csv`, `design_value_regions.csv`, `design_value_failures.csv`; per-cell paired draws in `output/results/design_value_*.csv`; `output/figures/design_value_grid.png`. Counts include full min/5th/median/95th/max distributions in the CSV. The Markdown report is tracked; CSVs and figures remain ignored.")
  writeLines(lines, "paper/design_value_diagnostic.md")
}

plot_design_value <- function(summary) {
  s <- summary[summary$stage == "grid", ]
  png("output/figures/design_value_grid.png", width=1200, height=620, res=130)
  on.exit(dev.off())
  par(mfrow=c(1,2), mar=c(4,4,3,1))
  colors <- c("#2166ac", "#b2182b", "#4d9221")
  for (metric in c("proportional", "absolute")) {
    if (metric == "proportional") {
      y <- 100*s$proportional_gain; low <- 100*s$proportional_lower; high <- 100*s$proportional_upper
    } else { y <- s$oracle_gain; low <- s$gain_lower; high <- s$gain_upper }
    plot(range(s$n)+c(-50,50), c(0,max(high, na.rm=TRUE)*1.05), type="n", xlab="Total sample size n",
      ylab=if (metric=="proportional") "Oracle gain (% of best fixed loss)" else "Oracle gain (squared-effect units)",
      main=if (metric=="proportional") "Attainable proportional gain" else "Attainable absolute gain")
    for (j in seq_along(config$design_value$p_G)) {
      i <- which(s$p_G==config$design_value$p_G[j]); i <- i[order(s$n[i])]
      x <- s$n[i]+(j-2)*12
      lines(x,y[i],type="b",pch=15+j,col=colors[j])
      arrows(x,low[i],x,high[i],angle=90,code=3,length=.04,col=colors[j])
    }
    legend("top",legend=paste("p_G =",config$design_value$p_G),col=colors,pch=16:18,lty=1,bty="n")
  }
}

# Prevent mixing cached stages from incompatible parameters/seeds/run sizes.
design_value_stage_config <- function(stage) {
  list(dgp=config$dgp, rng=config$rng, trim=config$trim,
       pairwise_tolerance=config$pairwise_tolerance, oracle_tolerance=config$oracle_tolerance,
       max_n=max(config$design_value$n),
       cells=if (stage=="grid") expand.grid(n=config$design_value$n, p_G=config$design_value$p_G) else config$design_value$confirmation_cells,
       B=if (stage=="grid") config$design_value$grid_reps else config$design_value$confirmation_reps,
       seed=if (stage=="grid") config$seeds$design_grid else config$seeds$design_confirmation)
}
