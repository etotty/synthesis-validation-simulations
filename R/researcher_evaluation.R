# Simulation evaluation ONLY. Truth must not be passed to the observed-data API.
attainable_fraction <- function(D_fixed, D_feasible, D_oracle) {
  denominator <- D_fixed-D_oracle
  if (denominator==0) return(NA_real_)
  (D_fixed-D_feasible)/denominator  # deliberately no clipping
}

researcher_replication <- function(r, p, theta_a) {
  max_n <- max(config$researcher$cases$n)
  seed <- config$seeds$researcher+r
  dat <- design_value_sample(p, design_value_draws(seed, max_n))
  # This explicit projection is the only input crossing the information boundary.
  observed <- dat[,c("Y","X","D","G")]
  feasible <- researcher_choice(observed, config$researcher$precision,
                                config$researcher$tie_tolerance)
  info <- feasible$information; decision <- feasible$decision
  all <- run_replication(r, p, theta_a, config$seeds$researcher, dat=dat)
  if (!decision$ok) stop("Feasible fit failed in replication ", r, ": ", decision$status)
  if (any(all$status!="ok")) stop("Menu fit failure in replication ",r,"; no silent exclusions")
  hc3 <- perceived_loss(info, "HC3", config$researcher$tie_tolerance)
  if (!hc3$ok) stop("HC3 diagnostic failed in replication ",r)
  pair <- pairwise_designs(all)
  stopifnot(pair$valid)
  # Prefix every truth-based diagnostic inherited from the old pairwise output.
  names(pair) <- paste0("evaluation_", names(pair))
  out <- cbind(pair, data.frame(
    gamma_hat=info$gamma_hat, gamma_se=info$gamma_se, sigma2_hat=info$sigma2_hat,
    c5_GX=info$c5_GX, estimated_bias_5=info$estimated_bias_5,
    perceived_loss_2=decision$perceived_loss_2, perceived_loss_5=decision$perceived_loss_5,
    perceived_bias_penalty_5=decision$perceived_bias_penalty_5,
    variance_proxy_2=decision$variance_proxy_2, variance_proxy_5=decision$variance_proxy_5,
    feasible_choice=decision$choice, feasible_tied=decision$tied,
    hc3_variance_proxy_2=hc3$variance_proxy_2, hc3_variance_proxy_5=hc3$variance_proxy_5,
    hc3_perceived_loss_2=hc3$perceived_loss_2, hc3_perceived_loss_5=hc3$perceived_loss_5,
    hc3_choice=hc3$choice, evaluation_true_gamma=p$gamma,
    evaluation_theta=target_effect(p),
    evaluation_oracle_gamma_variance=p$sigma^2*sum(info$gamma_h^2),
    evaluation_ovb_error=(info$T_5-info$T_augmented)-info$estimated_bias_5))
  # Independent outcome-noise evaluation, never an input to the choice.
  # Same Z and already-made choice; this is not sample splitting or a new regime.
  set.seed(config$seeds$researcher_fresh_noise+r)
  fresh_noise <- rnorm(max_n, sd=p$sigma)[seq_len(p$n)]
  for (a in c(2L,5L)) {
    row <- all[all$analysis==a, ]
    out[[paste0("oracle_loss_",a)]] <- row$L
    out[[paste0("oracle_variance_",a)]] <- row$V
    out[[paste0("oracle_bias_",a)]] <- row$conditional_bias
    out[[paste0("evaluation_theta_",a)]] <- row$theta_a
    out[[paste0("T_",a)]] <- row$T
    out[[paste0("evaluation_epsilon_",a)]] <- row$epsilon
    out[[paste0("evaluation_design_offset_",a)]] <- row$conditional_mean_T-row$theta_a
    out[[paste0("evaluation_noise_",a)]] <- row$T-row$conditional_mean_T
    h <- if (a==2L) info$h2 else info$h5
    fresh <- sum(h*fresh_noise)
    out[[paste0("evaluation_fresh_noise_",a)]] <- fresh
    out[[paste0("evaluation_fresh_T_",a)]] <- row$conditional_mean_T+fresh
    out[[paste0("evaluation_gamma_noise_cov_",a)]] <- p$sigma^2*sum(info$gamma_h*h)
  }
  out$oracle_choice <- if (out$oracle_loss_2<=out$oracle_loss_5) 2L else 5L
  out$oracle_bias_penalty_5 <- out$oracle_bias_5^2
  out$perceived_gap <- out$perceived_loss_2-out$perceived_loss_5
  out$oracle_loss_gap <- out$oracle_loss_2-out$oracle_loss_5
  out$evaluation_expected_penalty_inflation <- info$c5_GX^2*out$evaluation_oracle_gamma_variance
  list(pair=out, menu=all)
}

researcher_performance <- function(d, choice) {
  stopifnot(length(choice)==nrow(d), all(choice %in% c(2L,5L)))
  L2 <- d$oracle_loss_2; L5 <- d$oracle_loss_5
  fixed <- if (mean(L2)<=mean(L5)) L2 else L5
  selected <- ifelse(choice==2L,L2,L5)
  oracle <- pmin(L2,L5)
  numerator <- fixed-selected; denominator <- fixed-oracle
  fraction <- attainable_fraction(mean(fixed), mean(selected), mean(oracle))
  se <- function(x) stats::sd(x)/sqrt(length(x))
  fraction_se <- if (mean(denominator)>0)
    se((numerator-fraction*denominator)/mean(denominator)) else NA_real_
  data.frame(D_fixed=mean(fixed), D_feasible=mean(selected), D_oracle=mean(oracle),
    D_feasible_mcse=se(selected), attainable_gain=mean(denominator),
    attainable_gain_mcse=se(denominator), excess_over_oracle=mean(selected-oracle),
    excess_mcse=se(selected-oracle), improvement_over_fixed=mean(numerator),
    improvement_mcse=se(numerator), fraction_captured=fraction,
    fraction_mcse=fraction_se,
    denominator_fragile=mean(denominator)<=3*se(denominator),
    agreement=mean(choice==d$oracle_choice),
    agreement_mcse=se(as.numeric(choice==d$oracle_choice)), choice_5=mean(choice==5L))
}

researcher_policy_summary <- function(d) {
  policies <- list(always_2=rep(2L,nrow(d)), always_5=rep(5L,nrow(d)),
    oracle_design=d$oracle_choice, feasible=d$feasible_choice, HC3_diagnostic=d$hc3_choice)
  se <- function(x) stats::sd(x)/sqrt(length(x))
  do.call(rbind,lapply(names(policies),function(policy) {
    choice <- policies[[policy]]
    pick <- function(prefix) ifelse(choice==2L,d[[paste0(prefix,2)]],d[[paste0(prefix,5)]])
    systematic <- pick("evaluation_theta_")-d$evaluation_theta
    epsilon <- pick("evaluation_epsilon_")
    eta <- pick("evaluation_noise_")
    design <- pick("evaluation_design_offset_")
    T <- pick("T_")
    mse <- (T-d$evaluation_theta)^2
    fresh_mse <- (pick("evaluation_fresh_T_")-d$evaluation_theta)^2
    L <- pick("oracle_loss_")
    # Paired comparison with always-2 removes its mean-noise Monte Carlo error.
    noise_selection <- eta-d$evaluation_noise_2
    data.frame(policy=policy, researcher_performance(d,choice),
      systematic_distortion=mean(systematic), systematic_mcse=se(systematic),
      epsilon_distortion=mean(epsilon), epsilon_mcse=se(epsilon),
      design_offset=mean(design), design_offset_mcse=se(design),
      outcome_noise_distortion=mean(eta), outcome_noise_mcse=se(eta),
      noise_selection_vs_2=mean(noise_selection), noise_selection_mcse=se(noise_selection),
      total_bias=mean(T-d$evaluation_theta), total_bias_mcse=se(T-d$evaluation_theta),
      actual_mse=mean(mse), actual_mse_mcse=se(mse),
      independent_noise_mse=mean(fresh_mse), independent_noise_mse_mcse=se(fresh_mse),
      actual_mse_minus_design_loss=mean(mse-L), mse_minus_loss_mcse=se(mse-L),
      fresh_noise_distortion=mean(pick("evaluation_fresh_noise_")),
      fresh_noise_mcse=se(pick("evaluation_fresh_noise_")))
  }))
}

researcher_diagnostics <- function(d) {
  se <- function(x) stats::sd(x)/sqrt(length(x))
  cor <- function(x,y) if (stats::sd(x)>0 && stats::sd(y)>0) stats::cor(x,y) else NA_real_
  gamma_error <- d$gamma_hat-d$evaluation_true_gamma
  I5 <- as.numeric(d$feasible_choice==5L)
  data.frame(B=nrow(d), gamma_mean=mean(d$gamma_hat), gamma_mean_mcse=se(d$gamma_hat),
    gamma_sd=stats::sd(d$gamma_hat), gamma_rmse=sqrt(mean(gamma_error^2)),
    gamma_q025=unname(stats::quantile(d$gamma_hat,.025)),
    gamma_median=stats::median(d$gamma_hat), gamma_q975=unname(stats::quantile(d$gamma_hat,.975)),
    gamma_mean_se=mean(d$gamma_se), gamma_negative=mean(d$gamma_hat<0),
    gamma_coverage=mean(abs(gamma_error)<=1.96*d$gamma_se),
    sigma2_mean=mean(d$sigma2_hat),
    bias_hat_mean=mean(d$estimated_bias_5), oracle_bias_mean=mean(d$oracle_bias_5),
    bias_tracking_correlation=cor(d$estimated_bias_5,d$oracle_bias_5),
    bias_rmse=sqrt(mean((d$estimated_bias_5-d$oracle_bias_5)^2)),
    perceived_gap_mean=mean(d$perceived_gap), oracle_gap_mean=mean(d$oracle_loss_gap),
    loss_gap_correlation=cor(d$perceived_gap,d$oracle_loss_gap),
    mean_estimated_penalty=mean(d$perceived_bias_penalty_5),
    mean_oracle_penalty=mean(d$oracle_bias_penalty_5),
    mean_expected_penalty_inflation=mean(d$evaluation_expected_penalty_inflation),
    mean_variance_proxy_2=mean(d$variance_proxy_2), mean_variance_proxy_5=mean(d$variance_proxy_5),
    mean_oracle_variance_2=mean(d$oracle_variance_2), mean_oracle_variance_5=mean(d$oracle_variance_5),
    choice_cor_epsilon_2=cor(I5,d$evaluation_epsilon_2),
    choice_cor_epsilon_5=cor(I5,d$evaluation_epsilon_5),
    choice_cor_noise_2=cor(I5,d$evaluation_noise_2),
    choice_cor_noise_5=cor(I5,d$evaluation_noise_5),
    choice_cor_noise_difference=cor(I5,d$evaluation_noise_5-d$evaluation_noise_2),
    gamma_error_cor_noise_2=cor(gamma_error,d$evaluation_noise_2),
    gamma_error_cor_noise_5=cor(gamma_error,d$evaluation_noise_5),
    max_ovb_error=max(abs(d$evaluation_ovb_error)),
    max_bias_identity_error=max(abs(d$oracle_bias_5-d$evaluation_true_gamma*d$c5_GX)),
    feasible_ties=sum(d$feasible_tied),
    subgroup_mean=mean(d$evaluation_n_subgroup),
    subgroup_q05=unname(stats::quantile(d$evaluation_n_subgroup,.05)),
    subgroup_q95=unname(stats::quantile(d$evaluation_n_subgroup,.95)),
    treated_mean=mean(d$evaluation_n_subgroup_treated), control_mean=mean(d$evaluation_n_subgroup_control))
}

researcher_confusion <- function(d) {
  table <- as.data.frame(table(oracle=factor(d$oracle_choice,levels=c(2,5)),
                              feasible=factor(d$feasible_choice,levels=c(2,5))))
  names(table)[3] <- "count"
  table$fraction <- table$count/nrow(d)
  table
}

run_researcher_case <- function(case, p, B=config$researcher$reps) {
  population <- population_estimands(p)
  pair <- vector("list",B); menu <- vector("list",B)
  for (r in seq_len(B)) {
    result <- researcher_replication(r,p,population$theta_a)
    pair[[r]] <- result$pair; menu[[r]] <- result$menu
    if (r %% 1000L==0L) message(case, ": ",r,"/",B)
  }
  pair <- do.call(rbind,pair); menu <- do.call(rbind,menu)
  stopifnot(all(menu$epsilon==menu$T-menu$theta_a), all(menu$L>=0),
            max(abs(pair$evaluation_ovb_error))<1e-8,
            max(abs(pair$oracle_bias_5-pair$evaluation_true_gamma*pair$c5_GX))<1e-8)
  write_csv(pair,paste0("output/results/researcher_",case,"_draws.csv"))
  write_csv(menu,paste0("output/results/researcher_",case,"_menu.csv"))
  write_csv(population,paste0("output/tables/researcher_",case,"_population.csv"))
  write_csv(summarize_pilot(menu),paste0("output/tables/researcher_",case,"_menu_summary.csv"))
  list(policy=cbind(case=case,researcher_policy_summary(pair)),
       diagnostic=cbind(case=case,researcher_diagnostics(pair)),
       confusion=cbind(case=case,researcher_confusion(pair)))
}
