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
