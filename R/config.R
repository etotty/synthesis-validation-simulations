# Final DGP parameters, run sizes, and all random seeds live here.
# Calibration cases and numerical reporting tolerances are explicit below.
config <- list(
  dgp = list(n = 1500L, p_G = 0.12,
             alpha = c(-0.3, 0.4, 0.9, 1.0),
             beta = c(0, 0.7, 0.5, 0.3), # intercept, X, X^2, G
             gamma = 0.10, tau0 = 0.25, delta = 0.25, sigma = 1),
  trim = c(0.10, 0.90),
  gamma_grid = c(0, 0.05, 0.10, 0.15, 0.20),
  calibration_effects = list(v01 = c(tau0 = 0.33, delta = 0.17),
                             restored = c(tau0 = 0.25, delta = 0.25)),
  pairwise_tolerance = 1e-10,
  oracle_tolerance = 1e-12,
  pilot_reps = 2000L,
  calibration_reps = 300L,
  # v0.1.2 diagnostic only: never replaces dgp above.
  design_value = list(n = c(1000L, 1500L, 2000L), p_G = c(.08, .12, .16),
                      grid_reps = 3000L, confirmation_reps = 10000L,
                      confirmation_cells = data.frame(n = c(1500L, 2000L, 1000L),
                                                      p_G = c(.12, .08, .16))),
  # Replication r uses base + r. Stages have disjoint seed ranges.
  seeds = list(calibration = 100000L, pilot = 200000L,
               test = 300000L, oracle_check = 400000L,
               design_grid = 500000L, design_confirmation = 600000L),
  oracle_check_reps = 20000L,
  rng = c("Mersenne-Twister", "Inversion", "Rejection")
)
