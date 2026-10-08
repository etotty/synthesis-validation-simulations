source("R/simulation_helpers.R")
load_project()
expect_close <- function(actual, expected, tolerance = 1e-10) {
  stopifnot(length(actual) == length(expected), all(is.finite(actual)),
            all(is.finite(expected)), max(abs(actual - expected)) < tolerance)
}
passed <- character()
for (file in c("test_dgp.R", "test_analysis_menu.R", "test_oracle_loss.R", "test_reproducibility.R")) {
  source(file.path("tests", file), local = new.env(parent = .GlobalEnv))
  passed <- c(passed, file)
  message("PASS: ", file)
}
message("All ", length(passed), " test suites passed.")
