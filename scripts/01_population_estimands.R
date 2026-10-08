source("R/simulation_helpers.R")
load_project()
ensure_output()
population <- population_estimands()
refined <- population_estimands(bound = 12, tolerance = 1e-12)
population$refinement_difference <- refined$theta_a - population$theta_a
stopifnot(max(abs(population$refinement_difference)) < 1e-8)
write_csv(population, "output/tables/population_estimands.csv")
print(population, digits = 7, row.names = FALSE)
