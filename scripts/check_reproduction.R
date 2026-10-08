# Compare two independently generated output directories (CSV and PNG), plus
# the generated paper table. Ignore session_info: it describes each environment.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
paths <- sort(c(list.files(file.path(args[1], "output"), pattern = "\\.(csv|png|R)$",
                          recursive = TRUE, full.names = FALSE), "../paper/pilot_results.md"))
other <- sort(c(list.files(file.path(args[2], "output"), pattern = "\\.(csv|png|R)$",
                          recursive = TRUE, full.names = FALSE), "../paper/pilot_results.md"))
stopifnot(identical(paths, other))
a <- unname(tools::md5sum(file.path(args[1], "output", paths)))
b <- unname(tools::md5sum(file.path(args[2], "output", paths)))
if (!identical(a, b)) stop("Different files: ", paste(paths[a != b], collapse = ", "))
cat("PASS: identical hashes for", length(paths), "generated tables, draws, figures, config and report.\n")
