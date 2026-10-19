# Compare two independently generated output directories (CSV and PNG), plus
# the generated paper table. Ignore session_info: it describes each environment.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 2:3)
mode <- if (length(args)==3L) args[3] else "all"
stopifnot(mode %in% c("all", "design-value", "researcher"))
generated_paths <- function(root) {
  reports <- switch(mode, all=c("pilot_results.md", "design_value_diagnostic.md", "researcher_information_results.md"),
                    `design-value`="design_value_diagnostic.md", researcher="researcher_information_results.md")
  reports <- reports[file.exists(file.path(root, "paper", reports))]
  outputs <- list.files(file.path(root, "output"), pattern = "\\.(csv|png|R)$",
                        recursive = TRUE, full.names = FALSE)
  if (mode=="design-value") {
    outputs <- outputs[grepl("^design_value_", basename(outputs))]
    stopifnot(length(outputs)>0L, "design_value_diagnostic.md" %in% reports)
  }
  if (mode=="researcher") {
    outputs <- outputs[grepl("^researcher_", basename(outputs))]
    stopifnot(length(outputs)>0L, "researcher_information_results.md" %in% reports)
  }
  sort(c(outputs, paste0("../paper/", reports)))
}
paths <- generated_paths(args[1])
other <- generated_paths(args[2])
stopifnot(identical(paths, other))
a <- unname(tools::md5sum(file.path(args[1], "output", paths)))
b <- unname(tools::md5sum(file.path(args[2], "output", paths)))
stopifnot(!anyNA(a), !anyNA(b))
if (!identical(a, b)) stop("Different files: ", paste(paths[a != b], collapse = ", "))
cat("PASS: identical hashes for", length(paths), "generated tables, draws, figures, config and report.\n")
