# Standalone v0.2 benchmark, from the repository root. No calibration changes.
source("R/simulation_helpers.R")
load_project()
ensure_output()
args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args)) args[1] else "all"
stopifnot(mode %in% c("all",config$researcher$cases$case,"report"))
original_dgp <- config$dgp
if (mode!="report") for (i in seq_len(nrow(config$researcher$cases))) {
  case <- config$researcher$cases$case[i]
  if (mode!="all" && mode!=case) next
  p <- config$dgp
  p$n <- config$researcher$cases$n[i]; p$p_G <- config$researcher$cases$p_G[i]
  out <- run_researcher_case(case,p)
  for (name in names(out)) write_csv(out[[name]],
    paste0("output/tables/researcher_",case,"_",name,".csv"))
  dput(list(dgp=p, rule=config$researcher, seeds=config$seeds[c("researcher","researcher_fresh_noise")],
            rng=config$rng, trim=config$trim),
       file=paste0("output/results/researcher_",case,"_config.R"))
}
# Report-only mode does not rerun simulations; incompatible caches fail loudly.
for (case in config$researcher$cases$case) {
  path <- paste0("output/results/researcher_",case,"_config.R")
  if (!file.exists(path)) {
    message("Case ",case," is not yet complete; run the remaining case before report.")
    quit(save="no",status=0)
  }
  p <- config$dgp; i <- match(case,config$researcher$cases$case)
  p$n <- config$researcher$cases$n[i]; p$p_G <- config$researcher$cases$p_G[i]
  expected <- list(dgp=p,rule=config$researcher,
    seeds=config$seeds[c("researcher","researcher_fresh_noise")],rng=config$rng,trim=config$trim)
  stopifnot(identical(dget(path),expected))
}
for (name in c("policy","diagnostic","confusion")) {
  combined <- do.call(rbind,lapply(config$researcher$cases$case,function(case)
    read.csv(paste0("output/tables/researcher_",case,"_",name,".csv"))))
  write_csv(combined,paste0("output/tables/researcher_",name,".csv"))
  assign(name,combined)
}
stopifnot(identical(config$dgp,original_dgp))
write_researcher_report(policy,diagnostic,confusion)
plot_researcher_results(policy)
capture.output(sessionInfo(),file="output/results/researcher_session_info.txt")
message("Researcher-information benchmark complete; canonical DGP unchanged.")
