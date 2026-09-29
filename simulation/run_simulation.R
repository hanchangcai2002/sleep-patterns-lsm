# Run the simulation. From the repository root:
#
#   Rscript simulation/run_simulation.R                         # all scenarios, all replicates
#   Rscript simulation/run_simulation.R --scenarios 1,4 --reps 1-2   # pilot / subset
#
# Each replicate is saved to results/scenario_XX/rep_XXX.rds as soon as it
# finishes; replicates already on disk are skipped, so an interrupted run can
# simply be restarted.

source("simulation/00_settings.R")
source("simulation/01_generate_data.R")
source("simulation/02_fit_models.R")
library(parallel)

# --- Command-line arguments ---------------------------------------------------

parse_range <- function(x) {
  unlist(lapply(strsplit(x, ",")[[1]], function(p) {
    r <- as.integer(strsplit(p, "-")[[1]])
    if (length(r) == 2) seq(r[1], r[2]) else r
  }))
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) default else parse_range(args[i + 1])
}
run_scenarios <- get_arg("--scenarios", scenarios$scenario)
run_reps      <- get_arg("--reps",      seq_len(n_reps))

# --- Setup ----------------------------------------------------------------------

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(scenarios, file.path(results_dir, "scenarios.csv"), row.names = FALSE)

rep_seed <- function(scenario, rep) seed_base + 1000 * scenario + rep
rep_file <- function(scenario, rep) {
  file.path(results_dir, sprintf("scenario_%02d", scenario), sprintf("rep_%03d.rds", rep))
}

set.seed(seed_base)
message("Compiling the LSM ...")
lsm_template <- compile_lsm(
  simulate_data(63, 1, 1, truth, n_days, prop_sz)
)

# --- One replicate --------------------------------------------------------------

run_one <- function(sc, rep) {
  seed <- rep_seed(sc$scenario, rep)
  set.seed(seed)
  d <- simulate_data(sc$N, sc$between_var_ratio, sc$within_var_ratio, truth, n_days, prop_sz)

  res <- list(
    scenario  = sc,
    rep       = rep,
    seed      = seed,
    truth     = true_values(sc$between_var_ratio, sc$within_var_ratio, truth),
    data      = d,
    lme_ri    = fit_lme(d, random = ~ 1 | ID),
    lme_paper = fit_lme(d, random = ~ dxgrp_factor | ID),
    lsm       = fit_lsm(lsm_template, d, seed, mcmc)
  )
  saveRDS(res, rep_file(sc$scenario, rep))
  res$lsm$time_sec
}

# --- Main loop -------------------------------------------------------------------

for (s in run_scenarios) {
  sc <- scenarios[scenarios$scenario == s, ]
  dir.create(dirname(rep_file(s, 1)), showWarnings = FALSE)
  todo <- run_reps[!file.exists(rep_file(s, run_reps))]
  message(sprintf("Scenario %d (N = %d, between ratio = %g, within ratio = %g): %d replicates to run",
                  s, sc$N, sc$between_var_ratio, sc$within_var_ratio, length(todo)))
  if (length(todo) == 0) next

  t0 <- proc.time()[["elapsed"]]
  times <- mclapply(todo, function(r) run_one(sc, r), mc.cores = n_workers, mc.preschedule = FALSE)
  failed <- vapply(times, inherits, logical(1), "try-error")
  message(sprintf("  done in %.1f min; mean LSM fit %.1f s; %d worker errors",
                  (proc.time()[["elapsed"]] - t0) / 60,
                  mean(unlist(times[!failed])), sum(failed)))
}
