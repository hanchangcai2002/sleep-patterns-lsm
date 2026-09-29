# Simulation settings: true parameter values, scenario grid, MCMC settings.
#
# Goal (Section 1): check whether the location-scale model used in the
# paper -- location random effects (1 + dxgrp_factor | ID), scale random
# intercept (1 | ID) -- recovers group-specific between-subject variances
# when both groups share identical fixed effects, and compare it with the
# linear mixed-effects model.

# --- Design (Section 4) -----------------------------------------------------

n_days      <- 7          # T: days per subject, complete data
prop_sz     <- 25 / 63    # SZ share of subjects, as in the case study (38 HC : 25 SZ)
n_subjects  <- c(63, 150, 300)

# Ratios are SZ / HC on the VARIANCE scale; 1 = null (no group difference).
between_var_ratio <- c(1, 10)
within_var_ratio  <- c(1, 10)

scenarios <- expand.grid(
  N                 = n_subjects,
  between_var_ratio = between_var_ratio,
  within_var_ratio  = within_var_ratio
)
scenarios$scenario <- seq_len(nrow(scenarios))
scenarios <- scenarios[, c("scenario", "N", "between_var_ratio", "within_var_ratio")]

# --- True parameter values (Section 3) --------------------------------------
# Calibrated to the 7-day imputed case-study fit (brms), rounded.

truth <- list(
  beta0        = 370,               # HC mean TST on day 1 (min), ~ observed median
  beta_day     = rep(0, n_days),    # day effects, common to both groups (day 1 = 0)
  sigma_b_hc   = 60,                # HC between-subject SD (min)
  sigma_eps_hc = 60,                # HC within-subject SD (min), exp(tau0)
  sigma_omega  = 0.4                # SD of subject-level random scale effect omega_i
)
# Group fixed effects in the location model are zero: identical mean
# trajectories in both groups. Group differences enter only through
#   sigma_b_sz   = sigma_b_hc   * sqrt(between_var_ratio)   (Eq. 6: alpha1 = log ratio)
#   sigma_eps_sz = sigma_eps_hc * sqrt(within_var_ratio)    (tau2 = log(ratio) / 2)

# --- Simulation size and MCMC ------------------------------------------------

n_reps    <- 100          # replicates per scenario
seed_base <- 20260929

mcmc <- list(chains = 4, iter = 2000, warmup = 1000)

# Replicates run in parallel, one core per replicate (chains run sequentially
# within a replicate). Leave a few of the machine's cores free.
n_workers <- 10

# --- Paths --------------------------------------------------------------------

sim_dir     <- "simulation"
results_dir <- file.path(sim_dir, "results")
