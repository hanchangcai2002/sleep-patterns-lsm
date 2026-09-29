# Models fitted to each simulated dataset (Section 5):
#   1. LME, random intercept only (conventional homoscedastic model)
#   2. LME as in the paper, random = ~ dxgrp_factor | ID (Supplemental Tab. 2)
#   3. LSM as in the paper: location (1 + dxgrp_factor | ID), scale (1 | ID),
#      Day x diagnosis in both submodels (4-day model structure, 7 days),
#      without subject-level covariates.
#
# Each fit is reduced to a compact summary; the full brmsfit is not kept.

library(nlme)
library(brms)
library(posterior)

# --- LME -----------------------------------------------------------------------

summarise_lme <- function(fit) {
  tt <- summary(fit)$tTable
  vc <- getVarCov(fit)
  sd_int <- sqrt(vc[1, 1])
  out <- list(
    fixed     = data.frame(term = rownames(tt), estimate = tt[, "Value"], se = tt[, "Std.Error"],
                           df = tt[, "DF"], p = tt[, "p-value"], row.names = NULL),
    sigma_b_hc = sd_int,
    sigma_eps  = fit$sigma
  )
  if (nrow(vc) == 2) {
    out$sigma_b_sz <- sqrt(vc[1, 1] + vc[2, 2] + 2 * vc[1, 2])
    out$sd_slope   <- sqrt(vc[2, 2])
    out$cor        <- vc[1, 2] / (sd_int * sqrt(vc[2, 2]))
  }
  out
}

fit_lme <- function(d, random) {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch({
    fit <- lme(TST ~ Day * dxgrp_factor, random = random, method = "REML", data = d)
    c(list(ok = TRUE, error = NA_character_), summarise_lme(fit))
  }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  res$time_sec <- proc.time()[["elapsed"]] - t0
  res
}

# --- LSM (brms) ---------------------------------------------------------------

lsm_formula <- bf(
  TST   ~ 1 + Day * dxgrp_factor + (1 + dxgrp_factor | ID),
  sigma ~ Day * dxgrp_factor + (1 | ID)
)

# Fixed priors, set to what brms chooses by default on the case-study TST
# (median 368, MAD 112.7 min). brms would otherwise recompute them from each
# simulated dataset, which changes the Stan code and forces a recompile.
lsm_priors <- c(
  prior(student_t(3, 368, 112.7), class = Intercept),
  prior(student_t(3, 0, 112.7),   class = sd),
  prior(lkj(1),                   class = cor),
  prior(student_t(3, 0, 2.5),     class = Intercept, dpar = sigma),
  prior(student_t(3, 0, 2.5),     class = sd,        dpar = sigma)
)

# Compile once; every replicate reuses this model through update().
compile_lsm <- function(d) {
  brm(lsm_formula, data = d, prior = lsm_priors, sample_prior = TRUE,
      chains = 0, silent = 2, refresh = 0)
}

lsm_summary_vars <- function(fit) {
  v <- variables(fit)
  v[grepl("^(b_|sd_|cor_)", v)]
}

summarise_lsm <- function(fit) {
  draws <- as_draws_df(fit)

  sd0 <- draws$sd_ID__Intercept
  sd1 <- draws$sd_ID__dxgrp_factorSZ
  rho <- draws$cor_ID__Intercept__dxgrp_factorSZ

  # Group-specific between-subject SDs implied by (1 + dxgrp_factor | ID).
  derived <- draws_df(
    sigma_b_hc  = sd0,
    sigma_b_sz  = sqrt(sd0^2 + sd1^2 + 2 * rho * sd0 * sd1),
    alpha1      = log((sd0^2 + sd1^2 + 2 * rho * sd0 * sd1) / sd0^2),
    tau0        = draws$b_sigma_Intercept,
    tau2        = draws$b_sigma_dxgrp_factorSZ,
    sigma_omega = draws$sd_ID__sigma_Intercept
  )

  summ_fun <- list("mean", "median", "sd",
                   q2.5 = ~ quantile(.x, 0.025), q97.5 = ~ quantile(.x, 0.975),
                   "rhat", "ess_bulk", "ess_tail")
  params  <- summarise_draws(subset_draws(draws, variable = lsm_summary_vars(fit)), summ_fun)
  derived <- summarise_draws(derived, summ_fun)

  # Prior draws (sample_prior = TRUE), for prior-vs-posterior comparison of
  # the random-effect hyperparameters.
  pd <- prior_draws(fit)
  prior_q <- t(sapply(pd, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))

  np   <- nuts_params(fit)
  rh   <- rhat(fit)
  ess  <- neff_ratio(fit)
  list(
    params  = as.data.frame(params),
    derived = as.data.frame(derived),
    prior   = data.frame(variable = rownames(prior_q), prior_q, row.names = NULL, check.names = FALSE),
    diagnostics = list(
      n_divergent   = sum(np$Value[np$Parameter == "divergent__"]),
      n_max_treedepth = sum(np$Value[np$Parameter == "treedepth__"] >= 10),
      max_rhat      = max(rh, na.rm = TRUE),
      n_rhat_gt_101 = sum(rh > 1.01, na.rm = TRUE),
      min_neff_ratio = min(ess, na.rm = TRUE)
    )
  )
}

fit_lsm <- function(template, d, seed, mcmc) {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch({
    fit <- update(template, newdata = d, recompile = FALSE,
                  chains = mcmc$chains, iter = mcmc$iter, warmup = mcmc$warmup,
                  cores = 1, init = 0, seed = seed, refresh = 0, silent = 2)
    c(list(ok = TRUE, error = NA_character_), summarise_lsm(fit))
  }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  res$time_sec <- proc.time()[["elapsed"]] - t0
  res
}
