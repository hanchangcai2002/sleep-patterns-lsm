# Data-generating model (Section 2): a covariate-free version of the paper's
# location-scale model.
#
#   y_ik          = beta0 + beta_day[k] + b_i + eps_ik
#   b_i           ~ N(0, sigma_b[g(i)]^2)                  between-subject, by group
#   eps_ik        ~ N(0, sigma_ik^2)
#   log sigma_ik  = tau0 + tau2 * x_i + omega_i            within-subject
#   omega_i       ~ N(0, sigma_omega^2)                    shared by both groups
#
# with x_i = 1 for SZ, and identical fixed effects in both groups.

true_values <- function(between_var_ratio, within_var_ratio, truth) {
  list(
    beta0        = truth$beta0,
    beta_day     = truth$beta_day,
    sigma_b_hc   = truth$sigma_b_hc,
    sigma_b_sz   = truth$sigma_b_hc * sqrt(between_var_ratio),
    alpha0       = log(truth$sigma_b_hc^2),
    alpha1       = log(between_var_ratio),
    tau0         = log(truth$sigma_eps_hc),
    tau2         = log(within_var_ratio) / 2,
    sigma_omega  = truth$sigma_omega
  )
}

simulate_data <- function(N, between_var_ratio, within_var_ratio, truth, n_days, prop_sz) {
  tv <- true_values(between_var_ratio, within_var_ratio, truth)

  n_sz <- round(N * prop_sz)
  x    <- c(rep(0, N - n_sz), rep(1, n_sz))

  b     <- rnorm(N, 0, ifelse(x == 1, tv$sigma_b_sz, tv$sigma_b_hc))
  omega <- rnorm(N, 0, tv$sigma_omega)

  d <- expand.grid(time_k = seq_len(n_days), ID = seq_len(N))[, c("ID", "time_k")]
  xi        <- x[d$ID]
  sigma_eps <- exp(tv$tau0 + tv$tau2 * xi + omega[d$ID])

  d$TST          <- tv$beta0 + tv$beta_day[d$time_k] + b[d$ID] + rnorm(nrow(d), 0, sigma_eps)
  d$Day          <- factor(d$time_k)
  d$dxgrp_factor <- factor(xi, levels = c(0, 1), labels = c("HC", "SZ"))
  d
}
