# 7-day sensitivity model -- Section 4.3.3, Table 2

library(tidyverse)
library(readxl)
library(mice)
library(randomForest)
library(brms)

dat <- read_xlsx("data/SZ_obj_sleep_for_jinyuan.xlsx")

# --- Subject-level covariates: impute missing values ------------------------
# Random-forest imputation (mice, m = 5); the first completed dataset is used.

dat_covariates <- dat %>%
  select(subnum, psstot, phycomp, sapstot, sanstot, hba1c, agevisit, gender)

imp_covariates <- mice(dat_covariates, method = "rf", m = 5, maxit = 5, seed = 123)
completed_covariates <- complete(imp_covariates, 1) %>% rename(ID = subnum)

# Wide -> long, attach covariates, center step counts within/between subject.
prepare_long <- function(dat_tst, dat_steps, dxgrp, ID) {
  df_wide <- bind_cols(dat_tst, dat_steps)
  df_wide$dxgrp <- dxgrp
  df_wide$ID    <- ID

  df_wide %>%
    select(ID, dxgrp, starts_with("TST.")) %>%
    pivot_longer(starts_with("TST."), names_to = "time_k", values_to = "TST") %>%
    mutate(time_k = as.numeric(sub("TST\\.", "", time_k))) %>%
    left_join(
      df_wide %>%
        select(ID, starts_with("StepsCounts.")) %>%
        pivot_longer(-ID, names_to = "time_k", values_to = "StepsCounts") %>%
        mutate(time_k = as.numeric(sub("StepsCounts\\.", "", time_k))),
      by = c("ID", "time_k")
    ) %>%
    arrange(ID, time_k) %>%
    mutate(
      dxgrp_factor     = factor(dxgrp, levels = c(1, 2), labels = c("HC", "SZ")),
      Day              = time_k,
      StepsCounts_mean = ave(StepsCounts, ID, FUN = mean),
      StepsCounts_cw1  = StepsCounts - StepsCounts_mean,
      StepsCounts_cb1  = as.numeric(scale(StepsCounts, scale = FALSE)) - StepsCounts_cw1
    ) %>%
    left_join(completed_covariates, by = "ID")
}

# Location-scale model -- location Eq. 5, scale Eq. 8 -- with day modeled by
# a group-specific thin-plate spline; shared by the observed and imputed fits.
melsm_formula_7day <- bf(
  TST ~ 1 + s(Day, by = dxgrp_factor, k = 5) + dxgrp_factor +
    StepsCounts_cw1 + StepsCounts_cb1 +
    psstot + phycomp + sapstot + sanstot + hba1c + agevisit + gender +
    (1 + dxgrp_factor | ID),
  sigma ~ s(Day, by = dxgrp_factor, k = 5) + dxgrp_factor +
    StepsCounts_cw1 + StepsCounts_cb1 +
    sapstot + sanstot + (1 | ID)
)

# Comparison model: common (not group-specific) day spline.
melsm_formula_7day_noday <- bf(
  TST ~ 1 + s(Day, k = 5) + dxgrp_factor +
    StepsCounts_cw1 + StepsCounts_cb1 +
    psstot + phycomp + sapstot + sanstot + hba1c + agevisit + gender +
    (1 + dxgrp_factor | ID),
  sigma ~ s(Day, k = 5) + dxgrp_factor +
    StepsCounts_cw1 + StepsCounts_cb1 +
    sapstot + sanstot + (1 | ID)
)

fit_both <- function(df_long) {
  fit       <- brm(melsm_formula_7day,       data = df_long, init = 0, cores = 4,
                   sample_prior = TRUE, seed = 123)
  fit_noday <- brm(melsm_formula_7day_noday, data = df_long, init = 0, cores = 4,
                   sample_prior = TRUE, seed = 123)
  print(summary(fit))
  print(summary(fit_noday))
  print(loo_compare(loo(fit), loo(fit_noday)))
  print(bayes_R2(fit))
  print(bayes_R2(fit_noday))
  plot(conditional_effects(fit, effects = "Day:dxgrp_factor"))
  list(fit = fit, fit_noday = fit_noday)
}

dat_tst_7d   <- dat %>% select(starts_with("TST."))         %>% select(1:7)
dat_steps_7d <- dat %>% select(starts_with("StepsCounts.")) %>% select(1:7)

# --- Complete-case data -- Table 2, "Observed" columns ----------------------

complete_idx <- complete.cases(bind_cols(dat_tst_7d, dat_steps_7d))

df_long_complete <- prepare_long(
  dat_tst_7d[complete_idx, ], dat_steps_7d[complete_idx, ],
  dat$dxgrp[complete_idx], dat$subnum[complete_idx]
)

fits_complete <- fit_both(df_long_complete)

# --- Multiply-imputed data -- Table 2, "Imputed" columns --------------------
# TST and step counts imputed separately with mice (random forest, m = 5),
# with diagnosis as a predictor and using only earlier days to impute later
# days; the first completed dataset is used below.

impute_sequential <- function(dat_days, dxgrp, timepoints, seed = 123) {
  dat_days$dxgrp <- dxgrp
  init <- mice(dat_days, maxit = 0)
  meth <- init$method
  meth["dxgrp"] <- ""
  meth[names(meth)[meth != ""]] <- "rf"
  predM <- init$predictorMatrix
  for (i in seq_along(timepoints)[-length(timepoints)]) {
    predM[timepoints[i], timepoints[(i + 1):length(timepoints)]] <- 0
  }
  mice(dat_days, method = meth, predictorMatrix = predM, m = 5, maxit = 5, seed = seed)
}

imputed_tst   <- impute_sequential(dat_tst_7d,   dat$dxgrp, paste0("TST.", 1:7))
imputed_steps <- impute_sequential(dat_steps_7d, dat$dxgrp, paste0("StepsCounts.", 1:7))

df_long_imputed <- prepare_long(
  complete(imputed_tst, 1)   %>% select(starts_with("TST.")),
  complete(imputed_steps, 1) %>% select(starts_with("StepsCounts.")),
  dat$dxgrp, dat$subnum
)

fits_imputed <- fit_both(df_long_imputed)
