# 4-day model -- Section 4.2, Table 1

library(tidyverse)
library(readxl)
library(mice)
library(randomForest)
library(nlme)
library(brms)

dat <- read_xlsx("data/SZ_obj_sleep_for_jinyuan.xlsx")

# --- Subject-level covariates: impute missing values ------------------------
# Random-forest imputation (mice, m = 5); the first completed dataset is used.

dat_covariates <- dat %>%
  select(subnum, psstot, phycomp, sapstot, sanstot, hba1c, agevisit, gender)

imp_covariates <- mice(dat_covariates, method = "rf", m = 5, maxit = 5, seed = 123)
completed_covariates <- complete(imp_covariates, 1) %>% rename(ID = subnum)

# --- Prepare data: days 1-4 (all complete), wide -> long, centering ---------

dat_tst   <- dat %>% select(starts_with("TST."))         %>% select(1:4)
dat_steps <- dat %>% select(starts_with("StepsCounts.")) %>% select(1:4)

df_wide <- bind_cols(dat_tst, dat_steps)
df_wide$dxgrp <- dat$dxgrp
df_wide$ID    <- dat$subnum

df_long <- df_wide %>%
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
    Day              = factor(time_k),
    StepsCounts_mean = ave(StepsCounts, ID, FUN = mean),
    StepsCounts_cw1  = StepsCounts - StepsCounts_mean,
    StepsCounts_cb1  = as.numeric(scale(StepsCounts, scale = FALSE)) - StepsCounts_cw1
  ) %>%
  left_join(completed_covariates, by = "ID")

# --- Linear mixed-effects model (Eq. 2, Section 4.2.1) ---------------------

fit_lme <- lme(
  fixed  = TST ~ Day * dxgrp_factor + StepsCounts_cw1 + StepsCounts_cb1,
  random = ~ dxgrp_factor | ID,
  method = "REML",
  data   = df_long
)
summary(fit_lme)

# --- Location-scale model -- location Eq. 4, between-subject variance ------
# --- Eq. 6, scale Eq. 7 -- Table 1 ------------------------------------------

fit_lsm <- brm(
  bf(
    TST ~ 1 + Day * dxgrp_factor + StepsCounts_cw1 + StepsCounts_cb1 +
      psstot + phycomp + sapstot + sanstot + hba1c + agevisit + gender +
      (1 + dxgrp_factor | ID),
    sigma ~ Day * dxgrp_factor + StepsCounts_cw1 + StepsCounts_cb1 +
      sapstot + sanstot + (1 | ID)
  ),
  data         = df_long,
  init         = 0,
  cores        = 4,
  sample_prior = TRUE,
  seed         = 123
)
summary(fit_lsm)

# --- Comparison model without the day-by-diagnosis interaction -------------

fit_lsm_noday <- brm(
  bf(
    TST ~ 1 + Day + dxgrp_factor + StepsCounts_cw1 + StepsCounts_cb1 +
      psstot + phycomp + sapstot + sanstot + hba1c + agevisit + gender +
      (1 + dxgrp_factor | ID),
    sigma ~ Day + dxgrp_factor + StepsCounts_cw1 + StepsCounts_cb1 +
      sapstot + sanstot + (1 | ID)
  ),
  data         = df_long,
  init         = 0,
  cores        = 4,
  sample_prior = TRUE,
  seed         = 123
)
summary(fit_lsm_noday)

loo_compare(loo(fit_lsm), loo(fit_lsm_noday))
bayes_R2(fit_lsm)
bayes_R2(fit_lsm_noday)
plot(conditional_effects(fit_lsm, effects = "Day:dxgrp_factor"))
