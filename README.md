# Between- and Within-subject Variability in Actigraphy Data: Analysis Code

R code accompanying:

> Liu J, Cai H, Zhang Z, Wang J, Lee E, Zhang X. *Between- and
> Within-subject Variability in Actigraphy Data: A Case Study on Sleep
> Patterns in Schizophrenia.* (preprint, bioRxiv 2026.704499)

## Repository structure

```
data/
  README.md            # expected data file and schema (not published)
  synthetic_sleep_data.csv  # synthetic stand-in (eCDF-copula)
four_day_model.R        # Section 4.2 -> Table 1
seven_day_model.R       # Section 4.3.3 -> Table 2
simulation/             # simulation study
```

Run either R script from the repository root. Each script is self-contained
apart from the protected input dataset described below.

The simulation study, which checks recovery of the group-specific
between-subject variances, lives in [`simulation/`](simulation/). See
[`simulation/README.md`](simulation/README.md) for the design and how to
run it. It uses simulated data only and does not need the protected
dataset.

## Data

Not published, due to patient privacy restrictions — see
[`data/README.md`](data/README.md) for the expected file and its schema.

A synthetic dataset, generated with the eCDF-copula method
(<https://doi.org/10.64898/2026.08.03.742474>), is provided in
[`data/synthetic_sleep_data.csv`](data/synthetic_sleep_data.csv) so the
scripts can be run end to end.

## Correspondence with the manuscript

**Both scripts**
- Missing subject-level covariates (stress, physical health, SAPS, SANS,
  HbA1c, age, gender) imputed with `mice` (random forest, `m = 5`); the
  first completed dataset is used
- Daily step counts person-mean centered into a between-subject part
  (`StepsCounts_cb1`) and a within-subject daily deviation
  (`StepsCounts_cw1`)

**`four_day_model.R`** (days 1-4, all complete; N = 63, 252 observations)
- Linear mixed-effects model, Eq. 2, Section 4.2.1
- Location-scale model — location Eq. 4, between-subject variance Eq. 6,
  scale Eq. 7, with day as a categorical factor — Table 1

**`seven_day_model.R`** (days 1-7)
- Complete-case data (N = 40, 280 observations) — Table 2 "Observed"
- TST and step counts multiply imputed (`mice`, random forest, `m = 5`,
  earlier days only used to impute later days; N = 63, 441 observations)
  — Table 2 "Imputed"
- Location-scale model — location Eq. 5, scale Eq. 8, day modeled with a
  group-specific thin-plate spline `s(Day, by = dxgrp_factor, k = 5)`

Both scripts also fit a comparison model with a common (not
group-specific) day effect, and report `loo_compare()` and `bayes_R2()`
for the two models.

## Model

Both scripts fit a mixed-effects location-scale model with `brms::brm()`:

- **Location** (mean submodel): day, diagnosis, their interaction, the
  within-/between-subject step-count effects, and the remaining
  subject-level covariates, with random effects `(1 + dxgrp_factor | ID)`.
- **Scale** (`sigma ~ ...`, within-subject scale submodel, modeled on the
  log standard-deviation scale): the same day/group/step-count terms, plus
  positive and negative symptom severity, with a subject-level random
  scale intercept `(1 | ID)` (ω<sub>i</sub> in Eqs. 7-8).

Between-subject variance (Eq. 6) is absorbed into the random location
term `(1 + dxgrp_factor | ID)` rather than fit as a separate parameter.

## Software

R, with `brms` (Stan backend), `nlme`, `mice`, `randomForest`,
`tidyverse`, and `readxl`.

## Reading the model output

`summary(fit)` reports location and scale coefficients in the
**Regression Coefficients** block. Rows without a `sigma_` prefix belong
to the location submodel; rows with a `sigma_` prefix belong to the
scale submodel. `l-95% CI`/`u-95% CI` give the 95% Bayesian credible
interval. `Rhat`, `Bulk_ESS`, and `Tail_ESS` are sampling diagnostics
and should be checked before interpreting a fitted model.

In the **Multilevel Hyperparameters** block, `sd(Intercept)` is the
between-subject SD for HC. Because `dxgrp_factorSZ` enters as a random
slope, `sd(dxgrp_factorSZ)` is not itself the SZ between-subject SD; that
is `sqrt(sd(Intercept)^2 + sd(dxgrp_factorSZ)^2 + 2 * cor * sd(Intercept) * sd(dxgrp_factorSZ))`.
In the 7-day models, the `sds(s...)` rows under **Smoothing Spline
Hyperparameters** control the wiggliness of the day splines and are not
between-subject variance components.
