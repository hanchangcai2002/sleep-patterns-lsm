# Simulation Study

Simulation study checking whether the location-scale model (LSM) used in
the paper can recover group-specific between-subject variances, and
comparing it with the linear mixed-effects model (LME).

## Motivation

In the paper, diagnosis enters the location model as a subject-level
random slope, `(1 + dxgrp_factor | ID)`, although diagnosis does not vary
within subjects. Under this parameterization the between-subject variances
are

- HC: σ²<sub>b,HC</sub> = τ<sub>0</sub>²
- SZ: σ²<sub>b,SZ</sub> = τ<sub>0</sub>² + τ<sub>1</sub>² + 2ρτ<sub>0</sub>τ<sub>1</sub>

where τ<sub>0</sub> = `sd(Intercept)`, τ<sub>1</sub> = `sd(dxgrp_factorSZ)` and
ρ = `cor(Intercept, dxgrp_factorSZ)`. Because every subject belongs to only
one group, the data inform the two group-specific variances but not τ<sub>1</sub>
and ρ separately. The simulation therefore checks whether the model
recovers the two group-specific variances, which are the quantities of
interest in Eq. 6. It also examines how τ<sub>1</sub> and ρ behave.

The two groups share identical fixed effects and differ only in variance.
This is the setting where the LSM is expected to add information that the
LME cannot provide.

## Data-generating model

The model is a covariate-free version of the paper's location-scale model
(`01_generate_data.R`). For subject *i* on day *k*:

```
y_ik          = β0 + β_day[k] + b_i + ε_ik
b_i           ~ N(0, σ²_b,g(i))                  between-subject, group-specific
ε_ik          ~ N(0, σ²_ik)
log σ_ik      = τ0 + τ2·x_i + ω_i                within-subject
ω_i           ~ N(0, σ²_ω)                       shared by both groups, as in (1 | ID)
```

Here x<sub>i</sub> = 1 for SZ and 0 for HC. In Eq. 6 notation,
log σ²<sub>b,i</sub> = α<sub>0</sub> + α<sub>1</sub>x<sub>i</sub>.

- **Fixed effects are identical in both groups.** There is no diagnosis
  main effect and no day × diagnosis interaction, in either the location
  or the scale model.
- **Groups differ only in variance:** between-subject variance through α<sub>1</sub>,
  and within-subject variance through τ<sub>2</sub>.
- **Covariates are omitted.** Step counts and the subject-level
  covariates are left out, since they are not needed to answer the
  random-effects question.
- **Data are complete:** 7 days per subject, with no missing values.

## True parameter values

Values are calibrated to the 7-day imputed case-study fit and rounded
(`00_settings.R`).

| Parameter | Value | Source |
|---|---|---|
| β<sub>0</sub>, HC mean TST on day 1 | 370 min | observed median TST, 368 min |
| β<sub>day</sub>, day effects | 0 on all days | no strong day trend in the data |
| σ<sub>b,HC</sub>, HC between-subject SD | 60 min | `sd(Intercept)` = 60.6 |
| exp(τ<sub>0</sub>), HC within-subject SD | 60 min | `sigma_Intercept` = 4.08 → 59 min |
| σ<sub>ω</sub>, SD of the random scale effect | 0.4 | `sd(sigma_Intercept)` = 0.43 |
| SZ share of subjects | 25 / 63 | case study: 38 HC, 25 SZ |

The SZ parameters follow from the variance ratios in each scenario:

- σ<sub>b,SZ</sub> = σ<sub>b,HC</sub> · √(between-subject variance ratio), so α<sub>1</sub> = log(ratio)
- σ<sub>ε,SZ</sub> = σ<sub>ε,HC</sub> · √(within-subject variance ratio), so τ<sub>2</sub> = log(ratio) / 2

A variance ratio of 10 corresponds to an SD ratio of about 3.16. For
example, σ<sub>b,SZ</sub> ≈ 190 min when σ<sub>b,HC</sub> = 60 min. For comparison,
the case study implies a between-subject SD ratio of about 2.2, which is a
variance ratio of about 5.

## Scenarios

The design is fully crossed, with 12 scenarios:

| Factor | Levels |
|---|---|
| N, number of subjects | 63 (case-study size), 150, 300 |
| Between-subject variance ratio, SZ / HC | 1 (null), 10 |
| Within-subject variance ratio, SZ / HC | 1 (null), 10 |
| T, days per subject | 7 |

| Scenario | N | Between ratio | Within ratio |
|---|---|---|---|
| 1-3 | 63 / 150 / 300 | 1 | 1 |
| 4-6 | 63 / 150 / 300 | 10 | 1 |
| 7-9 | 63 / 150 / 300 | 1 | 10 |
| 10-12 | 63 / 150 / 300 | 10 | 10 |

Scenarios 1-3 are the full null, with no group difference of any kind.
Scenarios 4-6 match the pattern seen in the case study, where the group
difference is mainly between subjects. The design runs 100 replicates per
scenario.

## Models fitted to each dataset

Model fitting is in `02_fit_models.R`.

1. **LME, random intercept only.** This is the conventional homoscedastic
   model: `nlme::lme(TST ~ Day * dxgrp_factor, random = ~ 1 | ID)`.
2. **LME as in the paper.** The same model with `random = ~ dxgrp_factor | ID`
   (Supplemental Tab. 2).
3. **LSM as in the paper.** Fitted with `brms`, using the 4-day model
   structure extended to 7 days, without covariates:
   - location: `TST ~ 1 + Day * dxgrp_factor + (1 + dxgrp_factor | ID)`
   - scale: `sigma ~ Day * dxgrp_factor + (1 | ID)`

Day is a categorical factor in all three models. The fitted models include
the day × diagnosis terms even though their true values are zero, so the
fitted structure matches the paper.

**LSM priors.** The priors are fixed at the values `brms` selects by
default on the case-study TST (median 368 min, MAD 112.7 min):

- `student_t(3, 368, 112.7)` on the intercept
- `student_t(3, 0, 112.7)` on the location SDs
- `lkj(1)` on the correlation
- `student_t(3, 0, 2.5)` on the scale intercept and scale SD

Fixing the priors keeps the Stan code identical across datasets. The model
is therefore compiled once and reused for every replicate through
`update()`. `sample_prior = TRUE` is used, as in the paper, so the prior
and posterior of τ<sub>1</sub> and ρ can be compared.

**MCMC.** 4 chains with 2000 iterations each, of which 1000 are warmup, and
`init = 0`.

## Output

Each replicate is saved to `results/scenario_XX/rep_XXX.rds`, a list with
the following elements:

| Element | Content |
|---|---|
| `scenario`, `rep`, `seed` | scenario settings and random seed |
| `truth` | true parameter values for the scenario |
| `data` | the simulated dataset |
| `lme_ri`, `lme_paper` | fixed-effect estimates, SEs and p-values; between-subject SDs (the SZ SD is derived for `lme_paper`); residual SD |
| `lsm$params` | posterior summaries of all fixed effects and random-effect hyperparameters: mean, median, SD, 95% CrI, Rhat, ESS |
| `lsm$derived` | posterior summaries of σ<sub>b,HC</sub>, σ<sub>b,SZ</sub>, α<sub>1</sub>, τ<sub>0</sub>, τ<sub>2</sub> and σ<sub>ω</sub> |
| `lsm$prior` | prior-draw quantiles for the hyperparameters |
| `lsm$diagnostics` | divergent transitions, max-treedepth hits, max Rhat, number of parameters with Rhat > 1.01, minimum ESS ratio |

Each model element also records `ok` (whether the fit ran without error),
the error message if it did not, and the fitting time. Full `brmsfit`
objects are not kept, to save disk space. `results/scenarios.csv` lists
the scenario grid.

## Running

Run from the repository root:

```bash
Rscript simulation/run_simulation.R
```

To run a subset, for example as a pilot:

```bash
Rscript simulation/run_simulation.R --scenarios 1,4 --reps 1-2
```

- **Parallelism:** replicates run in parallel with `parallel::mclapply`,
  using `n_workers` workers (10 by default) and one core per replicate.
  The 4 chains within a replicate run sequentially.
- **Resuming:** replicates already saved are skipped, so an interrupted
  run can be restarted with the same command.
- **Seeds:** each replicate's seed is `seed_base + 1000 × scenario + rep`,
  so any single replicate can be reproduced on its own.

All settings are in `00_settings.R`: the scenario grid, true values,
number of replicates, MCMC settings and number of workers.

## Requirements

R with `brms` (either the `rstan` or `cmdstanr` backend), `posterior` and
`nlme`. The backend can be set with `options(brms.backend = "cmdstanr")`.

## Status

- Data generation and both LME fits have been checked.
- The `brms` fitting code has not been run yet.
- Performance evaluation is not written yet. Planned measures are bias,
  RMSE and coverage of the variance parameters, power and type I error
  for group differences in variance, fixed-effect coverage for the LSM
  versus the LME, and convergence rates.
