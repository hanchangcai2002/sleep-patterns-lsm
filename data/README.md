# Data

The raw dataset (`SZ_obj_sleep_for_jinyuan.xlsx`) is **not included** in
this repository, due to patient privacy restrictions, and is excluded
via `.gitignore`.

Every script expects it at, relative to the repository root:

```
data/SZ_obj_sleep_for_jinyuan.xlsx
```

One row per subject (N = 63: 38 HC, 25 SZ), with columns including
`subnum`, `dxgrp`, `agevisit`, `gender`, `psstot`, `phycomp`, `sapstot`,
`sanstot`, `hba1c`, and per-day actigraphy variables (`TST.1`-`TST.12`,
`StepsCounts.1`-`StepsCounts.12`, etc., up to 12 days of wear).

## Synthetic data

`synthetic_sleep_data.csv` is a synthetic stand-in for the raw dataset
(N = 63: 38 HC, 25 SZ), generated with the eCDF-copula method of
Cai et al. (<https://doi.org/10.64898/2026.08.03.742474>). It contains
`subnum` (synthetic IDs), `dxgrp`, the subject-level covariates used in the
models, and `TST.1`-`TST.7` / `StepsCounts.1`-`StepsCounts.7`. It does not
contain any real subject's record and is intended only for running the code;
results will not match the manuscript. To use it, replace the `read_xlsx()`
line in a script with:

```r
dat <- read.csv("data/synthetic_sleep_data.csv")
```
