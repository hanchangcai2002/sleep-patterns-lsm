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
