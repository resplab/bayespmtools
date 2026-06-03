## Submission notes

This is a bug-fix update (0.0.2) of an existing CRAN package. It corrects
several defects in the internal evidence-processing and correlation-induction
routines:

* `process_evidence()` no longer errors when native beta parameters are
  supplied by name (e.g. `beta(alpha = 2, beta = 3)`).
* Native `probitnorm` parameters now yield correct moments (an internal call
  used a mistyped distribution name).
* `infer_correlation()` previously mis-indexed the simulated data, corrupting
  the prevalence and mean-calibration entries of the correlation matrix used
  for correlation induction; this is now fixed.

A full list of changes is in NEWS.md. A `testthat` test suite has been added.

## Test environments

* Local: Windows 11, R 4.6.0 -- R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes

Note: vignettes were not rebuilt locally as Pandoc was unavailable in the
local environment; they are unchanged from the previous release.

## Reverse dependencies

There are no reverse dependencies.

Regards,
Mohsen Sadatsafavi
