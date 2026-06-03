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
* win-builder: R-devel, R-release (4.6.0), and R-oldrelease (4.5.3)

## R CMD check results

0 errors | 0 warnings | 0 notes

All three win-builder checks returned Status: OK (no errors, warnings, or
notes).

## Reverse dependencies

There are no reverse dependencies.

Regards,
Mohsen Sadatsafavi
