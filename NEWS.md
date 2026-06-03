# bayespmtools 0.0.2

Bug-fix release addressing several issues in the evidence-processing and
correlation-induction internals.

* `process_evidence()`: named native beta parameters (e.g.
  `beta(alpha = 2, beta = 3)`) no longer error with "non-numeric argument to
  binary operator"; parameters supplied by name are now coerced to a numeric
  vector before use.
* `process_evidence()`: native `probitnorm` parameters now produce correct
  moments (an internal call referenced a mistyped distribution name and
  returned `NULL` moments).
* `process_evidence()`: named native parameters `mu`/`sigma` are now accepted
  for `logitnorm` and `probitnorm`; a mix of named and unnamed parameters
  (e.g. `beta(0.4, var = 0.04)`) is now rejected as ambiguous rather than
  silently misread; moment names are standardised to `m`/`v` across all
  specification styles.
* `infer_correlation()`: fixed matrix indexing that caused the prevalence and
  mean-calibration columns to be computed from a single value rather than the
  simulated data, which corrupted the correlation matrix used for correlation
  induction in `bpm_valsamp()` and `bpm_valprec()`.
* Calibration specified via the observed-to-expected ratio (`cal_oe`) now works
  in both `bpm_valsamp()` and `bpm_valprec()`. Previously `bpm_valsamp()`
  handled only mean calibration when imputing the calibration intercept, and
  correlation induction failed for `cal_oe` evidence because
  `infer_correlation()` did not return an O/E column. All three calibration
  parameterisations (`cal_mean`, `cal_oe`, `cal_int`, each with `cal_slp`) are
  now supported by both functions.
* Documentation: expanded `process_evidence()` to document all supported
  evidence specification styles (native parameters, moments with a variance or
  standard deviation, and mean plus upper quantile).
* Added a `print()` method for `bpm_evidence` objects (returned by
  `process_evidence()`) that displays each component's distribution, native
  parameters, and implied mean and variance.
* `summary()` for `bpm_evidence` objects now also reports the 95% range
  (2.5th and 97.5th percentiles) of each evidence component.
* Added package URL and BugReports fields pointing to the GitHub repository.
* Added a `testthat` test suite covering the moment helpers, performance-metric
  calculators, evidence processing, and the main entry points.

# bayespmtools 0.0.1

* Initial CRAN submission.
