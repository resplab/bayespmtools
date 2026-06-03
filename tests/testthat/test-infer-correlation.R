# Tests for infer_correlation() in R/core.R, which simulates validation
# datasets and returns the correlation matrix of the performance metrics
# (prev, cstat, cal_mean, cal_oe, cal_int, cal_slp) used for correlation
# induction in bpm_valsamp()/bpm_valprec().
#
# Regression note: this function previously indexed the simulated data matrix
# with scalar subscripts -- mean(df[2]) / mean(df[2] - df[1]) -- instead of
# column subscripts, so the `prev` and `cal_mean` columns held a single
# arbitrary predicted probability rather than the simulated prevalence and
# mean calibration. The correlation matrix was therefore garbage. The cal_oe
# column was added so that evidence specified via the O/E ratio can be matched
# during correlation induction.

metric_names <- c("prev", "cstat", "cal_mean", "cal_oe", "cal_int", "cal_slp")

test_that("infer_correlation() returns a valid 6x6 correlation matrix", {
  skip_on_cran()
  set.seed(1)
  m <- bayespmtools:::infer_correlation("beta", c(2, 3),
                                        cal_int = 0, cal_slp = 1,
                                        n = 300, n_sim = 200)

  expect_true(is.matrix(m))
  expect_equal(dim(m), c(6L, 6L))
  expect_equal(rownames(m), metric_names)
  expect_equal(colnames(m), metric_names)

  # Hallmarks of a correlation matrix.
  expect_equal(diag(m), stats::setNames(rep(1, 6), metric_names))
  expect_true(isSymmetric(unname(m)))
  expect_true(all(m >= -1 - 1e-8 & m <= 1 + 1e-8))
})

test_that("infer_correlation() recovers the prev/cal_mean relationship (regression)", {
  skip_on_cran()
  set.seed(1)
  m <- bayespmtools:::infer_correlation("beta", c(2, 3),
                                        cal_int = 0, cal_slp = 1,
                                        n = 300, n_sim = 200)

  # cal_mean = prev - mean(predicted risk), so per-simulation prevalence and
  # mean calibration are strongly positively correlated. Under the old scalar
  # indexing bug this correlation collapsed to ~0.
  expect_true(m["prev", "cal_mean"] > 0.5)
})
