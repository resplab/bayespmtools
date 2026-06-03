# Tests for calc_riley_vars() and riley_samp() in R/core.R.
# These implement Riley's frequentist variance / sample-size formulae, several
# of which have exact closed forms.

test_that("calc_riley_vars() returns the expected structure and signs", {
  parms <- list(prev = 0.5, cstat = 0.75, dist_type = "beta",
                dist_parm1 = 1, dist_parm2 = 1, cal_int = 0, cal_slp = 1)
  v <- bayespmtools:::calc_riley_vars(N = 100, parms = parms)

  expect_named(v, c("prev", "cstat", "cal_mean", "cal_oe",
                    "cal_int", "cal_slp", "cov_int_slp"))
  # All variances must be strictly positive.
  for (nm in c("prev", "cstat", "cal_mean", "cal_oe", "cal_int", "cal_slp")) {
    expect_true(v[[nm]] > 0, info = nm)
  }
})

test_that("calc_riley_vars() uses the exact binomial variance for prevalence", {
  parms <- list(prev = 0.3, cstat = 0.75, dist_type = "beta",
                dist_parm1 = 1, dist_parm2 = 1, cal_int = 0, cal_slp = 1)
  v100 <- bayespmtools:::calc_riley_vars(N = 100, parms = parms)
  v400 <- bayespmtools:::calc_riley_vars(N = 400, parms = parms)

  expect_equal(v100$prev, 0.3 * 0.7 / 100)
  expect_equal(v400$prev, 0.3 * 0.7 / 400)
  # The mean-calibration variance is defined to equal the prevalence variance.
  expect_equal(v100$cal_mean, v100$prev)
  # Variances shrink with sample size.
  expect_true(v400$cstat < v100$cstat)
  expect_true(v400$cal_int < v100$cal_int)
})

test_that("riley_samp() reproduces the closed-form prevalence sample size", {
  parms <- list(prev = 0.5, cstat = 0.75, dist_type = "beta",
                dist_parm1 = 1, dist_parm2 = 1, cal_int = 0, cal_slp = 1)
  out <- bayespmtools:::riley_samp(list(prev = 0.05), parms = parms)

  v <- (0.05 / K95)^2
  expect_equal(out$fciw.prev, round(0.5 * 0.5 / v))
})

test_that("riley_samp() only returns the requested targets and is monotone", {
  parms <- list(prev = 0.5, cstat = 0.75, dist_type = "beta",
                dist_parm1 = 1, dist_parm2 = 1, cal_int = 0, cal_slp = 1)

  out <- bayespmtools:::riley_samp(list(cstat = 0.05), parms = parms)
  expect_named(out, "fciw.cstat")

  # A wider tolerated interval requires a smaller sample size.
  tight <- bayespmtools:::riley_samp(list(cstat = 0.05), parms = parms)
  loose <- bayespmtools:::riley_samp(list(cstat = 0.10), parms = parms)
  expect_true(loose$fciw.cstat < tight$fciw.cstat)
})
