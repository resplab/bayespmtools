# Tests for infer_cal_int_from_mean() and infer_cal_int_from_oe() in R/core.R.
# Under perfect calibration the recovered calibration intercept must be ~0,
# and the two routines should agree on the relationship between mean
# calibration and the O/E ratio.

test_that("a perfectly calibrated model yields a zero calibration intercept", {
  for (dist in list(list("beta", c(2, 3)), list("logitnorm", c(0, 1)))) {
    dt <- dist[[1]]; dp <- dist[[2]]
    prev <- bayespmtools:::moments(dt, dp)[[1]]

    # O/E = 1 (observed equals expected) <=> intercept 0 when slope = 1.
    int_oe <- bayespmtools:::infer_cal_int_from_oe(dt, dp, cal_oe = 1,
                                                   cal_slp = 1, prev = prev)
    expect_equal(int_oe, 0, tolerance = 1e-4)

    # Mean calibration error 0 <=> intercept 0 when slope = 1.
    int_mean <- bayespmtools:::infer_cal_int_from_mean(dt, dp, cal_mean = 0,
                                                       cal_slp = 1, prev = prev)
    expect_equal(int_mean, 0, tolerance = 1e-4)
  }
})

test_that("the inferred calibration intercept increases with the target O/E ratio", {
  # A larger O/E ratio (more observed events per predicted) requires a larger
  # calibration intercept, so the inferred intercept is monotone in O/E.
  dt <- "logitnorm"; dp <- c(0, 1)
  prev <- 0.5
  int_low_oe  <- bayespmtools:::infer_cal_int_from_oe(dt, dp, cal_oe = 0.8,
                                                      cal_slp = 1, prev = prev)
  int_high_oe <- bayespmtools:::infer_cal_int_from_oe(dt, dp, cal_oe = 1.2,
                                                      cal_slp = 1, prev = prev)
  expect_true(int_low_oe < int_high_oe)
})
