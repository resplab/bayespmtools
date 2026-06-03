# Tests for calc_cstat() and calc_se_sp() in R/core.R.
# Both have analytic values for simple inputs that we can pin down exactly.

test_that("calc_cstat() matches the closed form for a uniform risk distribution", {
  # Risks ~ Beta(1, 1) = Uniform(0, 1). The c-statistic integral evaluates to
  # (1 - 1/3 - 1/4) / (2 * 0.5 * 0.5) = 5/6.
  expect_equal(bayespmtools:::calc_cstat("beta", c(1, 1)), 5 / 6, tolerance = 1e-5)
})

test_that("calc_cstat() stays in (0.5, 1) and increases with spread", {
  c_lo <- bayespmtools:::calc_cstat("logitnorm", c(0, 0.5))
  c_mid <- bayespmtools:::calc_cstat("logitnorm", c(0, 1))
  c_hi <- bayespmtools:::calc_cstat("logitnorm", c(0, 2))

  expect_true(all(c(c_lo, c_mid, c_hi) > 0.5))
  expect_true(all(c(c_lo, c_mid, c_hi) < 1))
  # A wider spread of predicted risks means better discrimination.
  expect_true(c_lo < c_mid && c_mid < c_hi)
})

test_that("calc_se_sp() matches the closed form for the uniform/symmetric case", {
  # Beta(1,1), perfect calibration (int 0, slope 1), prevalence 0.5, threshold 0.5.
  # Threshold maps to hz = 0.5; TP = int_{0.5}^1 x dx = 0.375; se = sp = 0.75.
  res <- calc_se_sp("beta", c(1, 1), cal_int = 0, cal_slp = 1,
                    threshold = 0.5, prev = 0.5)
  expect_equal(unname(res[["se"]]), 0.75, tolerance = 1e-6)
  expect_equal(unname(res[["sp"]]), 0.75, tolerance = 1e-6)
})

test_that("calc_se_sp() always returns sensitivity and specificity within [0, 1]", {
  grid <- expand.grid(
    thr = c(0.05, 0.2, 0.5, 0.8),
    slp = c(0.7, 1, 1.3),
    int = c(-0.3, 0, 0.3)
  )
  for (i in seq_len(nrow(grid))) {
    res <- calc_se_sp("logitnorm", c(0, 1),
                      cal_int = grid$int[i], cal_slp = grid$slp[i],
                      threshold = grid$thr[i], prev = 0.5)
    expect_named(res, c("se", "sp"))
    expect_true(res[["se"]] >= 0 && res[["se"]] <= 1)
    expect_true(res[["sp"]] >= 0 && res[["sp"]] <= 1)
  }
})
