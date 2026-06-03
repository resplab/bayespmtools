# Tests for the pre-posterior CI-width engines in R/core.R:
#   calc_ciw_2s(), calc_ciw_sample(), calc_ciw_mc().
# These are stochastic, so we assert on structure, positivity, and the
# qualitative expectation that intervals narrow as the sample size grows.

ciw_metrics <- c("cstat", "cal_oe", "cal_mean", "cal_int", "cal_slp")

test_that("calc_ciw_2s() returns one CI width per requested sample size", {
  parms <- list(prev = 0.2, cstat = 0.75, dist_type = "beta",
                dist_parm1 = 1, dist_parm2 = 2, cal_int = 0, cal_slp = 1)
  set.seed(1)
  res <- bayespmtools:::calc_ciw_2s(N = c(100, 300, 500), parms = parms)

  expect_named(res, ciw_metrics)
  expect_true(all(vapply(res, length, integer(1)) == 3))
  expect_true(all(unlist(res) > 0))
})

test_that("calc_ciw_sample() returns positive widths that shrink with N", {
  parms <- list(prev = 0.2, dist_type = "beta", dist_parm1 = 2,
                dist_parm2 = 4, cal_int = 0, cal_slp = 1)
  set.seed(1)
  res <- bayespmtools:::calc_ciw_sample(N = c(100, 300, 1000), parms = parms)

  expect_named(res, ciw_metrics)
  expect_true(all(unlist(res) > 0))
  # Larger samples give narrower c-statistic intervals.
  expect_true(res$cstat[1] > res$cstat[3])
})

test_that("calc_ciw_mc() stacks per-draw results into matrices", {
  parms_sample <- data.frame(
    cstat = c(0.75, 0.8), prev = c(0.2, 0.3),
    dist_type = c("beta", "beta"),
    dist_parm1 = c(0.9, 1), dist_parm2 = c(1.5, 2),
    cal_int = c(0, 0.5), cal_slp = c(1, 0.9),
    stringsAsFactors = FALSE
  )
  set.seed(1)
  res <- bayespmtools:::calc_ciw_mc(N = c(100, 300, 500),
                                    parms_sample = parms_sample, method = "2s")

  expect_named(res, ciw_metrics)
  # One row per draw (2), one column per sample size (3).
  for (nm in ciw_metrics) {
    expect_equal(dim(res[[nm]]), c(2L, 3L), info = nm)
  }
})
