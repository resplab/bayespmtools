# Tests for find_n_mean() and find_n_quantile() in R/core.R, which invert a
# fitted (smoothed) curve of CI width versus sample size to find the N that
# achieves a target width. We feed them an exact convex/decreasing curve,
# ciw(N) = 1 / sqrt(N), whose inverse for a target t is N = (1 / t)^2.

test_that("find_n_mean() inverts a known decreasing convex curve", {
  N <- c(100, 200, 300, 400, 500)
  ciws <- matrix(1 / sqrt(N), nrow = 1)

  n_hat <- quietly(bayespmtools:::find_n_mean(target = 0.06, N = N, ciws = ciws))
  expect_true(is.finite(n_hat))
  expect_true(n_hat >= min(N) && n_hat <= max(N))
  # Exact answer is (1 / 0.06)^2 = 278; allow slack for the spline smoother.
  expect_equal(n_hat, (1 / 0.06)^2, tolerance = 0.1)
})

test_that("find_n_mean() needs a larger N for a tighter target", {
  N <- c(100, 200, 300, 400, 500)
  ciws <- matrix(1 / sqrt(N), nrow = 1)

  loose <- quietly(bayespmtools:::find_n_mean(target = 0.08, N = N, ciws = ciws))
  tight <- quietly(bayespmtools:::find_n_mean(target = 0.05, N = N, ciws = ciws))
  expect_true(tight > loose)
})

test_that("find_n_quantile() inverts the median of a noisy decreasing curve", {
  N <- c(100, 200, 300, 400, 500)
  base <- 1 / sqrt(N)
  ciws <- rbind(base, base * 1.02, base * 0.98)

  n_hat <- quietly(bayespmtools:::find_n_quantile(target = 0.06, N = N,
                                                  q = 0.5, ciws = ciws))
  expect_true(is.finite(n_hat))
  expect_true(n_hat >= min(N) && n_hat <= max(N))
  expect_equal(n_hat, (1 / 0.06)^2, tolerance = 0.15)
})
