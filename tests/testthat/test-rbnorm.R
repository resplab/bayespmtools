# Tests for rbnorm() in R/core.R: a bivariate normal sampler parameterised by
# marginal means/variances and a covariance.

test_that("rbnorm() returns an n x 2 matrix with the documented column names", {
  set.seed(1)
  x <- bayespmtools:::rbnorm(n = 10, mu1 = 0, mu2 = 0, var1 = 1, var2 = 1, cov = 0.5)
  expect_true(is.matrix(x))
  expect_equal(dim(x), c(10L, 2L))
  expect_equal(colnames(x), c("x1", "x2"))
})

test_that("rbnorm() recovers the target moments at large n", {
  set.seed(123)
  x <- bayespmtools:::rbnorm(n = 2e5, mu1 = 1, mu2 = 2, var1 = 1, var2 = 4, cov = 1)

  expect_equal(colMeans(x), c(x1 = 1, x2 = 2), tolerance = 0.02)
  expect_equal(unname(apply(x, 2, var)), c(1, 4), tolerance = 0.05)
  expect_equal(stats::cov(x)[1, 2], 1, tolerance = 0.05)
})

test_that("rbnorm() produces uncorrelated columns when cov = 0", {
  set.seed(7)
  x <- bayespmtools:::rbnorm(n = 1e5, mu1 = 0, mu2 = 0, var1 = 1, var2 = 1, cov = 0)
  expect_equal(stats::cor(x)[1, 2], 0, tolerance = 0.02)
})
