# Tests for the distribution moment helpers in R/core.R and R/evidence.R:
#   moments(), inv_moments(), inv_mean_quantile()
# These are deterministic, so we assert against closed-form values and verify
# that moments() and inv_moments() are inverses of each other.

test_that("moments() returns the textbook mean and variance for each family", {
  # Normal: mean = parms[1], variance = parms[2]^2
  expect_equal(unname(bayespmtools:::moments("norm", c(0, 1))), c(0, 1))
  expect_equal(unname(bayespmtools:::moments("norm", c(2, 3))), c(2, 9))

  # Beta(a, b): mean = a/(a+b), var = ab / ((a+b)^2 (a+b+1))
  m <- bayespmtools:::moments("beta", c(2, 3))
  expect_equal(m[[1]], 0.4)
  expect_equal(m[[2]], 2 * 3 / ((5^2) * 6))

  # logit-normal(0, 1) is symmetric about 0.5, so the mean is exactly 0.5.
  ml <- bayespmtools:::moments("logitnorm", c(0, 1))
  expect_equal(ml[[1]], 0.5, tolerance = 1e-6)
  expect_true(ml[[2]] > 0 && ml[[2]] < 0.25)

  # probit-normal(0, 1) is also symmetric about 0.5.
  mp <- bayespmtools:::moments("probitnorm", c(0, 1))
  expect_equal(mp[[1]], 0.5, tolerance = 1e-6)
  expect_equal(mp[[2]], 1 / 12, tolerance = 1e-4)
})

test_that("inv_moments() recovers the natural parameters (round trip)", {
  # Beta: method-of-moments should return the generating shapes.
  m <- bayespmtools:::moments("beta", c(2, 3))
  expect_equal(
    unname(bayespmtools:::inv_moments("beta", c(m[[1]], m[[2]]))),
    c(2, 3),
    tolerance = 1e-8
  )

  # Normal: mean unchanged, sd = sqrt(var).
  expect_equal(
    unname(bayespmtools:::inv_moments("norm", c(0.2, 0.04))),
    c(0.2, 0.2)
  )

  # logit-normal: optim-based inversion should round-trip within tolerance.
  ml <- bayespmtools:::moments("logitnorm", c(0.3, 0.8))
  back <- bayespmtools:::inv_moments("logitnorm", c(ml[[1]], ml[[2]]))
  expect_equal(unname(back), c(0.3, 0.8), tolerance = 1e-3)
})

test_that("inv_moments() rejects out-of-support means and variances", {
  expect_error(bayespmtools:::inv_moments("beta", c(1.5, 0.1)), "Mean must be")
  expect_error(bayespmtools:::inv_moments("beta", c(0.5, -0.1)), "Variance must be")
  expect_error(bayespmtools:::inv_moments("logitnorm", c(0, 0.1)), "Mean must be")
  expect_error(bayespmtools:::inv_moments("logitnorm", c(0.5, 0)), "Variance must be")
})

test_that("inv_mean_quantile() solves for parameters matching a target quantile", {
  # Normal with mean 0 and 97.5th percentile at 1.96 implies sd ~= 1.
  out <- bayespmtools:::inv_mean_quantile("norm", m = 0, q = 1.96, p = 0.975)
  expect_equal(out[["mean"]], 0)
  expect_equal(unname(out[2]), 1, tolerance = 1e-2)

  # Beta with mean 0.5 is symmetric, so shape1 == shape2.
  outb <- bayespmtools:::inv_mean_quantile("beta", m = 0.5, q = 0.7, p = 0.975)
  expect_equal(unname(outb[1]), unname(outb[2]), tolerance = 1e-4)
  # And the solved parameters should reproduce the requested quantile.
  expect_equal(stats::pbeta(0.7, outb[1], outb[2]), 0.975, tolerance = 1e-4)

  # logit-normal: solved parameters should reproduce the requested mean & quantile.
  outl <- bayespmtools:::inv_mean_quantile("logitnorm", m = 0.3, q = 0.5, p = 0.975)
  expect_equal(
    logitnorm::momentsLogitnorm(outl[1], outl[2])[[1]],
    0.3,
    tolerance = 1e-3
  )
})
