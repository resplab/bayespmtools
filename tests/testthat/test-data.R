# Tests for the bundled `isaric` example dataset (R/isaric.R, data/isaric.rda).

test_that("the isaric dataset loads with the documented shape", {
  isaric <- NULL
  data("isaric", package = "bayespmtools", envir = environment())

  expect_s3_class(isaric, "data.frame")
  expect_equal(nrow(isaric), 8L)
  expect_setequal(
    colnames(isaric),
    c("Region", "Sample_Size", "n", "n_events", "cstat", "cstat_l",
      "cal_mean", "cal_mean_l", "cal_slope", "cal_slope_l")
  )
})

test_that("isaric values are internally consistent", {
  isaric <- NULL
  data("isaric", package = "bayespmtools", envir = environment())

  # C-statistics live in (0.5, 1] and lower bounds sit below point estimates.
  expect_true(all(isaric$cstat > 0.5 & isaric$cstat <= 1))
  expect_true(all(isaric$cstat_l <= isaric$cstat))
  expect_true(all(isaric$cal_slope_l <= isaric$cal_slope))
  # Event and analysis counts are non-negative and ordered sensibly.
  expect_true(all(isaric$n_events >= 0))
  expect_true(all(isaric$n_events <= isaric$n))
})
