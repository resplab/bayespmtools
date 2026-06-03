# Integration tests for bpm_valsamp() (R/valsamp.R): solve for the sample size
# that meets precision / assurance targets via pre-posterior simulation and
# Robbins-Monro stochastic root finding. Stochastic and slow, so n_sim is kept
# small and the heavier net-benefit optimisation is exercised separately.

test_that("bpm_valsamp() returns an estimated sample size for a precision target", {
  skip_on_cran()
  set.seed(42)
  res <- quietly(bpm_valsamp(
    evidence = example_evidence(),
    targets = list(eciw.cstat = 0.1),
    n_sim = 40
  ))

  expect_type(res, "list")
  expect_true(all(c("results", "sample", "evidence", "targets") %in% names(res)))
  expect_true("eciw.cstat" %in% names(res$results))

  n <- res$results[["eciw.cstat"]]
  expect_true(is.finite(n) && n > 0)
  expect_s3_class(res$evidence, "bpm_evidence")
})

test_that("bpm_valsamp() carries through the supplied targets and a sample", {
  skip_on_cran()
  set.seed(7)
  targets <- list(eciw.cal_slp = 0.3)
  res <- quietly(bpm_valsamp(
    evidence = example_evidence(),
    targets = targets,
    n_sim = 40
  ))

  expect_identical(res$targets, targets)
  expect_s3_class(res$sample, "data.frame")
  expect_true(nrow(res$sample) == 40)
  # The simulated c-statistic draws are kept inside the valid (0.51, 0.99) band.
  expect_true(all(res$sample$cstat >= 0.51 & res$sample$cstat <= 0.99))
})

test_that("bpm_valsamp() accepts a previously simulated sample as evidence", {
  skip_on_cran()
  set.seed(11)
  first <- quietly(bpm_valsamp(
    evidence = example_evidence(),
    targets = list(eciw.cstat = 0.1),
    n_sim = 40
  ))
  # Re-using the returned $sample data frame should not re-process evidence.
  res <- quietly(bpm_valsamp(
    evidence = first$sample,
    targets = list(eciw.cstat = 0.1)
  ))
  expect_true("eciw.cstat" %in% names(res$results))
})
