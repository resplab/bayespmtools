# Integration tests for bpm_valprec() (R/valprec.R): given fixed sample sizes,
# compute precision (expected / quantile CI widths) and value-of-information
# quantities. These are stochastic and somewhat slow, so we keep n_sim small,
# skip on CRAN, and assert on structure plus qualitative monotonicity.

test_that("bpm_valprec() returns precision targets as a matrix indexed by N", {
  skip_on_cran()
  set.seed(42)
  res <- quietly(bpm_valprec(
    N = c(500, 1000),
    evidence = example_evidence(),
    targets = list(eciw.cstat = TRUE, qciw.cal_slp = 0.9),
    n_sim = 40
  ))

  expect_type(res, "list")
  expect_setequal(names(res), c("N", "evidence", "ciws", "results", "sample"))

  m <- res$results
  expect_true(is.matrix(m))
  expect_equal(rownames(m), c("500", "1000"))
  expect_true(all(c("eciw.cstat", "qciw.cal_slp") %in% colnames(m)))
  expect_true(all(m > 0))

  # Expected CI width must shrink as the validation sample size grows.
  expect_true(m["1000", "eciw.cstat"] < m["500", "eciw.cstat"])
})

test_that("bpm_valprec() computes net-benefit VoI quantities", {
  skip_on_cran()
  set.seed(42)
  res <- quietly(bpm_valprec(
    N = c(500, 1000),
    evidence = example_evidence(),
    targets = list(eciw.cstat = TRUE, voi.nb = TRUE, oa.nb = TRUE),
    threshold = 0.2,
    n_sim = 40
  ))

  m <- res$results
  expect_true(all(c("voi.evpi", "voi.evsi", "voi.nb", "oa.nb") %in% colnames(m)))
  # EVPI is a property of the prior, so it is constant across sample sizes.
  expect_equal(m["500", "voi.evpi"], m["1000", "voi.evpi"])
  # EVSI accrues with sample size and cannot exceed EVPI; the ratio is in [0, 1].
  expect_true(m["1000", "voi.evsi"] >= m["500", "voi.evsi"])
  expect_true(all(m[, "voi.nb"] >= 0 & m[, "voi.nb"] <= 1.0001))
})

test_that("bpm_valprec() errors when NB targets lack a threshold", {
  skip_on_cran()
  set.seed(1)
  expect_error(
    quietly(bpm_valprec(
      N = 500, evidence = example_evidence(),
      targets = list(voi.nb = TRUE), n_sim = 20
    )),
    "threshold"
  )
})
