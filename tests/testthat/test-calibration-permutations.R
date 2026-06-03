# The package contract is that calibration can be specified by the slope plus
# exactly ONE of: mean calibration (cal_mean), the observed-to-expected ratio
# (cal_oe), or the calibration intercept (cal_int). Both entry points must
# accept all three parameterisations.
#
# Regression note: bpm_valsamp() previously hard-coded the cal_mean path (a
# `#TODO` at the intercept-imputation loop) and its base-parameter block read
# evidence$cal_mean unconditionally, so cal_oe evidence errored. bpm_valprec()
# accepted cal_oe in its own branches but then failed during correlation
# induction because infer_correlation() emitted no cal_oe column. These tests
# exercise every (function x calibration spec) combination.

cal_specs <- list(
  cal_mean = list(cal_mean ~ norm(-0.009, 0.125)),
  cal_oe   = list(cal_oe   ~ norm(1.00, 0.05)),
  cal_int  = list(cal_int  ~ norm(0.00, 0.10))
)

evidence_with <- function(cal_spec) {
  c(
    list(prev ~ beta(116, 155), cstat ~ beta(3628, 1139),
         cal_slp ~ norm(0.995, 0.024)),
    cal_spec
  )
}

for (spec_name in names(cal_specs)) {
  local({
    nm <- spec_name
    ev <- evidence_with(cal_specs[[nm]])

    test_that(paste0("bpm_valprec() works with ", nm, " calibration"), {
      skip_on_cran()
      set.seed(1)
      res <- quietly(bpm_valprec(
        N = c(500, 1000), evidence = ev,
        targets = list(eciw.cstat = TRUE), n_sim = 30
      ))
      expect_true(is.matrix(res$results))
      expect_true("eciw.cstat" %in% colnames(res$results))
      expect_true(all(res$results[, "eciw.cstat"] > 0))
      # The intercept is always present in the returned sample (supplied or imputed).
      expect_true("cal_int" %in% colnames(res$sample))
      expect_true(all(is.finite(res$sample$cal_int)))
    })

    test_that(paste0("bpm_valsamp() works with ", nm, " calibration"), {
      skip_on_cran()
      set.seed(1)
      res <- quietly(bpm_valsamp(
        evidence = ev,
        targets = list(eciw.cstat = 0.1), n_sim = 30
      ))
      n <- res$results[["eciw.cstat"]]
      expect_true(is.finite(n) && n > 0)
      expect_true("cal_int" %in% colnames(res$sample))
      expect_true(all(is.finite(res$sample$cal_int)))
    })
  })
}

test_that("correlation induction succeeds for cal_oe evidence (regression)", {
  skip_on_cran()
  set.seed(1)
  # impute_cor = TRUE (the default) is the path that previously broke for
  # cal_oe because infer_correlation() produced no cal_oe column.
  res <- quietly(bpm_valprec(
    N = 800, evidence = evidence_with(cal_specs$cal_oe),
    targets = list(eciw.cal_oe = TRUE), n_sim = 30, impute_cor = TRUE
  ))
  expect_true("eciw.cal_oe" %in% colnames(res$results))
  expect_true(all(res$results[, "eciw.cal_oe"] > 0))
})
