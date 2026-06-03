# Tests for process_evidence() and the summary/print methods in R/evidence.R.
# This is the main user-facing input layer, so we cover formula input, named
# list input, distribution defaulting, and the validation error paths.

test_that("process_evidence() parses formula input into a bpm_evidence object", {
  ev <- example_evidence()
  pe <- process_evidence(ev)

  expect_s3_class(pe, "bpm_evidence")
  expect_setequal(names(pe), c("prev", "cstat", "cal_mean", "cal_slp"))

  # Parameters are carried through verbatim from the formula.
  expect_equal(unname(pe$prev$parms), c(116, 155))
  expect_equal(pe$prev$type, "beta")
  expect_equal(pe$cal_mean$type, "norm")
  expect_equal(unname(pe$cal_slp$parms), c(0.995, 0.024))

  # Moments are computed: Beta(116,155) has mean 116/271.
  expect_equal(pe$prev$moments[[1]], 116 / 271, tolerance = 1e-6)
  expect_equal(pe$cal_mean$moments[[1]], -0.009)
  expect_equal(pe$cal_mean$moments[[2]], 0.125^2, tolerance = 1e-8)
})

test_that("process_evidence() accepts named-list input and defaults distributions", {
  ev <- list(
    prev  = list(type = "beta", mean = 0.38, sd = 0.2),
    cstat = list(mean = 0.7, sd = 0.05),       # no type -> defaults to beta
    cal_int = list(mean = 0.2, sd = 0.2),      # no type -> defaults to norm
    cal_slp = list(mean = 0.8, sd = 0.3)
  )
  # Defaulting emits informative messages.
  expect_message(process_evidence(ev), "cstat")
  pe <- suppressMessages(process_evidence(ev))

  expect_s3_class(pe, "bpm_evidence")
  expect_equal(pe$cstat$type, "beta")
  expect_equal(pe$cal_int$type, "norm")
  # Beta fit by method of moments reproduces the requested mean & variance.
  expect_equal(pe$cstat$moments[[1]], 0.7, tolerance = 1e-8)
  expect_equal(pe$cstat$moments[[2]], 0.05^2, tolerance = 1e-8)
  expect_equal(sum(pe$cstat$parms) > 0, TRUE)
})

test_that("process_evidence() recognises each calibration parameterisation", {
  base <- list(prev ~ beta(116, 155), cstat ~ beta(3628, 1139),
               cal_slp ~ norm(0.995, 0.024))

  pe_mean <- process_evidence(c(base, list(cal_mean ~ norm(-0.009, 0.125))))
  expect_true("cal_mean" %in% names(pe_mean))

  pe_oe <- suppressMessages(
    process_evidence(c(base, list(cal_oe ~ norm(1.0, 0.1))))
  )
  expect_true("cal_oe" %in% names(pe_oe))

  pe_int <- process_evidence(c(base, list(cal_int ~ norm(0.0, 0.1))))
  expect_true("cal_int" %in% names(pe_int))
})

test_that("process_evidence() errors on malformed evidence", {
  expect_error(process_evidence("not a list"), "list")

  # Missing required components.
  expect_error(
    process_evidence(list(cstat ~ beta(3, 1), cal_slp ~ norm(1, .1),
                          cal_mean ~ norm(0, .1))),
    "prev"
  )
  expect_error(
    process_evidence(list(prev ~ beta(1, 1), cal_slp ~ norm(1, .1),
                          cal_mean ~ norm(0, .1))),
    "cstat"
  )
  # Calibration slope present but no intercept/mean/oe.
  expect_error(
    process_evidence(list(prev ~ beta(1, 1), cstat ~ beta(3, 1))),
    "calibration"
  )
  # Duplicate parameter specification.
  expect_error(
    process_evidence(list(prev ~ beta(1, 1), prev ~ beta(2, 1),
                          cstat ~ beta(3, 1), cal_slp ~ norm(1, .1),
                          cal_mean ~ norm(0, .1))),
    "Duplicate"
  )
})

test_that("summary.bpm_evidence() builds a tidy one-row-per-component table", {
  pe <- process_evidence(example_evidence())
  s <- summary(pe)

  expect_s3_class(s, "summary.bpm_evidence")
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), length(pe))
  expect_setequal(
    colnames(s),
    c("component", "type", "parameter1", "parameter2", "mean", "variance",
      "q2.5", "q97.5")
  )
  expect_setequal(s$component, c("prev", "cstat", "cal_mean", "cal_slp"))
  # The mean column matches the per-component first moment.
  prev_row <- s[s$component == "prev", ]
  expect_equal(prev_row$mean, 116 / 271, tolerance = 1e-6)

  # The print method runs and returns its argument invisibly.
  expect_output(print(s), "Summary of BPM evidence")
  expect_invisible(print(s))
})

test_that("summary.bpm_evidence() reports a sensible 95% range", {
  pe <- process_evidence(example_evidence())
  s <- summary(pe)

  # The 2.5% / 97.5% quantiles bracket the mean for every component.
  expect_true(all(s$q2.5 < s$mean))
  expect_true(all(s$mean < s$q97.5))

  # For a normal component the range is the textbook mean +/- 1.96 sd.
  slp <- s[s$component == "cal_slp", ]
  expect_equal(slp$q2.5, 0.995 + qnorm(0.025) * 0.024, tolerance = 1e-6)
  expect_equal(slp$q97.5, 0.995 + qnorm(0.975) * 0.024, tolerance = 1e-6)

  # When a beta component is specified by its upper bound (cih), q97.5 recovers it.
  pe2 <- process_evidence(list(
    prev ~ beta(116, 155), cstat ~ beta(mean = 0.761, cih = 0.773),
    cal_mean ~ norm(-0.009, 0.125), cal_slp ~ norm(0.995, 0.024)
  ))
  s2 <- summary(pe2)
  expect_equal(s2[s2$component == "cstat", "q97.5"], 0.773, tolerance = 1e-3)
})

test_that("print.bpm_evidence shows native parameters and moments", {
  pe <- process_evidence(example_evidence())

  expect_output(print(pe), "Bayesian prediction-model evidence")
  expect_invisible(print(pe))

  out <- capture.output(print(pe))
  # Every component is listed.
  for (nm in names(pe)) {
    expect_true(any(grepl(nm, out, fixed = TRUE)), info = nm)
  }
  # Native parameter names are shown (beta shapes and normal mean/sd), not just
  # the generic parameter1/parameter2 used by summary().
  expect_true(any(grepl("shape1", out)))
  expect_true(any(grepl("sd", out)))
  # Moment columns are present.
  expect_true(any(grepl("mean", out)) && any(grepl("variance", out)))
})


# ---------------------------------------------------------------------------
# Equivalence matrix: the core "flexible characterization" contract.
#
# process_evidence() must accept each element specified as either native
# distribution parameters OR moments (with the second moment given as a
# variance or an SD). Every spec style that describes the *same* underlying
# distribution must therefore produce the *same* parms and moments.
# ---------------------------------------------------------------------------

# Process a single `prev` specification against a fixed, valid remainder and
# return the processed `prev` element.
proc_prev <- function(prev_spec) {
  base <- list(
    cstat    ~ beta(3628, 1139),
    cal_mean ~ norm(-0.009, 0.125),
    cal_slp  ~ norm(0.995, 0.024)
  )
  suppressMessages(process_evidence(c(prev_spec, base)))$prev
}

test_that("all beta spec styles agree (native params <-> moments)", {
  # Target: Beta(2, 3) -> mean 0.4, var 0.04, sd 0.2.
  specs <- list(
    "positional shapes" = list(prev ~ beta(2, 3)),
    "named alpha/beta"  = list(prev ~ beta(alpha = 2, beta = 3)),
    "named mean/var"    = list(prev ~ beta(mean = 0.4, var = 0.04)),
    "named m/v"         = list(prev ~ beta(m = 0.4, v = 0.04)),
    "named mean/sd"     = list(prev ~ beta(mean = 0.4, sd = 0.2)),
    "list mean/var"     = list(prev = list(type = "beta", mean = 0.4, var = 0.04))
  )

  for (label in names(specs)) {
    p <- proc_prev(specs[[label]])
    expect_equal(p$type, "beta", info = label)
    expect_equal(unname(p$parms), c(2, 3), tolerance = 1e-5, info = label)
    expect_equal(unname(p$moments), c(0.4, 0.04), tolerance = 1e-6, info = label)
    expect_equal(names(p$moments), c("m", "v"), info = label)
  }
})

test_that("all normal spec styles agree (native params <-> moments)", {
  # Target: Normal(mean 0.4, sd 0.1) -> var 0.01.
  specs <- list(
    "positional mean/sd" = list(prev ~ norm(0.4, 0.1)),
    "named mean/sd"      = list(prev ~ norm(mean = 0.4, sd = 0.1)),
    "named mean/var"     = list(prev ~ norm(mean = 0.4, var = 0.01)),
    "named m/v"          = list(prev ~ norm(m = 0.4, v = 0.01)),
    "named m/sd"         = list(prev ~ norm(m = 0.4, sd = 0.1)),
    "list mean/var"      = list(prev = list(type = "norm", mean = 0.4, var = 0.01))
  )

  for (label in names(specs)) {
    p <- proc_prev(specs[[label]])
    expect_equal(p$type, "norm", info = label)
    expect_equal(unname(p$parms), c(0.4, 0.1), tolerance = 1e-6, info = label)
    expect_equal(unname(p$moments), c(0.4, 0.01), tolerance = 1e-8, info = label)
    expect_equal(names(p$moments), c("m", "v"), info = label)
  }
})

test_that("all logit-normal spec styles agree on moments", {
  # logit-normal has no named-native alias, so we check positional native
  # against the two moment forms. Target mean 0.3, var 0.01.
  ml <- proc_prev(list(prev ~ logitnorm(mean = 0.3, var = 0.01)))
  ms <- proc_prev(list(prev ~ logitnorm(mean = 0.3, sd = 0.1)))

  expect_equal(unname(ml$moments), c(0.3, 0.01), tolerance = 1e-4)
  expect_equal(unname(ms$moments), c(0.3, 0.01), tolerance = 1e-4)
  expect_equal(ml$parms, ms$parms, tolerance = 1e-4)
  expect_equal(names(ml$moments), c("m", "v"))
})

test_that("mean + upper-quantile (cih) spec reproduces the requested quantile", {
  # Normal: mean 0, 97.5th percentile at 1.96 -> sd ~ 1.
  pn <- proc_prev(list(prev ~ norm(mean = 0, cih = 1.96)))
  expect_equal(unname(pn$parms), c(0, 1), tolerance = 1e-2)

  # Beta: solved shapes must put 0.975 mass below the supplied cih.
  pb <- proc_prev(list(prev ~ beta(mean = 0.3, cih = 0.5)))
  expect_equal(pb$moments[["m"]], 0.3, tolerance = 1e-4)
  expect_equal(stats::pbeta(0.5, pb$parms[[1]], pb$parms[[2]]), 0.975,
               tolerance = 1e-3)
})

# ---------------------------------------------------------------------------
# Regression tests for two bugs found in process_evidence_element():
#   (1) named native beta params [beta(alpha=, beta=)] passed a list to
#       moments(), causing "non-numeric argument to binary operator".
#   (2) the probitnorm native branch called moments("probit", ...) -- a typo
#       for "probitnorm" -- yielding NULL moments.
# ---------------------------------------------------------------------------

test_that("named native beta parameters no longer error (regression)", {
  expect_no_error(proc_prev(list(prev ~ beta(alpha = 2, beta = 3))))
  p_named <- proc_prev(list(prev ~ beta(alpha = 2, beta = 3)))
  p_pos   <- proc_prev(list(prev ~ beta(2, 3)))
  expect_equal(p_named$parms, p_pos$parms)
  expect_equal(p_named$moments, p_pos$moments)
  expect_type(unclass(p_named$parms), "double")
})

test_that("probitnorm native parameters yield finite moments (regression)", {
  p <- proc_prev(list(prev ~ probitnorm(0, 1)))
  expect_equal(p$type, "probitnorm")
  expect_false(is.null(p$moments))
  expect_true(all(is.finite(p$moments)))
  # probitnorm(0, 1) is symmetric about 0.5 with variance 1/12.
  expect_equal(unname(p$moments), c(0.5, 1 / 12), tolerance = 1e-4)
  expect_equal(names(p$moments), c("m", "v"))
})

test_that("named native mu/sigma works for logit- and probit-normal", {
  # Previously only positional native parameters were accepted for these.
  l_named <- proc_prev(list(prev ~ logitnorm(mu = 0, sigma = 1)))
  l_pos   <- proc_prev(list(prev ~ logitnorm(0, 1)))
  expect_equal(l_named$parms, l_pos$parms)
  expect_equal(l_named$moments, l_pos$moments)

  p_named <- proc_prev(list(prev ~ probitnorm(mu = 0, sigma = 1)))
  p_pos   <- proc_prev(list(prev ~ probitnorm(0, 1)))
  expect_equal(p_named$parms, p_pos$parms)
  expect_equal(p_named$moments, p_pos$moments)
})

test_that("mixing named and unnamed parameters is rejected", {
  # beta(0.4, var = 0.04) is ambiguous: is 0.4 a shape or a mean?
  expect_error(
    proc_prev(list(prev ~ beta(0.4, var = 0.04))),
    "all named or all unnamed"
  )
  expect_error(
    proc_prev(list(prev ~ norm(0, sd = 0.1))),
    "all named or all unnamed"
  )
})
