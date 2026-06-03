# Shared helpers / constants for the test suite.

# Width multiplier for a two-sided 95% normal CI: 2 * z_{0.975}.
K95 <- 2 * stats::qnorm(0.975)

# A small, well-behaved evidence list (formula form) used by integration tests.
# Mirrors the example evidence shipped in the package documentation.
example_evidence <- function() {
  list(
    prev     ~ beta(116, 155),
    cstat    ~ beta(3628, 1139),
    cal_mean ~ norm(-0.009, 0.125),
    cal_slp  ~ norm(0.995, 0.024)
  )
}

# Run an expression with all of testthat's noisy console output (messages,
# warnings, and cat()-based chatter from cobs/quantreg) silenced.
quietly <- function(expr) {
  out <- NULL
  invisible(utils::capture.output(
    suppressWarnings(suppressMessages(
      out <- force(expr)
    ))
  ))
  out
}
