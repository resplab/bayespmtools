# Tests for the plotting helpers in R/core.R: plot_cal_instability() and
# plot_cal_distance(). These are side-effect (base graphics) functions, so we
# just confirm they execute without error on a valid sample, drawing to a
# throwaway device.

cal_sample <- function() {
  data.frame(
    dist_type  = rep("beta", 3),
    dist_parm1 = c(1, 2, 3),
    dist_parm2 = c(3, 4, 5),
    cal_int    = c(0, 0.05, 0.1),
    cal_slp    = c(1, 0.9, 0.8),
    stringsAsFactors = FALSE
  )
}

with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

test_that("plot_cal_instability() runs without error (loess and line methods)", {
  set.seed(1)
  with_null_device({
    expect_no_error(plot_cal_instability(N = 200, sample = cal_sample()))
    expect_no_error(
      plot_cal_instability(N = 200, sample = cal_sample(), method = "line")
    )
  })
})

test_that("plot_cal_distance() runs without error (loess and line methods)", {
  set.seed(1)
  with_null_device({
    expect_no_error(plot_cal_distance(N = 200, sample = cal_sample()))
    expect_no_error(
      plot_cal_distance(N = 200, sample = cal_sample(), method = "line")
    )
  })
})
