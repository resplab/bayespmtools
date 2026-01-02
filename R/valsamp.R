#' Bayesian Sample Size Calculator for External Validation
#'
#' @description
#' Bayesian sample size calculation for external validation studies of
#' clinical risk prediction models. The function evaluates sample sizes
#' required to meet precision-, assurance-, or decision-based targets
#' using pre-posterior simulation.
#'
#' @param evidence
#' A named list containing prior evidence components for model performance
#' parameters (e.g., prevalence, discrimination, calibration).
#' Alternatively, `evidence` may be a data frame of pre-posterior
#' draws (element \code{$sample}) returned by a previous call to this
#' function or to \code{bpm_valprec()}, in which case those draws are used
#' directly.
#'
#' @param targets
#' A named list specifying sample size targets.
#'
#' Supported targets include:
#' \itemize{
#'   \item Precision-based targets using expected 95\% interval widths
#'   (prefix \code{eciw}).
#'   \item Assurance-based targets specifying the probability that the
#'   95\% interval width does not exceed a given value
#'   (prefix \code{qciw}).
#'   \item Net benefit targets, including optimality assurance
#'   (\code{oa.nb}) and value-of-information ratios
#'   (\code{voi.nb = EVSI / EVPI}).
#' }
#'
#' For example, \code{eciw.cstat = 0.1} targets an expected interval width
#' of 0.1 for the c-statistic, while
#' \code{qciw.cal_slp = c(0.90, 0.22)} targets a 90 percent assurance that the
#' calibration slope interval width does not exceed 0.22.
#' Finally, \code{oa.nb = 0.80} targets a sample size that would correspond to
#' 80 percent assurance that the strategy with the highest NB in the sample will
#' be the strategy with the highest NB in the  population.
#'
#' @param n_sim
#' Number of Monte Carlo simulations used to generate the pre-posterior
#' distribution. If evidence is a data frame from previous calls to relevant functions,
#' n_sim will automatically be set to the number of rows of the data frame.
#'
#' @param method
#' Method used to compute the pre-posterior distribution of 95\% intervals.
#' One of \code{"sample"} (simulation-based) or \code{"2s"} (two-stage
#' approximation). Default is \code{"sample"}.
#'
#' @param threshold
#' Risk threshold used for decision-analytic quantities and net benefit
#' calculations. Required if \code{oa.nb} or \code{voi.nb} targets are
#' specified.
#'
#' @param dist_type
#' Distribution assumed for calibrated risks. Default is
#' \code{"logitnorm"}.
#'
#' @param impute_cor
#' Logical indicating whether correlation between performance measures
#' should be induced when simulating from marginal evidence distributions.
#' Default is \code{TRUE}.
#'
#' @param ex_args
#' Optional list of additional arguments passed to internal simulation or
#' root-finding routines (experimental feature).
#'
#' @return
#' A list with the following components:
#' \itemize{
#'   \item \code{results}: Estimated sample sizes required to meet each target.
#'   \item \code{sample}: Data frame of pre-posterior simulation draws.
#'   \item \code{evidence}: Processed evidence object used in the analysis.
#'   \item \code{trace}: Trace output from the stochastic root-finding algorithm.
#'   \item \code{targets}: The targets argument supplied to the function.
#' }
#'
#' @examples
#' \donttest{
#' evidence <- list(
#'  prev~beta(mean=0.428, sd=0.030),
#'  cstat~beta(mean=0.761, cih=0.773),
#'  cal_mean~norm(-0.009, 0.125),  #mean and SD
#'  cal_slp~norm(0.995, 0.024)     #mean and SD
#' )
#'
#' targets <- list(
#'   eciw.cstat = 0.1,
#'   qciw.cstat = c(0.9, 0.1),
#'   oa.nb      = 0.8
#' )
#'
#' samp <- bpm_valsamp(
#'   evidence  = evidence,
#'   targets   = targets,
#'   n_sim     = 1000,
#'   threshold = 0.2
#' )
#'
#' samp$results
#' }
#'
#' @export
bpm_valsamp <- function(
  evidence,
  targets,
  n_sim = NULL,
  method = "sample",
  threshold = NULL,
  dist_type = "logitnorm",
  impute_cor = TRUE,
  ex_args = NULL
) {
  out <- list()

  tmp <- as.data.frame(lapply(names(targets), strsplit, "[.]"))
  target_rules <- unname(unlist(tmp[1, ]))
  target_metrics <- unname(unlist(tmp[2, ]))
  target_values <- (targets)

  #Process requested stuff. Mainly to remove those if the rule=F
  to_remove <- c()
  for (i in 1:length(target_values)) {
    if (isFALSE(target_values[[i]])) {
      to_remove <- c(to_remove, i)
    }
  }
  if (length(to_remove) > 0) {
    target_values <- target_values[-to_remove]
    target_metrics <- target_metrics[-to_remove]
    target_rules <- target_rules[-to_remove]
  }

  #Check for some basic validation
  for (i in 1:length(target_values)) {
    if (target_rules[i] == "eciw") {
      if (
        !target_metrics[i] %in%
          c("cstat", "cal_mean", "cal_int", "cal_slp", "cal_oe")
      ) {
        stop(paste0("Target metric ", target_metrics[i], " not recognized."))
      }
      if (length(target_values[[i]]) != 1) {
        stop(paste0(
          "For eciw rules, I expect a single (scalar) values. This is not the case for ",
          names(targets)[i],
          "."
        ))
      }
    }
    if (target_rules[i] == "qciw") {
      if (
        !target_metrics[i] %in%
          c("cstat", "cal_mean", "cal_int", "cal_slp", "cal_oe")
      ) {
        stop(paste0("Target metric ", target_metrics[i], " not recognized."))
      }
      if (length(target_values[[i]]) != 2) {
        stop(paste0(
          "For qciw rules, I expect a vector of size 2 including the quantile (assurance) and desired CI width (in the same order, e.g. c(0.9, 0.1)). This is not the case for ",
          names(targets)[i],
          "."
        ))
      }
    }
  }

  if (is.function(ex_args$f_progress)) {
    f_progress <- ex_args$f_progress
  } else {
    f_progress <- base::message
  }

  ##Step 1: Process evidence
  f_progress("Processing evidence...")
  evidence <- process_evidence(evidence)
  out$evidence <- evidence

  ##Generate base parms
  base <- list()
  base$dist_type <- dist_type
  base$prev <- evidence$prev$moments[[1]]
  base$cstat <- evidence$cstat$moments[[1]]
  tmp <- mcmapper::mcmap(c(base$prev, base$cstat), type = dist_type)$valu
  base$dist_parm1 <- tmp[1]
  base$dist_parm2 <- tmp[2]
  base$cal_slp <- evidence$cal_slp$moments[[1]]
  #Cal intercept is not provided. Need to derive it for base params for correlation induction
  if (is.na(match("cal_int", names(evidence)))) {
    base$cal_int <- infer_cal_int_from_mean(
      base$dist_type,
      c(base$dist_parm1, base$dist_parm2),
      cal_mean = evidence$cal_mean$moments[[1]],
      cal_slp = base$cal_slp,
      prev = base$prev
    )
  } else {
    base$cal_int <- evidence$cal_int$moments[[1]]
  }

  #Step 2: generate sample of marginals
  f_progress("Generating Monte Carlo sample...")
  sample <- NULL
  for (element in evidence) {
    sample <- cbind(
      sample,
      do.call(
        paste0("r", element$type),
        args = as.list(c(n = n_sim, element$parms))
      )
    )
  }
  colnames(sample) <- names(evidence)

  n_bads <- 0
  repeat {
    bads <- which(sample[, 'cstat'] < 0.51 | sample[, 'cstat'] > 0.99)
    if (length(bads) == 0) {
      break
    }
    n_bads <- n_bads + length(bads)
    subsample <- NULL
    for (element in evidence) {
      subsample <- cbind(
        subsample,
        do.call(
          paste0("r", element$type),
          args = as.list(c(n = length(bads), element$parms))
        )
      )
    }
    sample[bads, ] <- subsample
  }
  if (n_bads > 0) {
    warning(paste(
      "in step 'Generating MOnte Carlo sample' - ",
      n_bads,
      "observations were replaced due to bad value of c-statistic."
    ))
  }

  #Step 3: induce correlation (if asked)
  if (impute_cor) {
    eff_n <- round(
      evidence$prev$moments[[1]] *
        (1 - evidence$prev$moments[[1]]) /
        evidence$prev$moments[[2]],
      0
    )
    f_progress(paste(
      "Imputing correlation, based on effective sample size:",
      eff_n,
      "..."
    ))
    base_cor <- infer_correlation(
      base$dist_type,
      c(base$dist_parm1, base$dist_parm2),
      base$cal_int,
      base$cal_slp,
      eff_n,
      1000
    )
    good_rows <- match(colnames(sample), colnames(base_cor))
    base_cor <- base_cor[good_rows, good_rows]
    sample <- mc2d::cornode(sample, target = base_cor)
  }

  #Step 4: if intercept is missing, impute it for the whole sample
  f_progress("Infering calibration intercept...")
  sample <- as.data.frame(sample)
  sample$dist_type <- dist_type
  sample$dist_parm1 <- 0
  sample$dist_parm2 <- 0

  if (is.na(match("cal_int", names(evidence)))) {
    sample$cal_int <- NA
    for (i in 1:nrow(sample)) {
      prev <- unname(sample[i, 'prev'])
      cstat <- unname(sample[i, 'cstat'])
      cal_mean <- unname(sample[i, 'cal_mean']) #TODO
      cal_slp <- unname(sample[i, 'cal_slp'])

      parms <- mcmap(c(prev, cstat), dist_type)$value
      sample$dist_parm1[i] <- parms[1]
      sample$dist_parm2[i] <- parms[2]

      cal_int <- infer_cal_int_from_mean(
        dist_type = dist_type,
        dist_parms = parms,
        cal_mean = cal_mean,
        cal_slp = cal_slp,
        prev = prev
      )

      sample[i, 'cal_int'] <- cal_int
    }
  }

  # Step 5: Bayesian Riley
  f_progress("Computing CI sample size...")

  if (!is.na(match("fciw", target_rules))) {
    #Frequentist CIWs
    target_ciws <- list()
    for (i in which(target_rules == 'fciw')) {
      target_ciws[target_metrics[[i]]] <- target_values[[i]]
    }
    N <- riley_samp(target_ciws, parms = base)
    out$N <- N
  }

  if (
    !is.na(match("eciw", target_rules)) | !is.na(match("qciw", target_rules))
  ) {
    indices <- which(target_rules == "eciw" | target_rules == "qciw")
    res <- find_n_RM(
      sample,
      target_rules[indices],
      target_metrics[indices],
      target_values[indices],
      base_parms = base
    )
    out$N <- c(out$N, res$N)
    out$trace <- res$trace
  }

  # Step 6: Calculate se and sp
  b_voi <- !is.na(match("voi", target_rules)) & !isTRUE(targets$voi.nb)
  b_oa <- !is.na(match("oa", target_rules)) &
    !isTRUE(targets$oa.nb)
  if (b_voi | b_oa) {
    if (is.null(threshold)) {
      stop("NB-related stuff was requested but threshold is not specified")
    } #Todo: move earlier to avoid this late error
    f_progress("Computing se/sp...")
    sample$sp <- sample$se <- NA
    for (i in 1:nrow(sample)) {
      se_sp <- calc_se_sp(
        sample$dist_type[i],
        c(sample$dist_parm1[i], sample$dist_parm2[i]),
        sample$cal_int[i],
        sample$cal_slp[i],
        threshold,
        sample$prev[i]
      )
      sample[i, c('se', 'sp')] <- se_sp
    }

    f_progress("VoI / NB assuraance...")

    #require(evsiexval)
    #res <- evsiexval::EVSI_gf(sample[,c('prev','se','sp')], future_sample_sizes=N,  ignore_prior=TRUE, z=threshold)

    tnb1 <- sample$prev *
      sample$se -
      (1 - sample$prev) * (1 - sample$sp) * threshold / (1 - threshold)
    tnb2 <- sample$prev - (1 - sample$prev) * threshold / (1 - threshold)
    tnbs <- cbind(0, tnb1, tnb2)
    maxnbs <- apply(tnbs, 1, max)
    maxenb <- max(colMeans(tnbs))
    oa0 <- mean(apply(tnbs, 1, which.max) == which.max(colMeans(tnbs)))
    evpi <- mean(maxnbs) - maxenb

    if (is.null(n_sim)) {
      n_sim <- nrow(sample)
    }

    f <- function(x, voi) {
      n <- c(round(x))

      nd <- rbinom(rep(n, n_sim), n, sample$prev)
      ntp <- rbinom(rep(n, n_sim), nd, sample$se)
      nfp <- rbinom(rep(n, n_sim), n - nd, 1 - sample$sp)

      nb1 <- ntp / n - nfp / n * threshold / (1 - threshold)
      nb2 <- nd / n - (1 - nd / n) * threshold / (1 - threshold)

      winners <- apply(cbind(0, nb1, nb2), 1, which.max)
      winnertnbs <- tnbs[cbind(1:n_sim, winners)]
      evsi <- mean(winnertnbs) - maxenb
      oa <- mean(winnertnbs == maxnbs)

      if (voi) {
        return((evsi / evpi - target)^2)
      } else {
        return((oa - target)^2)
      }
    }

    if (b_voi) {
      target <- targets$voi.nb
      res <- OOR::StoSOO(
        c(1000),
        f,
        lower = 100,
        upper = 10^5,
        nb_iter = 1000,
        voi = TRUE
      )
      out$N <- unlist(c(out$N, voi.nb = round(res$par)))
    }
    if (b_oa) {
      target <- targets$oa.nb
      res <- OOR::StoSOO(
        c(1000),
        f,
        lower = 100,
        upper = 10^5,
        nb_iter = 1000,
        voi = FALSE
      )
      out$N <- unlist(c(out$N, oa.nb = round(res$par)))
    }
  }

  out$results <- out$N
  out$N <- NULL
  out$targets <- targets
  out$sample <- sample
  out
}
