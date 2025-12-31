#' @export
summary.bpm_evidence <- function(object, ..., digits = getOption("digits")) {
  stopifnot(inherits(object, "bpm_evidence"))

  # ---- extract core elements ----
  # (adjust names to your actual list structure)
  n_models <- length(object$models)
  weights <- object$weights

  res <- list(
    n_models = n_models,
    weights = weights,
    call = object$call
  )

  class(res) <- "summary.bpm_evidence"
  res
}


#' @export
print.summary.bpm_evidence <- function(x, ..., digits = getOption("digits")) {
  cat("Summary of BPM evidence\n")
  cat("-----------------------\n\n")

  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }

  cat("Number of models:", x$n_models, "\n\n")

  if (!is.null(x$weights)) {
    cat("Model weights:\n")
    print(round(x$weights, digits))
  }

  invisible(x)
}


#'Transforms Evidence Into Standardized Format
#'
#'@description Verifies evidence object has correct members, and standardizes it
#'@param evidence named list of evidence elements including:
#' prev: prevalence
#' cstat: c-statistic
#' cal_slp: calibration slope and,
#' one of cal_mean (mean calibration), cal_oe (observed to expected ratio), or cal_int (calibration intercept)
#'@return Modified evidence object that has been standardized and restructured
#'@examples
#'evidence <- list(
#' prev=list(type="beta", mean=0.38, sd=0.2),
#' cstat=list(mean=0.7, sd=0.05),
#' cal_int=list(mean=0.2, sd=0.2),
#' cal_slp=list(mean=0.8, sd=0.3))
#'process_evidence(evidence=evidence)
#'@export
process_evidence <- function(evidence) {
  #Transform formulas to list
  if (!is.list(evidence)) {
    stop("evidence should be specified as a list of formulas")
  }
  evidence2 <- list()
  n <- length(evidence)
  for (i in 1:n) {
    out <- list()
    this <- evidence[[i]]
    if (any(!is.na(match(class(this), "formula")))) {
      if (length(this) != 3) {
        stop("Unusual evidence specification:", this)
      }
      if (this[[1]] != "~") {
        stop("Evidence element not of formay x~dist(parm1, parm2):", this)
      }
      parm <- as.character(this[[2]])
      if (any(!is.na(match(names(evidence2), parm)))) {
        stop("Duplicate parameter specification for:", parm)
      }
      out[[parm]] <- list()
      out[[parm]]["type"] <- as.character(this[[3]][[1]])
      nm1 <- names(this[[3]][2])
      if (length(nm1) == 0) {
        out[[parm]][[2]] <- eval(this[[3]][[2]])
      } else {
        out[[parm]][nm1] <- eval(this[[3]][[2]])
      }
      nm2 <- names(this[[3]][3])
      if (length(nm2) == 0) {
        out[[parm]][[3]] <- eval(this[[3]][[3]])
      } else {
        out[[parm]][nm2] <- eval(this[[3]][[3]])
      }
      evidence2[parm] <- out
    } else {
      #evidence is not a formula, so must be a typed list
      parm <- names(evidence)[i]
      if (length(parm) == 0) {
        stop("Evidence element ", i, " is a list but does not have a 'type'")
      }
      if (any(!is.na(match(names(evidence2), parm)))) {
        stop("Duplicate parameter specification for:", parm)
      }
      evidence2[parm] <- list(this)
    }
  }

  evidence <- evidence2

  if (is.null(evidence$prev)) {
    stop("evidence object must have a prev (prevalence) member")
  }
  if (is.null(evidence$prev$type)) {
    evidence$prev$type <- "beta"
    message("Assuming prev has a beta distribution")
  }
  evidence$prev <- process_evidence_element(evidence$prev)

  if (is.null(evidence$cstat)) {
    stop("evidence object must have a cstat (c-statistic) member")
  }
  if (is.null(evidence$cstat$type)) {
    evidence$cstat$type <- "beta"
    message("Assuming cstat has a beta distribution")
  }
  evidence$cstat <- process_evidence_element(evidence$cstat)

  nms <- names(evidence)
  possible_args <- list(
    c("cal_mean", "cal_slp"),
    c("cal_oe", "cal_slp"),
    c("cal_int", "cal_slp")
  )
  renamed_args <- list(
    c("cal_mean", "cal_slp"),
    c("cal_oe", "cal_slp"),
    c("cal_int", "cal_slp")
  )
  cal_parms <- c()
  for (i in 1:length(possible_args)) {
    res <- match(possible_args[[i]], nms)
    if (!any(is.na(res))) {
      if (length(cal_parms) > 0) {
        stop("Multiple arguments matched")
      }
      cal_parms <- evidence[res]
      names(cal_parms) <- renamed_args[[i]]
    }
  }
  if (length(cal_parms) == 0) {
    stop(
      "No valid parameter specification for calibration (calibration slope AND at lease one of intercept, mean calibration, or O/E ratio are required)"
    )
  }

  cal_mean <- match(c("cal_mean"), names(cal_parms))
  cal_int <- match(c("cal_int"), names(cal_parms))
  cal_slp <- match(c("cal_slp"), names(cal_parms))
  cal_oe <- match(c("cal_oe"), names(cal_parms))

  if (!is.na(cal_mean)) {
    if (is.null(cal_parms[[cal_mean]]$type)) {
      cal_parms[[cal_mean]]$type <- "norm"
      message("Assuming normal distirbution for calibration mean")
    }
    evidence$cal_mean <- process_evidence_element(cal_parms[[cal_mean]])
  }
  if (!is.na(cal_int)) {
    if (is.null(cal_parms[[cal_int]]$type)) {
      cal_parms[[cal_int]]$type <- "norm"
      message("Assuming normal distirbution for calibration interccept")
    }
    evidence$cal_int <- process_evidence_element(cal_parms[[cal_int]])
  }
  if (!is.na(cal_slp)) {
    if (is.null(cal_parms[[cal_slp]]$type)) {
      cal_parms[[cal_slp]]$type <- "norm"
      message("Assuming normal distirbution for calibration slope")
    }
    evidence$cal_slp <- process_evidence_element(cal_parms[[cal_slp]])
  }
  if (!is.na(cal_oe)) {
    if (is.na(cal_parms[[cal_oe]]$type)) {
      cal_parms[[cal_oe]]$type <- "norm"
      message("Assuming normal distirbution for O/E ratio")
    }
    evidence$cal_oe <- process_evidence_element(cal_parms[[cal_oe]])
  }

  return(structure(evidence, class = "bpm_evidence"))
}


process_evidence_element <- function(element) {
  e <- list()
  e$type <- element$type
  element$type <- NULL
  if (any(nchar(names(element)) == 0)) {
    #We should expect two unnamed parameters
    if (e$type == "norm") {
      e$parms <- c(mean = element[[1]], sd = element[[2]])
      e$moments <- c(m = element[[1]], v = element[[2]]^2)
    }
    if (e$type == "beta") {
      e$parms <- c(shape1 = element[[1]], shape2 = element[[2]])
      e$moments <- moments("beta", e$parms)
    }
    if (e$type == "logitnorm") {
      e$parms <- c(mean = element[[1]], sd = element[[2]])
      e$moments <- moments("logitnorm", e$parms)
    }
    if (e$type == "probitnorm") {
      e$parms <- c(mean = element[[1]], sd = element[[2]])
      e$moments <- moments("probit", e$parms)
    }
  } else {
    ##Should have any of the following members: (mu, var), (mu, sd), (alpha, beta)
    nms <- names(element)
    possible_args <- list(
      c("mean", "var"),
      c("m", "v"),
      c("mean", "sd"),
      c("m", "sd"),
      c("alpha", "beta"),
      c("mean", "cih"),
      c("m", "cih")
    )
    renamed_args <- list(
      c("m", "v"),
      c("m", "v"),
      c("m", "s"),
      c("m", "s"),
      c("shape1", "shape2"),
      c("m", "cih"),
      c("m", "cih")
    )
    parms <- c()
    for (i in 1:length(possible_args)) {
      res <- match(possible_args[[i]], nms)
      if (!any(is.na(res))) {
        if (length(parms) > 0) {
          stop("Multiple arguments matched")
        }
        parms <- element[res]
        names(parms) <- renamed_args[[i]]
      }
    }
    if (length(parms) == 0) {
      stop("No valid parameter specification")
    }
    ##IF (m,s), apply the method of moments. Note that for normal it is a stupid tail chasing but OK!
    if (names(parms)[1] == 'm') {
      m <- parms[[1]]
      if (names(parms)[2] == "cih") {
        e$parms <- inv_mean_quantile(e$type, m, q = parms[[2]], p = 0.975)
        e$moments <- c(m = m, v = 0)
        e$moments[2] <- moments(e$type, e$parms)[[2]]
      } else {
        v <- parms[[2]]
        if (names(parms)[2] == "s") {
          v <- v^2
        }
        e$moments <- c(m = m, v = v)
        e$parms <- inv_moments(e$type, list(m, v))
      }
    } else {
      e$parms <- parms
      e$moments <- moments(e$type, parms)
    }
  }

  e
}
