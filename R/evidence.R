#' @export
summary.bpm_evidence <- function(object, ..., digits = getOption("digits")) {
  stopifnot(inherits(object, "bpm_evidence"))

  tab <- data.frame(
    component = names(object),
    type = character(length(object)),
    parameter1 = double(length(object)),
    parameter2 = double(length(object)),
    mean = double(length(object)),
    variance = double(length(object)),
    q2.5 = double(length(object)),
    q97.5 = double(length(object)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (i in seq_along(object)) {
    el <- object[[i]]
    tab$type[i] <- el$type
    tab$parameter1[i] <- el$parms[1]
    tab$parameter2[i] <- el$parms[2]
    tab$mean[i] <- el$moments[1]
    tab$variance[i] <- el$moments[2]
    tab[i, "q2.5"] <- dist_quantile(el$type, el$parms, 0.025)
    tab[i, "q97.5"] <- dist_quantile(el$type, el$parms, 0.975)
  }

  class(tab) <- c("summary.bpm_evidence", "data.frame")
  tab
}


#' @export
print.summary.bpm_evidence <- function(x, ..., digits = getOption("digits")) {
  cat("Summary of BPM evidence\n")
  cat("-----------------------\n\n")

  tab <- x
  num_cols <- c("parameter1", "parameter2", "mean", "variance", "q2.5", "q97.5")
  for (cn in intersect(num_cols, colnames(tab))) {
    tab[[cn]] <- round(tab[[cn]], digits)
  }

  print.data.frame(tab, row.names = FALSE)

  invisible(x)
}


#' @export
print.bpm_evidence <- function(x, ..., digits = 4) {
  cat("Bayesian prediction-model evidence\n")
  cat("(", length(x), " component", if (length(x) != 1) "s" else "", ")\n\n",
      sep = "")

  tab <- data.frame(
    component = names(x),
    distribution = vapply(x, function(e) e$type, character(1)),
    parameters = vapply(
      x,
      function(e) {
        paste(names(e$parms), signif(unname(e$parms), digits),
              sep = " = ", collapse = ", ")
      },
      character(1)
    ),
    mean = vapply(x, function(e) signif(e$moments[[1]], digits), numeric(1)),
    variance = vapply(x, function(e) signif(e$moments[[2]], digits), numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  print.data.frame(tab, row.names = FALSE, right = FALSE)
  invisible(x)
}


#'Transforms Evidence Into Standardized Format
#'
#'@description
#'Verifies that an evidence object has the required members and standardizes it
#'into a `bpm_evidence` object. Each element's distribution is recorded together
#'with both its native parameters (`$parms`) and its first two moments
#'(`$moments`, named `m` and `v`).
#'
#'@param evidence
#'A named list of evidence elements. The required members are:
#'\itemize{
#'  \item `prev`: outcome prevalence (defaults to a `beta` distribution),
#'  \item `cstat`: c-statistic (defaults to a `beta` distribution),
#'  \item `cal_slp`: calibration slope (defaults to `norm`), and
#'  \item exactly one of `cal_mean` (mean calibration), `cal_oe`
#'    (observed-to-expected ratio), or `cal_int` (calibration intercept),
#'    each defaulting to `norm`.
#'}
#'
#'Each element may be given either as a formula, `name ~ dist(par1, par2)`, or
#'as a named list, `name = list(type = "dist", ...)`. The supported
#'distributions (`type`) are `"norm"`, `"beta"`, `"logitnorm"`, and
#'`"probitnorm"`.
#'
#'@details
#'The two parameters of each element may be characterized flexibly, as either
#'native distribution parameters or summary moments. The parameters must be
#'**either all unnamed or all named** (a mix such as `beta(0.4, var = 0.04)` is
#'ambiguous and raises an error).
#'
#'**Unnamed (positional)** parameters are taken as the native parameters of the
#'distribution:
#'\itemize{
#'  \item `norm(mean, sd)`,
#'  \item `beta(shape1, shape2)`,
#'  \item `logitnorm(mu, sigma)`,
#'  \item `probitnorm(mu, sigma)`.
#'}
#'
#'**Named** parameters are matched against the following aliases (pick one pair
#'per element):
#'\itemize{
#'  \item moments with a variance: `mean`/`var` (or `m`/`v`),
#'  \item moments with a standard deviation: `mean`/`sd` (or `m`/`sd`),
#'  \item a mean and an upper 97.5\% quantile bound: `mean`/`cih` (or `m`/`cih`),
#'  \item native `beta` parameters: `alpha`/`beta`,
#'  \item native `logitnorm`/`probitnorm` parameters: `mu`/`sigma`.
#'}
#'When moments are supplied, the native parameters are obtained by the method of
#'moments (or, for `cih`, by matching the requested quantile).
#'
#'@return A `bpm_evidence` object: the standardized, restructured evidence list.
#'@examples
#'# Formula form, mixing native parameters and moments:
#'evidence <- list(
#'  prev     ~ beta(116, 155),       # native beta parameters
#'  cstat    ~ beta(mean = 0.76, sd = 0.006),
#'  cal_mean ~ norm(-0.009, 0.125),
#'  cal_slp  ~ norm(0.995, 0.024))
#'process_evidence(evidence = evidence)
#'
#'# Equivalent named-list form:
#'evidence <- list(
#'  prev=list(type="beta", mean=0.38, sd=0.2),
#'  cstat=list(mean=0.7, sd=0.05),
#'  cal_int=list(mean=0.2, sd=0.2),
#'  cal_slp=list(mean=0.8, sd=0.3))
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
        stop("Evidence element not of format x~dist(parm1, parm2):", this)
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
      "No valid parameter specification for calibration (calibration slope AND at least one of intercept, mean calibration, or O/E ratio are required)"
    )
  }

  cal_mean <- match(c("cal_mean"), names(cal_parms))
  cal_int <- match(c("cal_int"), names(cal_parms))
  cal_slp <- match(c("cal_slp"), names(cal_parms))
  cal_oe <- match(c("cal_oe"), names(cal_parms))

  if (!is.na(cal_mean)) {
    if (is.null(cal_parms[[cal_mean]]$type)) {
      cal_parms[[cal_mean]]$type <- "norm"
      message("Assuming normal distribution for calibration mean")
    }
    evidence$cal_mean <- process_evidence_element(cal_parms[[cal_mean]])
  }
  if (!is.na(cal_int)) {
    if (is.null(cal_parms[[cal_int]]$type)) {
      cal_parms[[cal_int]]$type <- "norm"
      message("Assuming normal distribution for calibration intercept")
    }
    evidence$cal_int <- process_evidence_element(cal_parms[[cal_int]])
  }
  if (!is.na(cal_slp)) {
    if (is.null(cal_parms[[cal_slp]]$type)) {
      cal_parms[[cal_slp]]$type <- "norm"
      message("Assuming normal distribution for calibration slope")
    }
    evidence$cal_slp <- process_evidence_element(cal_parms[[cal_slp]])
  }
  if (!is.na(cal_oe)) {
    if (is.null(cal_parms[[cal_oe]]$type)) {
      cal_parms[[cal_oe]]$type <- "norm"
      message("Assuming normal distribution for O/E ratio")
    }
    evidence$cal_oe <- process_evidence_element(cal_parms[[cal_oe]])
  }

  return(structure(evidence, class = "bpm_evidence"))
}


process_evidence_element <- function(element) {
  e <- list()
  e$type <- element$type
  element$type <- NULL

  # Decide how the two parameters were supplied. They must be *either* all
  # unnamed (positional native parameters) or all named (matched by alias).
  # A mix is ambiguous (e.g. beta(0.4, var = 0.04)) and is rejected.
  nms <- names(element)
  if (is.null(nms)) {
    nms <- rep("", length(element))
  }
  n_empty <- sum(nchar(nms) == 0)

  if (n_empty == length(element)) {
    # Two unnamed parameters: interpret as the native parameters of the type.
    if (e$type == "norm") {
      e$parms <- c(mean = element[[1]], sd = element[[2]])
      e$moments <- c(m = element[[1]], v = element[[2]]^2)
    }
    if (e$type == "beta") {
      e$parms <- c(shape1 = element[[1]], shape2 = element[[2]])
      e$moments <- moments("beta", e$parms)
    }
    if (e$type == "logitnorm") {
      e$parms <- c(mu = element[[1]], sigma = element[[2]])
      e$moments <- moments("logitnorm", e$parms)
    }
    if (e$type == "probitnorm") {
      e$parms <- c(mu = element[[1]], sigma = element[[2]])
      e$moments <- moments("probitnorm", e$parms)
    }
  } else if (n_empty == 0) {
    # Named parameters. Accepted aliases (per element, pick one pair):
    #   moments: (mean, var) / (m, v), (mean, sd) / (m, sd)
    #   mean + upper 97.5% quantile: (mean, cih) / (m, cih)
    #   native: (alpha, beta) for beta; (mu, sigma) for logit-/probit-normal
    possible_args <- list(
      c("mean", "var"),
      c("m", "v"),
      c("mean", "sd"),
      c("m", "sd"),
      c("alpha", "beta"),
      c("mu", "sigma"),
      c("mean", "cih"),
      c("m", "cih")
    )
    renamed_args <- list(
      c("m", "v"),
      c("m", "v"),
      c("m", "s"),
      c("m", "s"),
      c("shape1", "shape2"),
      c("mu", "sigma"),
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
      stop(
        "No valid parameter specification. Named parameters must be one of: ",
        "(mean, var), (mean, sd), (mean, cih), (alpha, beta), or (mu, sigma)."
      )
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
      # Native parameters supplied by name (e.g. beta(alpha, beta) or
      # logitnorm(mu, sigma)). `parms` is a list here, so coerce to a named
      # numeric vector before use; passing a list to moments() triggers
      # "non-numeric argument to binary operator".
      e$parms <- unlist(parms)
      e$moments <- moments(e$type, e$parms)
    }
  } else {
    stop(
      "Evidence parameters must be either all named or all unnamed; ",
      "a mix of named and unnamed values is ambiguous."
    )
  }

  # Normalise moment names across all specification styles. Different code
  # paths above otherwise return c("mean","var") or c("m.shape1", ...) etc.,
  # which is fragile for any name-based downstream access.
  if (!is.null(e$moments)) {
    names(e$moments) <- c("m", "v")
  }

  e
}
