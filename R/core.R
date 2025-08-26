logit <- function(x) {log(x/(1-x))}
expit <- function(x) {1/(1+exp(-x))}

#' Mean and Variance Calculator
#' 
#' @description Calculates the first two moments (mean and variance) of the given model type and parameters.
#' @param type The distribution type, one of c("norm", "beta", "logitnorm", "probitnorm").
#' @param parms A numeric vector containing parameters relevant to the model.
#' @return A numeric vector representing the mean and variance.
#' @examples
#' moments("norm", c(0, 1))
#' moments("beta", c(1, 1))
#'@export
moments <- function(type, parms) {
  # Check if the input values are valid
  if(type=="norm")
  {
    return(c(m=parms[1], v=(parms[2])^2))
  }
  if(type=="beta")
  {
    a <- parms[1]; b <- parms[2]
    return(c(m=a/(a+b), v=a*b/(((a+b)^2)*(a+b+1))))
  }
  if(type=="logitnorm")
  {
    return(momentsLogitnorm(parms[[1]], parms[[2]]))
  }
  if(type=="probitnorm")
  {
    mu <- parms[1]; sigma <- parms[2]
    m <- pnorm(mu/sqrt(1+sigma^2))
    ex2 <- integrate(function(x) {x^2*mcmapper::dprobitnorm(x,mu,sigma)}, 0,1)$value
    return(c(m=m, v=ex2-m^2))
  }
}


#'Calculates the Model Parameters Given Moments
#'
#'@description Calculates the model parameters of interest given the first two moments. 
#'@param type The distribution type, one of c("norm", "beta", "logitnorm").
#'@param moments A numeric vector containing the first two moments of the model
#'@return Returns the two parameters for each model.
#'  mean and sd for norm
#'  mu and sigma for logitnorm
#'  shape1 (alpha) and shape2 (beta) for beta
#'@examples 
#'inv_moments("norm", c(0,1))
#'@export
inv_moments <- function(type, moments) {
  # Check if the input values are valid
  m <- moments[[1]]
  v <- moments[[2]]
  
  if(type=="norm")
  {
    sd <- sqrt(v)
    return(c(mean=m, sd=sd))
  }
  
  if(type=="logitnorm")
  {
    if (m <= 0 || m >= 1) {
      stop("Mean must be between 0 and 1 (exclusive).")
    }
    if (v <= 0) {
      stop("Variance must be positive.")
    }
    
    objective <- function(params) {
      mu <- params[1]
      sigma <- params[2]
      
      tmp <- logitnorm::momentsLogitnorm(mu, sigma)
      mean_computed <- tmp[1]
      var_computed <- tmp[2]
      
      # Sum of squared differences between target and computed values
      (mean_computed - m)^2 + (var_computed - v)^2
    }
    
    # Initial guesses for mu and sigma
    initial_guess <- c(0, 1)
    
    # Use optim() to minimize the objective function
    result <- optim(initial_guess, objective, method = "L-BFGS-B", lower = c(-Inf, 1e-6), upper = c(Inf, Inf))
    
    # Extract mu and sigma
    mu <- result$par[1]
    sigma <- result$par[2]
    
    return(c(mu = mu, sigma = sigma))
  }
  
  if(type=="beta")
  {
    if (m <= 0 || m >= 1) {
      stop("Mean must be between 0 and 1 (exclusive).")
    }
    if (v <= 0) {
      stop("Variance must be positive.")
    }
    
    # Calculate alpha and beta
    alpha <- m * ((m * (1 - m) / v) - 1)
    beta <- (1 - m) * ((m * (1 - m) / v) - 1)
    
    # Return the alpha and beta as a named list
    return(c(shape1=alpha, shape2=beta))
  }
}

#'Calculates the Model Parameters Given Quantile
#' 
#'@description
#'Calculate the model parameters given the distribution type, mean, quantile, and percentile.
#'@param type The distribution type, one of c("norm", "beta", "logitnorm", "probitnorm").
#'@param m Mean of the of distribution.
#'@param q The quantile value.
#'@param p The percentile at which the quantile occurs.
#'@return The model parameters of the given type.
#'@examples inv_mean_quantile("beta", 0.5, 0.25, 0.25)
#'@export
inv_mean_quantile <- function(type, m, q, p)
{
  if(type=="logitnorm")
  {
    res <- tryCatch(
      {logitnorm::twCoefLogitnormE(mean=m, quant=q, perc=p)},
      error=(function(cond) {logitnorm::twCoefLogitnormE(mean=m, quant=q, perc=p,theta0=logitnorm::twCoefLogitnorm(m,q,p))}))
    out <- c(mu=res[1], sigma=res[2])
  }
  if(type=="beta") 
  {
    res <- uniroot(function(x) {pbeta(q,x,x*(1-m)/m)-p}, interval=c(0.0001,10000))
    out <- c(shape1=res$root, shape2=res$root*(1-m)/m)
  }
  if(type=="probitnorm") #TODO
  {
    # res <- uniroot(function(x) {mcmapper:::pprobitnorm(q,x,x*(1-m)/m)-p}, interval=c(0.0001,10000))
    # out <- c(shape1=res$root, shape2=res$root*(1-m)/m)
  }
  if(type=="norm")
  {
    res <- uniroot(function(x) {pnorm(q,m,x)-p}, interval=c(0.0001,10))
    out <- c(mean=m, sd=res$root)
  }
  
  out
}

#' Calculates the C-statistic of Model
#' 
#' @description Calculates the c-statistic given the model type and parameters.
#' @param type A character string; one of c("beta", "logitnorm", "probitnorm") indicating the model type.
#' @param parms A numeric vector containing parameters relevant to the model.
#' @param m Mean, default is NULL
#' @return The C-statistic
#' @examples
#' calc_cstat("norm", c(0,2))
#' @export
calc_cstat <- function(type, parms, m=NULL) #For now we assume we know m
{
  if(type=="logitnorm")
  {
    if(is.null(m)) {m <- mcmapper:::elogitnorm(parms[1], parms[2])}
    C <- (1-integrate(f=function(x){mcmapper::plogitnorm(x,parms[1],parms[2])^2}, lower=0, upper=1)$value-m^2)/(2*m*(1-m))
  }
  if(type=="beta")
  {
    if(is.null(m)) {m <- parms[1]/(parms[1]+parms[2])}
    C <- (1-integrate(f=function(x){pbeta(x,parms[1],parms[2])^2}, lower=0, upper=1)$value-m^2)/(2*m*(1-m))
  }
  if(type=="probitnorm")
  {
    if(is.null(m)) {m <- pnorm(parms[1]/sqrt(1+parms[2]^2))} #TODO: check
    C <- (1-integrate(f=function(x){mcmapper::pprobitnorm(x,parms[1],parms[2])^2}, lower=0, upper=1)$value-m^2)/(2*m*(1-m))
  }
  
  C
}


#'Infer Calibration Intercept from Mean Calibration
#'
#'@description Infer calibration intercept from mean calibration given a fixed calibration slope and a given distribution for calibrated risks
#'@param dist_type The distribution type, one of c("logitnorm", "probitnorm", "beta").
#'@param dist_parms The two parameters that index the type.
#'@param cal_mean The mean calibration.
#'@param cal_slp The calibration slope.
#'@param prev Outcome prevalence. Optional; if not provided, estimate is as the expected value of the distribution of calibrated risks.
#'@return The estimated calibration intercept
#'@examples
#'infer_cal_int_from_mean("beta", c(1,1), 1 ,1.1, 0.25)
#'@export
infer_cal_int_from_mean <- function(dist_type, dist_parms, cal_mean, cal_slp, prev=NULL) #TODO: prev
{
  if(dist_type=="logitnorm")
  {
    f <- function(intercept)
    {
      integrate(function(x, intercept, slope){expit((logit(x)-intercept)/cal_slp)*mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])}, 0, 1,intercept=intercept, slope=cal_slp)$value
    }
    res <- optim(0, fn=function(x) {((prev-f(x))-cal_mean)^2}, lower=-2, upper=2, method="Brent")
  }

  if(dist_type=="beta")
  {
    f <- function(intercept)
    {
      integrate(function(x, intercept, slope){expit((logit(x)-intercept)/cal_slp)*dbeta(x,dist_parms[1],dist_parms[2])}, 0, 1,intercept=intercept, slope=cal_slp)$value
    }
    res <- optim(0, fn=function(x) {((prev-f(x))-cal_mean)^2}, lower=-2, upper=2, method="Brent")
  }

  if(dist_type=="probitnorm")
  {
    f <- function(intercept)
    {
      integrate(function(x, intercept, slope){expit((logit(x)-intercept)/cal_slp)*mcmapper::dprobitnorm(x,dist_parms[1],dist_parms[2])}, 0, 1,intercept=intercept, slope=cal_slp)$value
    }
    res <- optim(0, fn=function(x) {((prev-f(x))-cal_mean)^2}, lower=-4, upper=4, method="Brent")
  }

  unname(res$par)
}


#'Infer Calibration Intercept from O/E ratio
#' 
#'@description Infer calibration intercept from observed-to-expected outcome ratio given a fixed calibration slope and a given distribution for calibrated risks
#'@param dist_type The distribution type, one of c("logitnorm", "probitnorm", "beta").
#'@param dist_parms The two parameters that index the type.
#'@param cal_oe The observed-to-expected outcome ratio.
#'@param cal_slp The calibration slope.
#'@param prev Outcome prevalence. Optional; if not provided, estimate is as the expected value of the distribution of calibrated risks.
#'@return The estimated calibration intercept
#'@examples
#'infer_cal_int_from_oe("beta", c(1,1), 0.9, 1.1, 0.25)
#'@export
infer_cal_int_from_oe <- function(dist_type, dist_parms, cal_oe, cal_slp, prev=NULL) #TODO: prev
{
  if(dist_type=="logitnorm")
  {
    f <- function(intercept)
    {
      integrate(function(x, intercept, slope){expit((logit(x)-intercept)/cal_slp)*mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])}, 0, 1,intercept=intercept, slope=cal_slp)$value
    }
    res <- optim(0, fn=function(x) {((prev/f(x))-cal_oe)^2}, lower=-2, upper=2, method="Brent")
  }
  
  if(dist_type=="beta")
  {
    f <- function(intercept)
    {
      integrate(function(x, intercept, slope){expit((logit(x)-intercept)/cal_slp)*dbeta(x,dist_parms[1],dist_parms[2])}, 0, 1,intercept=intercept, slope=cal_slp)$value
    }
    res <- optim(0, fn=function(x) {((prev/f(x))-cal_oe)^2}, lower=-2, upper=2, method="Brent")
  }
  
  if(dist_type=="probitnorm")
  {
    f <- function(intercept)
    {
      integrate(function(x, intercept, slope){expit((logit(x)-intercept)/cal_slp)*mcmapper::dprobitnorm(x,dist_parms[1],dist_parms[2])}, 0, 1,intercept=intercept, slope=cal_slp)$value
    }
    res <- optim(0, fn=function(x) {((prev/f(x))-cal_oe)^2}, lower=-4, upper=4, method="Brent")
  }
  
  unname(res$par)
}

#'Calculates the Sensitivity and Specificity
#'
#'@description Calculate the sensitivity and specificity of the model at given threshold
#'@param dist_type The distribution type, one of c("logitnorm", "beta", "probitnorm").
#'@param dist_parms Vector of the two parameters of interest given the distribution.
#'@param cal_int The calibration intercept.
#'@param cal_slp The calibration slope.
#'@param threshold The risk threshold
#'@param prev The outcome prevalence, the expectation of the model
#'@return A vector containing sensitivity and specificity
#'@examples
#'calc_se_sp("beta", c(1,1), 0.9, 0.75, 0.5, 0.5)
#'@export
calc_se_sp <- function(dist_type, dist_parms, cal_int, cal_slp, threshold, prev) #TODO: prev should be optional
{

  hz <- expit(logit(threshold)*cal_slp+cal_int)

  if(dist_type=="logitnorm")
  {
    tp <- integrate(f=function(x) {x*mcmapper::dlogitnorm(x, dist_parms[1], dist_parms[2])}, hz, 1)$value
    se <- unname(tp/prev)
    sp <- unname((mcmapper::plogitnorm(hz, dist_parms[1], dist_parms[2])-(prev-tp))/(1-prev))
  }

  if(dist_type=="beta")
  {
    tp <- integrate(f=function(x) {x*dbeta(x, dist_parms[1], dist_parms[2])}, hz, 1)$value
    se <- unname(tp/prev)
    sp <- unname((pbeta(hz, dist_parms[1], dist_parms[2])-(prev-tp))/(1-prev))
  }

  if(dist_type=="probitnorm")
  {
    tp <- integrate(f=function(x) {x*mcmapper::dprobitnorm(x, dist_parms[1], dist_parms[2])}, hz, 1)$value
    se <- unname(tp/prev)
    sp <- unname((mcmapper::pprobitnorm(hz, dist_parms[1], dist_parms[2])-(prev-tp))/(1-prev))
  }

  return(c(se=max(0,min(se,1)), sp=max(0,min(sp,1))))
}

#' Calculates Correlation
#'
#'@description Calculates correlation based on simulated data
#'@param dist_type The distribution type
#'@param dist_parms The two parameters of interest for the given distribution type
#'@param cal_int The calibration intercept.
#'@param cal_slp The calibration slope.
#'@param n number of observations for each simulation.
#'@param n_sim number of simulations
#'@return correlation among the simulated data
#'@examples
#'infer_correlation("beta", c(1,1), 1, 0.5, 10, 10)
#'@export
infer_correlation <- function(dist_type, dist_parms, cal_int, cal_slp, n, n_sim)
{
  # require(pROC)
  # require(mcmapper)

  out <- matrix(NA, nrow=n_sim, ncol=5)
  colnames(out) <- c("prev", "cstat", "cal_mean", "cal_int", "cal_slp")

  for(i in 1:n_sim)
  {
    p <- do.call(paste0("r",dist_type), args=as.list(c(n,unname(dist_parms))))
    Y <- rbinom(n,1,p)
    pi <- expit((logit(p)-cal_int)/cal_slp)
    df <- cbind(pi=pi,Y=Y)
    out[i,]<-c(mean(df[2]),
               pROC::roc(df[,2]~df[,1], quiet=TRUE)$auc,
               mean(df[2]-df[1]),
               coef(glm(df[,2]~logit(df[,1]), family=binomial(link="logit")))[1:2])
  }

  cor(out)
}


#'
# express_evidence <- function(evidence, round_digits=3)
# {
#   out <- "";
#   for(i in 1:length(evidence))
#   {
#    nm <- names(evidence)[i]
#     item <- evidence[[i]]
#     out <- paste0(out, nm, "~", item$type, "(", round(item$parms[1], round_digits), ",", round(item$parms[2], round_digits),")\n")
#   }
#   
#   paste(out, collapse="\n")
# }


#'Calculates Approximate Variances and Covariance for Performance Metrics
#'
#'@description Calculates approximate variances performance metrics and covariance of calibration intercept and slope using the Riley framework
#'@param N sample size of the validation dataset
#'@param parms list containing model and distribution parameters:
#'  prev: expected prevalence
#'  cstat: c-statistic of the model
#'  dist_type: one of ("logitnorm", "beta, "probitnorm")
#'  dist_parm1: first parameter of the distribution
#'  dist_parm2: second parameter of the distribution
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'@return list of approximate variances and covariance of the performance metrics.
#'@examples 
#'parameters <- list(
#'   prev=0.5, 
#'   cstat=0.75, 
#'   dist_type="beta", 
#'   dist_parm1=1,
#'   dist_parm2=1, 
#'   cal_int=0, 
#'   cal_slp=1)
#'    
#'calc_riley_vars(N=100, parms=parameters)
#'@export
calc_riley_vars <- function(N, parms)
{
  prev <- parms$prev
  C <- parms$cstat
  dist_type <- parms$dist_type
  dist_parms <- c(parms$dist_parm1,parms$dist_parm2)
  cal_int <- parms$cal_int
  cal_slp <- parms$cal_slp
  bad_place <- NULL

  if(dist_type=="logitnorm")
  {
    #tryCatch({
      E_pi <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(expit((logit(x)-cal_int)/cal_slp))},0,1)$value
      I_a <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(x*(1-x))}, 0, 1)$value
      I_b <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)^2*(x*(1-x)))}, 0, 1)$value
      I_ab <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)*(x*(1-x)))}, 0, 1)$value
    #}, error=function(e){bad_place$parms <<- parms})
  }
  if(dist_type=="beta")
  {
    tryCatch({
    E_pi <- integrate(function(x){dbeta(x,dist_parms[1],dist_parms[2])*(expit((logit(x)-cal_int)/cal_slp))},0,1)$value
    I_a <- integrate(function(x){ dbeta(x,dist_parms[1],dist_parms[2])*(x*(1-x))}, 0, 1)$value
    I_b <- integrate(function(x){ dbeta(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)^2*(x*(1-x)))}, 0, 1)$value
    I_ab <- integrate(function(x){dbeta(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)*(x*(1-x)))}, 0, 1)$value
    }, error=function(e){bad_place$parms <<- parms})
    
  }
  if(dist_type=="probitnorm")
  {
    E_pi <- integrate(function(x){mcmapper::dprobitnorm(x,dist_parms[1],dist_parms[2])*(expit((logit(x)-cal_int)/cal_slp))},0,1)$value
    I_a <- integrate(function(x){mcmapper::dprobitnorm(x,dist_parms[1],dist_parms[2])*(x*(1-x))}, 0, 1)$value
    I_b <- integrate(function(x){mcmapper::dprobitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)^2*(x*(1-x)))}, 0, 1)$value
    I_ab <- integrate(function(x){mcmapper::dprobitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)*(x*(1-x)))}, 0, 1)$value
  }

  v_prev <- prev*(1-prev)/N
  v_cstat <- (C*(1-C)*(1+(N/2-1)*((1-C)/(2-C))+(N/2-1)*C/(1+C)))/(N*N*prev*(1-prev))
  v_cal_mean <- v_prev #TODO
  v_cal_oe <- (1-prev)/(N*E_pi)
  det <- N*(I_a*I_b-I_ab^2)
  v_cal_int <- I_b/det
  v_cal_slp <- I_a/det
  cov_int_slp <- -I_ab/det

  list(prev=v_prev,cstat=v_cstat, cal_mean=v_cal_mean, cal_oe=v_cal_oe,  cal_int=v_cal_int, cal_slp=v_cal_slp, cov_int_slp=cov_int_slp)
}

#'Calculates Sample Size that Achieves Target CI Widths
#'
#'@description Calculates sample size that achieves target confidence interval widths using Riley's framework
#'@param target_ciws Named list containing target confidence interval width for at least one of:
#'  prev: prevalence
#'  cstat: c-statistic
#'  cal_mean: mean calibration
#'  cal_oe: observed to expected outcome ratio
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'@param parms List containing model parameters and distribution:
#'  prev: expected prevalence
#'  cstat: c-statistic of the model
#'  dist_type: one of ("logitnorm", "beta, "probitnorm")
#'  dist_parm1: first parameter of the distribution
#'  dist_parm2: second parameter of the distribution
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'@return A named list of estimated sample sizes that achieve target confidence interval widths:
#'  fciw.prev, fciw.cstat, fciw.cal_mean, fciw.cal_oe, fciw.cal_int, fciw.cal_slp
#'@examples 
#'parameters <- list(
#'   prev=0.5, 
#'   cstat=0.75, 
#'   dist_type="beta", 
#'   dist_parm1=1, 
#'   dist_parm2=1, 
#'   cal_int=0, 
#'   cal_slp=1)
#'   
#'riley_samp(target_ciws=list(prev=0.05, cstat=0.05) ,parms=parameters)
#'@export
riley_samp <- function(target_ciws, parms)
{
  out <- list()
  K <- 2*qnorm(0.975)
  
  prev <- parms$prev
  C <- parms$cstat
  dist_type <- parms$dist_type
  dist_parms <- c(parms$dist_parm1,parms$dist_parm2)
  cal_int <- parms$cal_int
  cal_slp <- parms$cal_slp
  E_pi <- expit((logit(prev)-cal_int)/cal_slp)
  
  if(!is.null(target_ciws[['prev']]))
  {
    v <- (target_ciws[['prev']]/K)^2
    out$fciw.prev <- round(prev*(1-prev)/v)
  }
  
  if(!is.null(target_ciws[['cstat']]))
  {
    v <- (target_ciws[['cstat']]/K)^2
    A <- C * (1 - C)
    D <- ((1 - C) / (2 - C)) + (C / (1 + C))
    a <- v * prev * (1 - prev)
    b <- -A * D / 2
    c <- A * D - A
    discriminant <- b^2 - 4 * a * c
    N1 <- (-b + sqrt(discriminant)) / (2 * a)
    N2 <- (-b - sqrt(discriminant)) / (2 * a)
    out$fciw.cstat <- round(max(N1,N2))
  }
  
  if(!is.null(target_ciws[['cal_mean']]))
  {
    v <- (target_ciws[['cal_mean']]/K)^2 #TODO
    out$fciw.cal_mean <- round(prev*(1-prev)/v)
  }
  
  if(!is.null(target_ciws[['cal_oe']]))
  {
    v <- (target_ciws[['cal_oe']]/K)^2
    out$fciw.cal_oe <- round((1-prev)/(v*E_pi))
  }
  
  
  if(!is.null(target_ciws[['cal_int']]) | !is.null(target_ciws[['cal_slp']]))
  {
    v1 <- (target_ciws[['cal_int']]/K)^2
    v2 <- (target_ciws[['cal_slp']]/K)^2
    
    if(dist_type=="logitnorm")
    {
      #tryCatch({
      #E_pi <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(expit((logit(x)-cal_int)/cal_slp))},0,1)$value
      I_a <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(x*(1-x))}, 0, 1)$value
      I_b <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)^2*(x*(1-x)))}, 0, 1)$value
      I_ab <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)*(x*(1-x)))}, 0, 1)$value
      #}, error=function(e){bad_place$parms <<- parms})
    }
    if(dist_type=="beta")
    {
      #E_pi <- integrate(function(x){x*dbeta(expit(logit(x)*cal_slp+cal_int),dist_parms[1],dist_parms[2])}, 0, 1,intercept=cal_int, slope=cal_slp)$value
      # I_a <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*(exp(cal_int+x*cal_slp)/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
      # I_b <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*((x^2*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
      # I_ab <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*((x*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
      I_a <- integrate(function(x){dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*(exp(cal_int+x*cal_slp)/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
      I_b <- integrate(function(x){dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*((x^2*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
      I_ab <- integrate(function(x){dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*((x*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
    }
    if(dist_type=="probitnorm")
    {
      #E_pi <- integrate(function(x){x*mcmapper::dprobitnorm(expit(logit(x)*cal_slp+cal_int),dist_parms[1],dist_parms[2])}, 0, 1,intercept=cal_int, slope=cal_slp)$value
      I_a <- integrate(function(x){mcmapper::dprobitnorm(expit(logit(x)), dist_parms[1], dist_parms[2])*(exp(cal_int+x*cal_slp)/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
      I_b <- integrate(function(x){mcmapper::dprobitnorm(expit(logit(x)), dist_parms[1], dist_parms[2])*((x^2*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
      I_ab <- integrate(function(x){mcmapper::dprobitnorm(expit(logit(x)), dist_parms[1], dist_parms[2])*((x*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
    }
    
    if(length(v1)>0)  out$fciw.cal_int <- round(I_b/(v1*(I_a*I_b-I_ab^2)))
    if(length(v2)>0)  out$fciw.cal_slp <- round(I_a/(v2*(I_a*I_b-I_ab^2)))
  }
  out
}


#'Calculates Pre-Posterior Distribution of 95% CI Widths Using Two-step Method
#'
#'@description Calculates pre-posterior distribution of 95% CI widths using two-step method. 
#'@param N A vector of sample sizes
#'@param parms Parameters for the distribution containing:
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'  prev: prevalence
#'  dist_type: distribution type
#'  cstat: c-statistic 
#'  dist_type: one of ("logitnorm", "beta, "probitnorm")
#'  dist_parm1: first parameter of the distribution
#'  dist_parm2: second parameter of the distribution
#'@return List of length N, of vectors containing 95% confidence interval width for each of:
#'  cstat: c-statistic
#'  cal_oe: observed to expected ratio
#'  cal_mean: mean calibration
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'@examples
#'parameters <- list(
#' prev = 0.2, 
#' cstat = 0.75, 
#' dist_type="beta", 
#' dist_parm1=1, 
#' dist_parm2=2, 
#' cal_int = 0, 
#' cal_slp = 1)
#' 
#'calc_ciw_2s(N=c(100,300,500), parms=parameters)
#'@export
calc_ciw_2s <- function(N, parms)
{
  ciw_cstat <- ciw_cal_oe <- ciw_cal_mean <- ciw_cal_int <- ciw_cal_slp <- rep(NA, length(N))

  l <- length(N)

  s1 <- calc_riley_vars(N, parms=parms)
  
  hat_cals <- rbnorm(l,
                rep(parms$cal_int,l),
                rep(parms$cal_slp,l),
                s1$cal_int,
                s1$cal_slp,
                s1$cov_int_slp)

  k <- 2*qnorm(0.975)
  
  for(i in 1:l)
  {
    hat_prev <- inv_moments("beta",c(parms$prev, s1$prev[i]))
    hat_cstat <- inv_moments("beta",c(parms$cstat, s1$cstat[i]))

    f <- function()
    {
      repeat{
        x<-rbeta(1, hat_cstat[1], hat_cstat[2])
        if(x>=0.51 & x<0.99) break;
        }
      x
    }

    new_parms <- list(prev=rbeta(1, hat_prev[1], hat_prev[2]),
                    cstat=f(),
                    cal_int=hat_cals[i,1],
                    cal_slp=hat_cals[i,2])


    dist_parms <- mcmapper::mcmap_logitnorm(c(new_parms$prev, new_parms$cstat)) #TODO


    s2_vars <- calc_riley_vars(N[i], parms=list(prev=new_parms$prev,
                                        cstat=new_parms$cstat,
                                        cal_int=new_parms$cal_int,
                                        cal_slp=new_parms$cal_slp,
                                        dist_type="logitnorm",
                                        dist_parm1=dist_parms[1],
                                        dist_parm2=dist_parms[2]
                                        )
                        )
    ciw_cstat[i] <- k*sqrt(s2_vars$cstat)
    ciw_cal_oe[i] <- k*sqrt(s2_vars$cal_oe)
    ciw_cal_mean[i] <- k*sqrt(s2_vars$cal_mean)
    ciw_cal_int[i] <- k*sqrt(s2_vars$cal_int)
    ciw_cal_slp[i] <- k*sqrt(s2_vars$cal_slp)
  }

  if(is.infinite(sum(ciw_cal_oe))) {browser()}

  list(cstat=ciw_cstat,
       cal_oe=ciw_cal_oe,
       cal_mean=ciw_cal_mean,
       cal_int=ciw_cal_int,
       cal_slp=ciw_cal_slp)
}

#'Calculates Pre-Posterior Distribution of 95% CI Widths Using Sampling-based Simulation 
#'
#'@description Calculates pre-posterior distribution of 95% CI widths using sampling-based simulation 
#'@param N A vector of sample sizes
#'@param parms Parameters for the distribution containing:
#'  prev: prevalence
#'  dist_type: distribution type
#'  dist_parm1: first parameter of distribution
#'  dist_parm2: second parameter of distribution
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'@return List of length N, of vectors containing 95% confidence interval width for each of:
#'  cstat: c-statistic
#'  cal_oe: observed to expected ratio
#'  cal_mean: mean calibration
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'@examples
#'parameters<-list(prev = 0.2, dist_type="beta", dist_parm1=2, dist_parm2=4, cal_int = 0, cal_slp = 1)
#'calc_ciw_sample(N=c(100,300,500), parms=parameters)
#'@export
calc_ciw_sample <- function(N, parms)
{
  # require(fastLogisticRegressionWrap)
  se_cstat <- se_cal_oe <- se_cal_mean <- se_cal_int <- se_cal_slp <- rep(NA, length(N))

  prev <- parms$prev
  #C <- parms$cstat
  dist_type <- parms$dist_type
  dist_parms <- c(parms$dist_parm1, parms$dist_parm2)
  cal_int <- parms$cal_int
  cal_slp <- parms$cal_slp

  N_max <- max(N)
  p <- do.call(paste0("r",dist_type), as.list(c(N_max, dist_parms)))
  pi <- expit((logit(p)-cal_int)/cal_slp)
  #logit_pi <- logit(pi)
  logit_pi <- cbind(1,logit(pi))
  Y <- rbinom(N_max, 1, p)

  for(i in 1:length(N))
  {
    n <- N[i]
    O <- sum(Y[1:n])
    se_cal_oe[i] <- sqrt((1-O/n)/O)
    se_cal_mean[i] <- t.test(Y[1:n]-pi[1:n])$stderr
    reg <- fastLogisticRegressionWrap::fast_logistic_regression(logit_pi[1:n,],Y[1:n],do_inference_on_var="all")
    se_cal_int[i] <- reg$se[1]
    se_cal_slp[i] <- reg$se[2]

    #tmp <- ci.auc(Y[1:n]~logit_pi[1:n], quiet=TRUE)
    tmp <- ci.auc(Y[1:n]~logit_pi[1:n,2], quiet=TRUE)
    se_cstat[i] <- (tmp[3]-tmp[2])/qnorm(0.975)
  }

  k <- 2*qnorm(0.975)
  list(cstat=k*se_cstat, cal_oe=k*se_cal_oe, cal_mean=k*se_cal_mean, cal_int=k*se_cal_int, cal_slp=k*se_cal_slp)
}





#'#'Calculates Pre-Posterior Distribution of 95% CI Widths Based on Given Method
#'
#'@description Calculates pre-posterior distribution of 95% CI widths based on given method
#'@param N A vector of sample sizes
#'@param parms_sample Matrix of parameters for the distribution each row with appropriate parameters:
#'  cstat: c-statistic
#'  prev: prevalence
#'  dist_type: distribution type
#'  dist_parm1: first parameter of distribution
#'  dist_parm2: second parameter of distribution
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'@param method Method to calculate 95% confident interval width, one of sample, 2s
#'@return List of matrices each with dimension (number of rows in parms_sample x length N) containing 95% confidence interval width for each of:
#'  cstat: c-statistic
#'  cal_oe: observed to expected ratio
#'  cal_mean: mean calibration
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'@examples
#'parameters <- data.frame(cstat=c(0.75, 0.8), 
#'                         prev=c(0.2,0.3), 
#'                         dist_type=c("beta", "beta"), 
#'                         dist_parm1=c(0.9, 1),
#'                         dist_parm2=c(1.5,2),
#'                         cal_int=c(0, 0.5), 
#'                         cal_slp=c(1,0.9))
#'calc_ciw_mc(N=c(100,300,500), parms_sample=parameters, method="2s")
#'@export
calc_ciw_mc <- function(N, parms_sample, method)
{
  n <- nrow(parms_sample)

  ciw_cstat <- ciw_cal_oe <- ciw_cal_mean <- ciw_cal_int <- ciw_cal_slp <- matrix(NA, nrow=n, ncol=length(N))

  if(method=="sample")
  {
    for(i in 1:n)
    {
      tmp <- calc_ciw_sample(N, parms=parms_sample[i,])
      ciw_cstat[i, ] <- tmp$cstat
      ciw_cal_oe[i, ] <- tmp$cal_oe
      ciw_cal_mean[i, ] <- tmp$cal_mean
      ciw_cal_int[i, ] <- tmp$cal_int
      ciw_cal_slp[i, ] <- tmp$cal_slp
    }
  } else if(method=="2s")
  {
    for(i in 1:n)
    {
      tmp <- calc_ciw_2s(N, parms=parms_sample[i,])
      ciw_cstat[i, ] <- tmp$cstat
      ciw_cal_oe[i, ] <- tmp$cal_oe
      ciw_cal_mean[i, ] <- tmp$cal_mean
      ciw_cal_int[i, ] <- tmp$cal_int
      ciw_cal_slp[i, ] <- tmp$cal_slp
    }
  }

  list(cstat=ciw_cstat, cal_oe=ciw_cal_oe, cal_mean=ciw_cal_mean, cal_int=ciw_cal_int, cal_slp=ciw_cal_slp)
}

#'Calculates Sample Size Given Target Mean CI
#'
#'@description Calculates sample size N, so that the mean confidence interval is equal to given target, assumes function is decreasing and convex
#'@param target The target mean confidence interval width
#'@param N Sample sizes corresponding to each row of ciws,=
#'@param ciws Matrix of confidence intervals widths, each row corresponding to N
#'@param decreasing Logical. Constraining function to decreasing
#'@param convex Logical. Constraining function to convex
#'@return Integer. Estimated sample size needed to achieve the target
#'@examples
#'ciw_matrix <- matrix(c(0.09, 0.07, 0.06, 0.05), nrow=1)
#'find_n_mean(target=0.055, N=c(100, 200, 300, 400), ciws=ciw_matrix)
#'@export
find_n_mean <- function(target, N, ciws, decreasing=T, convex=T)
{
  # require(cobs)
  
  X <- rep(N, each=nrow(ciws))
  Y <- as.vector(ciws)
  
  constraint <- c(ifelse(decreasing,"decreasing",""),ifelse(convex,"convex",""))
  
  S <- cobs::cobs(X,Y, constraint=c("decrease","convex"))
  #plot(predict(S, seq(min(N),max(N),length.out=100)), type='l')
  
  f <- function(x) {predict(S, z=x)[,2] - target}
  round(uniroot(f, c(min(N),max(N)))$root,0)
}


#'Calculates Sample Size Given Target Quantile
#'
#'@description Find sample size N, so that the specified quantile is equal to given target
#'@param target The desired quantile target value
#'@param N Sample sizes corresponding to each row of ciws
#'@param q Desired quantile level, between 0 and 1.
#'@param ciws A matrix of confidence intervals widths, each row corresponding to N
#'@return Estimated sample size needed to achieve the target
#'@importFrom quantreg rqss qss
#'@examples
#'ciw_matrix <- matrix(c(
#'0.09, 0.08, 0.07, 0.06,
#'0.085, 0.075, 0.065, 0.055,
#'0.095, 0.085, 0.075, 0.065
#'), nrow = 3, byrow = TRUE)
#'
#'N_vals <- c(100, 200, 300, 400)
#'
#'find_n_quantile(target = 0.065, N = N_vals, q = 0.5, ciws = ciw_matrix)
#'@export
find_n_quantile <- function(target, N, q, ciws)
{
  # require(quantreg)
  # library(quantreg)
  
  X <- rep(N, each=nrow(ciws))
  Y <- as.vector(ciws)
  data <- data.frame(X=X, Y=Y)
  
  fit <- quantreg::rqss(Y ~ qss(X, lambda = 1), tau=q, data=data)
  
  # plot(X, Y, main = "Non-parametric Quantile Regression", xlab = "X", ylab = "Y")
  # for (i in 1:length(tau)) {
  #   lines(sort(X), predict(fit, newdata = data.frame(X = sort(X)), tau = tau[i]), col = i + 1)
  # }
  # legend("topleft", legend = paste("tau =", tau), col = 2:(length(tau) + 1), lty = 1)
  
  f <- function(x) {predict(fit, newdata=data.frame(X=x)) - target}
  round(uniroot(f, c(min(N),max(N)))$root,0)
}


#'Generates Samples From Normal Distribution
#'
#'@description generates samples from a normal distribution using marginal means, variances, and covariance
#'@param n Number of samples to be generated
#'@param mu1 Mean of first variable
#'@param mu2 Mean of second variable
#'@param var1 Variance of first variable
#'@param var2 Variance of second variable
#'@param cov Covariance between the two variables
#'@return Matrix of nx2 where
#'  column 1 contains samples for the first variable, and
#'  column 2 contains samples for the second variable conditioned on the first
#'@examples
#'rbnorm(n = 1000, mu1 = 0, mu2 = 0, var1 = 1, var2 = 1, cov = 0.5)
#'@export
rbnorm <- function(n, mu1, mu2, var1, var2, cov) {

  # Calculate standard deviations and correlation from variances and covariance
  sigma1 <- sqrt(var1)
  sigma2 <- sqrt(var2)
  rho <- cov / (sigma1 * sigma2)

  x1 <- rnorm(n, mu1, sigma1)

  conditional_mean <- mu2 + rho * (sigma2 / sigma1) * (x1 - mu1)
  conditional_sd <- sigma2 * sqrt(1 - rho^2)

  x2 <- rnorm(n, conditional_mean, conditional_sd)

  return(cbind(x1, x2))
}





#'Creates a Standardized List Including type, parms, moments
#'@description Standardizes element names, and creates a clean list of names
#'@param element Named list with type (e.g "beta", "norm", logitnorm"), and one of valid combination of parameters:
#'  - "mean" and "var"
#'  - "m" and "v"
#'  - "mean" and "sd"
#'  - "m" and "sd"
#'  - "alpha" and "beta"
#'  - "mean" and "cih"
#'  - "m" and "cih"
#'@return an element with type, params and moments.
#'@examples 
#'process_evidence_element(list(type="norm", mean=0, sd=1))
#'@export
process_evidence_element <- function(element)
{
  e <- list()
  e$type <- element$type
  
  if(any(nchar(names(element))==0)) {stop("Unnamed objects found in ...")}
  ##Should have any of the following members: (mu, var), (mu, sd), (alpha, beta)
  nms <- names(element)
  possible_args <- list(c("mean","var"), c("m","v"), c("mean","sd"), c("m","sd"), c("alpha","beta"), c("mean","cih"), c("m","cih"))
  renamed_args <-  list(c("m","v"),   c("m","v"), c("m","s"), c("m","s"), c("shape1","shape2"), c("m","cih"), c("m","cih"))
  parms <- c()
  for(i in 1:length(possible_args))
  {
    res <- match(possible_args[[i]], nms)
    if(!any(is.na(res)))
    {
      if(length(parms)>0) stop("Multiple arguments matched")
      parms <- element[res]
      names(parms) <- renamed_args[[i]]
    }
  }
  if(length(parms)==0) stop("No valid parameter specification")
  ##IF (m,s), apply the method of moments. Note that for normal it is a stupid tail chasing but OK!
  if(names(parms)[1]=='m')
  {
    m <- parms[[1]]
    if(names(parms)[2]=="cih")
    {
      e$parms <- inv_mean_quantile(element$type, m, q=parms[[2]], p=0.975)
      e$moments <- c(m=m, v=0)
      e$moments[2] <- moments(element$type, e$parms)[[2]]
    }
    else
    {
      v <- parms[[2]]
      if(names(parms)[2]=="s") v <- v^2
      e$moments <- c(m=m,v=v)
      e$parms <- inv_moments(element$type, list(m,v))
    }
  }
  else
  {
    e$parms <- parms
    e$moments <- moments(element$type, parms)
  }

  e
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
process_evidence <- function(evidence)
{
  if(is.null(evidence$prev)) {stop("evidence object must have a prev (prevalence) member")}
  if(is.null(evidence$prev$type)) {evidence$prev$type<-"beta"; message("Assuming prev has a beta distribution") }
  evidence$prev <- process_evidence_element(evidence$prev)
  
  if(is.null(evidence$cstat)) {stop("evidence object must have a cstat (c-statistic) member")}
  if(is.null(evidence$cstat$type)) {evidence$cstat$type<-"beta"; message("Assuming cstat has a beta distribution") }
  evidence$cstat <- process_evidence_element(evidence$cstat)
  
  nms <- names(evidence)
  possible_args <- list(c("cal_mean","cal_slp"), c("cal_oe","cal_slp"), c("cal_int","cal_slp"))
  renamed_args <-  list(c("cal_mean","cal_slp"), c("cal_oe","cal_slp"), c("cal_int","cal_slp"))
  cal_parms <- c()
  for(i in 1:length(possible_args))
  {
    res <- match(possible_args[[i]], nms)
    if(!any(is.na(res)))
    {
      if(length(cal_parms)>0) stop("Multiple arguments matched")
      cal_parms <- evidence[res]
      names(cal_parms) <- renamed_args[[i]]
    }
  }
  if(length(cal_parms)==0) stop("No valid parameter specification for calibration (calibration slope AND at lease one of intercept, mean calibration, or O/E ratio are required)")
  
  cal_mean <- match(c("cal_mean"),names(cal_parms))
  cal_int <- match(c("cal_int"),names(cal_parms))
  cal_slp <- match(c("cal_slp"),names(cal_parms))
  cal_oe <- match(c("cal_oe"),names(cal_parms))
  
  if(!is.na(cal_mean))
  {
    if(is.null(cal_parms[[cal_mean]]$type)) { cal_parms[[cal_mean]]$type<-"norm"; message("Assuming normal distirbution for calibration mean")}
    evidence$cal_mean <- process_evidence_element(cal_parms[[cal_mean]])
  }
  if(!is.na(cal_int))
  {
    if(is.null(cal_parms[[cal_int]]$type)) { cal_parms[[cal_int]]$type<-"norm"; message("Assuming normal distirbution for calibration interccept")}
    evidence$cal_int <- process_evidence_element(cal_parms[[cal_int]])
  }
  if(!is.na(cal_slp))
  {
    if(is.null(cal_parms[[cal_slp]]$type)) { cal_parms[[cal_slp]]$type<-"norm"; message("Assuming normal distirbution for calibration slope")}
    evidence$cal_slp <- process_evidence_element(cal_parms[[cal_slp]])
  }
  if(!is.na(cal_oe))
  {
    if(is.na(cal_parms[[cal_oe]]$type)) { cal_parms[[cal_oe]]$type<-"norm"; message("Assuming normal distirbution for O/E ratio")}
    evidence$cal_oe <- process_evidence_element(cal_parms[[cal_oe]])
  }
  
  evidence
}

#'Plots Calibration Instability from Simulated Calibration Curves
#'
#'@description Simulates calibration curves based on given method, and uses plot to visualize calibration instability.
#'@param N Number of observations to simulate in each sample
#'@param sample Data frame with columns:
#'  dist_type: distribution type
#'  dist_parm1: first distribution parameter (e.g. mean, alpha, shape1)
#'  dist_parm2: second distribution parameter (e.g. sd, beta, shape2)
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'@param method One of loess or line, on default is loess
#'@param X Vector of predicted probabilities, on default is 0.01 to 0.99
#'@return Plot of simulated calibration curves
#'@examples
#'sample <- data.frame(
#' dist_type = rep("beta", 3),
#' dist_parm1 = c(1,2,3),
#' dist_parm2 = c(3,4,5),
#' cal_int = c(0, 0.05, 0.1),
#' cal_slp = c(1, 0.9, 0.8))
#'plot_cal_instability(N=200, sample=sample)
#'@export
plot_cal_instability <- function(N, sample, method="loess", X=(1:99)/100)
{
  out <- matrix(NA, nrow=nrow(sample), ncol=length(X))
  for(i in 1:nrow(sample))
  {
    p <- do.call(paste0("r",sample$dist_type[i]), list(N, sample$dist_parm1[i], sample$dist_parm2[i]))
    Y <- rbinom(N, 1, p)
    pi <- expit((logit(p)-sample$cal_int[i])/sample$cal_slp[i])
    if(method=="loess")
    {
      fcl <- loess(Y~pi)
      out[i,] <- predict(fcl, X)
    }
    if(method=="line")
    {
      reg <- glm(Y~logit(pi), family=binomial())
      out[i,] <- predict(reg, newdata=data.frame(pi=X), type='response')
    }
  }
  
  ci <- apply(out, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
  ylim <- c(0,1) #c(min(out, na.rm=TRUE), max(out, na.rm=TRUE)) 
  plot(X, out[1,], ylim=ylim, xlab="Simulated predicted risks", ylab="Simulated observede risks", type='l', col='grey')
  for(i in 2:nrow(out))
  {
    lines(X, out[i,], col='grey')
  }
  #lines(X, ci[1,])
  #lines(X, ci[2,])
  lines(c(0,1),c(0,1))
}

#'Plots Calibration Distance from Simulation Curves
#'
#'@description simulates calibration curves based on given method, and uses plot to visualize calibration distance (difference between predicted and observed)
#'@param N Number of observations to simulate in each sample
#'@param sample Data frame with columns:
#'  dist_type: distribution type
#'  dist_parm1: first distribution parameter (e.g. mean, alpha, shape1)
#'  dist_parm2: second distribution parameter (e.g. sd, beta, shape2)
#'  cal_int: calibration intercept
#'  cal_slp: calibration slope
#'@param method One of loess or line, on default is loess
#'@param X Vector of predicted probabilities, on default is 0.01 to 0.99
#'@return Plot of simulated calibration curves
#'@examples
#'sample <- data.frame(
#' dist_type = rep("beta", 3),
#' dist_parm1 = c(1,2,3),
#' dist_parm2 = c(3,4,5),
#' cal_int = c(0, 0.05, 0.1),
#' cal_slp = c(1, 0.9, 0.8))
#'plot_cal_distance(N=200, sample=sample)
#'@export
plot_cal_distance <- function(N, sample, method="loess", X=(1:99)/100)
{
  out <- matrix(NA, nrow=nrow(sample), ncol=length(X))
  
  for(i in 1:nrow(sample))
  {
    p <- do.call(paste0("r",sample$dist_type[i]), list(N, sample$dist_parm1[i], sample$dist_parm2[i]))
    Y <- rbinom(N, 1, p)
    pi <- expit((logit(p)-sample$cal_int[i])/sample$cal_slp[i])
    
    if(method=="loess")
    {
      fcl <- loess(Y~pi)
      out[i,] <- predict(fcl, X)-expit(logit(X)*sample$cal_slp[i]+sample$cal_int[i])
    }
    if(method=="line")
    {
      reg <- glm(Y~logit(pi), family=binomial())
      out[i,] <- predict(reg, newdata=data.frame(pi=X), type="response")-expit(logit(X)*sample$cal_slp[i]+sample$cal_int[i])
    }
  }
  
  ci <- apply(out, 2, quantile, c(0.025, 0.975), na.rm=TRUE)
  ylim <- c(-1,1) #c(min(out, na.rm=TRUE), max(out, na.rm=TRUE)) 
  for(i in 1:nrow(out))
  {
    if(i==1)
    {
      plot(X, out[i,], ylim=ylim, xlab="Simulated predicted risks", ylab="Simulated calibration error", type='l', col='grey')
    }
    else
    {
      lines(X, out[i,], col='grey')
    }
  }
  lines(X, ci[1,], lty=3)
  lines(X, ci[2,], lty=3)
  lines(c(0,1),c(0,0))
}
