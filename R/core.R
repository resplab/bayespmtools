
logit <- function(x) {log(x/(1-x))}

expit <- function(x) {1/(1+exp(-x))}




#'@export
moments <- function(type, parms) {
  # Check if the input values are valid
  if(type=="norm")
  {
    return(c(m=parms[1], v=sqrt(parms[2])))
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
    #TODO
  }
  if(type=="norm")
  {
    res <- uniroot(function(x) {pnorm(q,m,x)-p}, interval=c(0.0001,10))
    out <- c(mean=m, sd=res$root)
  }
  
  out
}



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





#'@export
calc_se_sp <- function(dist_type, dist_parms, cal_int, cal_slp, threshold, prev=NULL)
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










#'@export
infer_correlation <- function(dist_type, dist_parms, cal_int, cal_slp, n, n_sim)
{
  require(pROC)
  require(mcmapper)

  out <- matrix(NA, nrow=n_sim, ncol=5)
  colnames(out) <- c("prev", "cstat", "cal_mean", "cal_int", "cal_slp")

  for(i in 1:n_sim)
  {
    p <- do.call(paste0("r",dist_type), args=as.list(c(n,unname(dist_parms))))
    Y <- rbinom(n,1,p)
    pi <- expit((logit(p)-cal_int)/cal_slp)
    df <- cbind(pi=pi,Y=Y)
    out[i,]<-c(mean(df[2]),
               roc(df[,2]~df[,1], quiet=TRUE)$auc,
               mean(df[2]-df[1]),
               coef(glm(df[,2]~logit(df[,1]), family=binomial(link="logit")))[1:2])
  }

  cor(out)
}









#' @export
express_evidence <- function(evidence, round_digits=3)
{
  out <- "";
  for(i in 1:length(evidence))
  {
    nm <- names(evidence)[i]
    item <- evidence[[i]]
    out <- paste0(out, nm, "~", item$type, "(", round(item$parms[1], round_digits), ",", round(item$parms[2], round_digits),")\n")
  }
  
  paste(out, collapse="\n")
}






#N can be vector!
#' @export
calc_riley_vars <- function(N, parms)
{
  prev <- parms$prev
  C <- parms$cstat
  dist_type <- parms$dist_type
  dist_parms <- c(parms$dist_parm1,parms$dist_parm2)
  cal_int <- parms$cal_int
  cal_slp <- parms$cal_slp

  if(dist_type=="logitnorm")
  {
    #tryCatch({
      E_pi <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(expit((logit(x)-cal_int)/cal_slp))},0,1)$value
      I_a <- integrate(function(x){ mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(x*(1-x))}, 0, 1)$value
      I_b <- integrate(function(x){ mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)^2*(x*(1-x)))}, 0, 1)$value
      I_ab <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)*(x*(1-x)))}, 0, 1)$value
    #}, error=function(e){bad_place$parms <<- parms})
  }
  if(dist_type=="beta")
  {
    E_pi <- integrate(function(x){x*dbeta(expit(logit(x)*cal_slp+cal_int),dist_parms[1],dist_parms[2])}, 0, 1,intercept=cal_int, slope=cal_slp)$value
    I_a <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*(exp(cal_int+x*cal_slp)/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
    I_b <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*((x^2*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
    I_ab <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*((x*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
  }
  if(dist_type=="probitnorm")
  {
    E_pi <- integrate(function(x){x*mcmapper::dprobitnorm(expit(logit(x)*cal_slp+cal_int),dist_parms[1],dist_parms[2])}, 0, 1,intercept=cal_int, slope=cal_slp)$value
    I_a <- integrate(function(x){mcmapper::dprobitnorm(expit(logit(x)), dist_parms[1], dist_parms[2])*(exp(cal_int+x*cal_slp)/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
    I_b <- integrate(function(x){mcmapper::dprobitnorm(expit(logit(x)), dist_parms[1], dist_parms[2])*((x^2*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
    I_ab <- integrate(function(x){mcmapper::dprobitnorm(expit(logit(x)), dist_parms[1], dist_parms[2])*((x*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
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





#' @export
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
      I_a <- integrate(function(x){ mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(x*(1-x))}, 0, 1)$value
      I_b <- integrate(function(x){ mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)^2*(x*(1-x)))}, 0, 1)$value
      I_ab <- integrate(function(x){mcmapper::dlogitnorm(x,dist_parms[1],dist_parms[2])*(((logit(x)-cal_int)/cal_slp)*(x*(1-x)))}, 0, 1)$value
      #}, error=function(e){bad_place$parms <<- parms})
    }
    if(dist_type=="beta")
    {
      #E_pi <- integrate(function(x){x*dbeta(expit(logit(x)*cal_slp+cal_int),dist_parms[1],dist_parms[2])}, 0, 1,intercept=cal_int, slope=cal_slp)$value
      I_a <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*(exp(cal_int+x*cal_slp)/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
      I_b <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*((x^2*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
      I_ab <- integrate(function(x){mcmapper::dbeta(expit(logit(x)), dist_parms[1], dist_parms[2])*((x*exp(cal_int+x*cal_slp))/(1+exp(cal_int+x*cal_slp))^2)}, 0, 1)$value
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






#N can be vector
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



#N can be vector
#'@export
calc_ciw_sample <- function(N, parms)
{
  require(fastLogisticRegressionWrap)
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
    reg <- fast_logistic_regression(logit_pi[1:n,],Y[1:n],do_inference_on_var="all")
    se_cal_int[i] <- reg$se[1]
    se_cal_slp[i] <- reg$se[2]

    #tmp <- ci.auc(Y[1:n]~logit_pi[1:n], quiet=TRUE)
    tmp <- ci.auc(Y[1:n]~logit_pi[1:n,2], quiet=TRUE)
    se_cstat[i] <- (tmp[3]-tmp[2])/qnorm(0.975)
  }

  k <- 2*qnorm(0.975)
  list(cstat=k*se_cstat, cal_oe=k*se_cal_oe, cal_mean=k*se_cal_mean, cal_int=k*se_cal_int, cal_slp=k*se_cal_slp)
}







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


#Convext monotonically decreasing root finder
find_n_mean <- function(target, N, ciws, decreasing=T, convex=T)
{
  require(cobs)
  
  X <- rep(N, each=nrow(ciws))
  Y <- as.vector(ciws)
  
  constraint <- c(ifelse(decreasing,"decreasing",""),ifelse(convex,"convex",""))
  
  S <- cobs(X,Y, constraint=c("decrease","convex"))
  #plot(predict(S, seq(min(N),max(N),length.out=100)), type='l')
  
  f <- function(x) {predict(S, z=x)[,2] - target}
  round(uniroot(f, c(min(N),max(N)))$root,0)
}



#Convext monotonically decreasing root finder
find_n_quantile <- function(target, N, q, ciws)
{
  require(quantreg)
  
  X <- rep(N, each=nrow(ciws))
  Y <- as.vector(ciws)
  
  
  fit <- rqss(Y ~ qss(X, lambda = 1), tau=q)
  
  # plot(X, Y, main = "Non-parametric Quantile Regression", xlab = "X", ylab = "Y")
  # for (i in 1:length(tau)) {
  #   lines(sort(X), predict(fit, newdata = data.frame(X = sort(X)), tau = tau[i]), col = i + 1)
  # }
  # legend("topleft", legend = paste("tau =", tau), col = 2:(length(tau) + 1), lty = 1)
  
  f <- function(x) {predict(fit, newdata=data.frame(X=x)) - target}
  round(uniroot(f, c(min(N),max(N)))$root,0)
}






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



















# Creates a clean list including type, parms, moments
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


#' @export
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
    if(is.na(cal_parms[[cal_mean]]$type)) { cal_parms[[cal_mean]]$type<-"normal"; message("Assuming normal distirbution for calibration mean")}
    evidence$cal_mean <- process_evidence_element(cal_parms[[cal_mean]])
  }
  if(!is.na(cal_int))
  {
    if(is.na(cal_parms[[cal_int]]$type)) { cal_parms[[cal_int]]$type<-"normal"; message("Assuming normal distirbution for calibration mean")}
    evidence$cal_int <- process_evidence_element(cal_parms[[cal_int]])
  }
  if(!is.na(cal_slp))
  {
    if(is.na(cal_parms[[cal_slp]]$type)) { cal_parms[[cal_slp]]$type<-"normal"; message("Assuming normal distirbution for calibration mean")}
    evidence$cal_slp <- process_evidence_element(cal_parms[[cal_slp]])
  }
  if(!is.na(cal_oe))
  {
    if(is.na(cal_parms[[cal_oe]]$type)) { cal_parms[[cal_oe]]$type<-"normal"; message("Assuming normal distirbution for calibration mean")}
    evidence$cal_oe <- process_evidence_element(cal_parms[[cal_oe]])
  }
  
  evidence
}



#' @export
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





#' @export
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
