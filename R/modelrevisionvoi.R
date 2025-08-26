#'@export
bpm_modelrevisionvoi <- function(N, sample, threshold, method=c("int","int+slp")[2], ex_args=NULL)
{
  out <- list(N=N)
  
  if(is.function(ex_args$f_progress))
  {
    f_progress <- ex_args$f_progress
  }else
  {
    f_progress <- base::message
  }
  
  f_progress("VoI and NB assuraance for model revision...")
  
  enb1 <- sample$prev*sample$se - (1-sample$prev)*(1-sample$sp)*threshold/(1-threshold)
  enb2 <- sample$prev - (1-sample$prev)*threshold/(1-threshold)
  enbs <- cbind(0,enb1,enb2)
  maxenbs <- mean(apply(enbs, 1, max))
  
  n_sim <- nrow(sample)
  assurance <- evsi <- rep(0, length(N))
  
  for(i in 1:length(N))
  {
    n <- N[i]
    
    for(j in 1:nrow(sample))
    {
      this_row <- sample[j,]
      p <- do.call(paste0("r",this_row$dist_type), as.list(unname(c(n,this_row$dist_parm1, this_row$dist_parm2))))
      pi <- expit((logit(p)-this_row$cal_int)/this_row$cal_slp)
      logit_pi <- cbind(1,logit(pi))
      Y <- rbinom(n, 1, p)
        
      reg <- fast_logistic_regression(logit_pi[1:n,],Y[1:n])
      est_cal_int <- reg$coefficients[1]
      est_cal_slp <- reg$coefficients[2]
      this_row$cal_int <- this_row$cal_int-(this_row$cal_slp/est_cal_slp)*est_cal_int
      this_row$cal_slp <- this_row$cal_slp/est_cal_slp
            
      updated_se_sp <- calc_se_sp(this_row$dist_type, 
                                  c(this_row$dist_parm1, this_row$dist_parm2), 
                                  this_row$cal_int, 
                                  this_row$cal_slp, 
                                  threshold, this_row$prev)
      
      sample$se[j] <- updated_se_sp[1]
      sample$sp[j] <- updated_se_sp[2]
      
      #cat('.')
    }
    
    unb1 <- sample$prev*sample$se - (1-sample$prev)*(1-sample$sp)*threshold/(1-threshold)
    unb2 <- sample$prev - (1-sample$prev)*threshold/(1-threshold)
    unbs <- cbind(0,unb1,unb2)
    maxunbs <- apply(unbs, 1, max)
    
    evsi[i] <- evsi[i]+mean(maxunbs)
    assurance[i] <- NA
  }
  
  out$evsi <- evsi-maxenbs
  out$assurance <- assurance
  out
}
