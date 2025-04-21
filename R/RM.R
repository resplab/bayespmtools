
find_n_RM <- function(sample, target_rules, target_metrics, target_values, N0=1000, base_parms=NULL)
{
  n <- length(target_rules)
  if(is.null(base_parms))
  {
    base_parms <- list()
    base_parms$dist_type <- sample[1,'dist_type']
    base_parms$prev <- mean(sample$prev)
    base_parms$cstat <- mean(sample$cstat)
    tmp <- mcmapper::mcmap(c(base_parms$prev, base_parms$cstat), type=dist_type)$valu
    base_parms$dist_parm1 <- tmp[1]; base_parms$dist_parm2 <- tmp[2]
    base_parms$cal_slp <- mean(sample$cal_slp)
    base_parms$cal_int <- mean(sample$cal_int)
  }

  base_ciws <- bayescpm:::calc_ciw_2s(N0, base_parms)
  working_Ns <- rep(N0,n)
  names(working_Ns) <- paste0(target_rules, ".", target_metrics)
  for(i in 1:n)
  {
    K <- base_ciws[[target_metrics[i]]]/target_values[[i]][1]
    working_Ns[i] <- working_Ns[i]*K^2
  }
  n_sim <- nrow(sample)
  steps <- working_Ns
  
  X <- Y <- matrix(NA, nrow = n_sim, ncol=n)
  colnames(X) <- names(working_Ns) 
  colnames(Y) <- names(working_Ns) 

  for(i in 1:n_sim)
  {
    y <- bayescpm:::calc_ciw_sample(round(working_Ns,0), sample[i,])
    
    for(j in 1:n)
    {
      if(target_rules[j]=="eciw")
      {
        working_Ns[j] <- working_Ns[j]+steps[j]/i*(y[[target_metrics[j]]][j]-target_values[[j]])
      }
      if(target_rules[j]=="qciw")
      {
        working_Ns[j] <- working_Ns[j]+steps[j]/i*ifelse(y[[target_metrics[j]]][j]>target_values[[j]][1], target_values[[j]][2], -(1-target_values[[j]][2]))
      }
      X[i,j] <- working_Ns[j]
      Y[i,j] <- y[[target_metrics[j]]][j]
    }
  }
  return(list(N=round(working_Ns), trace=cbind(X,Y)))
}




