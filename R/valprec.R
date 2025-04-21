#' @export
bpm_valprec <- function(N, evidence, targets=list(fciw.cstat=T, fciw.cal_slp=T, eciw.cstat=T, eciw.cal_oe=T, qciw.cal_oe=c(0.9), assurance.nb=T, voi.nb=T),
                     n_sim=1000, 
                     method="sample", 
                     threshold=NULL, 
                     dist_type="logitnorm",
                     impute_cor=TRUE,
                     ex_args=NULL)
{
  out <- list(N=N)
  
  tmp <- as.data.frame(lapply(names(targets), strsplit, "[.]"))
  target_rules <- unname(unlist(tmp[1,]))
  target_metrics <- unname(unlist(tmp[2,]))
  target_values <- (targets)
  
  if(is.function(ex_args$f_progress))
  {
    f_progress <- ex_args$f_progress
  }else
  {
    f_progress <- base::message
  }
  
  
  ##Step 1: Process evidence
  if(!is.matrix(evidence) & !is.data.frame(evidence))
  {
    f_progress("Processing evidence...")
    evidence <- process_evidence(evidence)
    out$evidence <- evidence
    
    ##Generate base parms
    base <- list()
    base$dist_type <- dist_type
    base$prev <- evidence$prev$moments[[1]]
    base$cstat <- evidence$cstat$moments[[1]]
    tmp <- mcmapper::mcmap(c(base$prev, base$cstat), type=dist_type)$valu
    base$dist_parm1 <- tmp[1]; base$dist_parm2 <- tmp[2]
    base$cal_slp <- evidence$cal_slp$moments[[1]]
    #Cal intercept is not provided. Need to derive it for base params for correlation induction
    if(is.na(match("cal_int", names(evidence))))
    {
      base$cal_int <- infer_cal_int_from_mean(base$dist_type, c(base$dist_parm1, base$dist_parm2), cal_mean=evidence$cal_mean$moments[[1]], cal_slp=base$cal_slp, prev=base$prev)
    }else
    {
      base$cal_int <- evidence$cal_int$moments[[1]]
    }
    
    
    
    #Step 2: generate sample of marginals
    f_progress("Generating Monte Carlo sample...")
    sample <- NULL
    for(element in evidence)
    {
      sample <- cbind(sample, do.call(paste0("r",element$type), args=as.list(c(n=n_sim,element$parms))))
    }
    colnames(sample) <- names(evidence)
    
    #TODO: replace bad c-statistic values
    n_bads <- 0
    repeat{
      bads <- which(sample[,'cstat']<0.51 |  sample[,'cstat']>0.99)
      if(length(bads)==0) break;
      n_bads <- n_bads+length(bads)
      subsample <- NULL
      for(element in evidence)
      {
        subsample <- cbind(subsample, do.call(paste0("r",element$type), args=as.list(c(n=length(bads),element$parms))))
      }
      sample[bads,] <- subsample
    }
    if(n_bads>0) warning(paste("in step 'Generating MOnte Carlo sample' - ", n_bads, "observations were replaced due to bad value of c-statistic."))
    
    #Step 3: induce correlation (if asked)
    if(impute_cor)
    {
      f_progress("Imputing correlation...")
      eff_n <- round(evidence$prev$moments[[1]]*(1-evidence$prev$moments[[1]])/evidence$prev$moments[[2]],0)
      f_progress(paste("Based on effective sample size:", eff_n))
      base_cor <- infer_correlation(base$dist_type, c(base$dist_parm1,base$dist_parm2), base$cal_int, base$cal_slp, eff_n, 1000)
      good_rows <- match(colnames(sample), colnames(base_cor))
      base_cor <- base_cor[good_rows,good_rows]
      sample <- mc2d::cornode(sample, target=base_cor)
    }
    
    sample <- as.data.frame(sample)
    sample$dist_type <- dist_type
    sample$dist_parm1 <- 0
    sample$dist_parm2 <- 0
  }else
  {
    sample <- evidence
    base <- as.data.frame(t(colMeans(sample[,which(colnames(sample)!="dist_type")])))
    base$dist_type <- sample[1,'dist_type']
  }
  
  
  
  #Step 4: if intercept is missing, impute it for the whole sample
  f_progress("Infering calibration intercept...")
  
  if(is.na(match("cal_int",colnames(sample))))
  {
    sample$cal_int <- NA
    for(i in 1:nrow(sample))
    {
      prev <- unname(sample[i,'prev'])
      cstat<- unname(sample[i,'cstat'])
      cal_mean <- unname(sample[i,'cal_mean']) #TODO
      cal_slp <- unname(sample[i,'cal_slp'])
      
      parms <- mcmap(c(prev, cstat), dist_type)$value
      sample$dist_parm1[i] <- parms[1]
      sample$dist_parm2[i] <- parms[2]
      
      cal_int <- infer_cal_int_from_mean(dist_type=dist_type, dist_parms=parms, cal_mean=cal_mean, cal_slp=cal_slp, prev=prev)
      
      sample[i,'cal_int'] <- cal_int
    }
  }
  
  
  
  # Step 5: Freq & Bayesian Riley
  
  fciws <- which(target_rules=="fciw" & !isFALSE(target_values))
  if(length(fciws>0)) #Frequentist CIWs
  {
    fv <- calc_riley_vars(N, parms=base)
    
    for(item in names(fv))
    {
      if(is.na(match(item,target_metrics[fciws])))
      {
        fv[[item]] <- NULL
      }
    }
  
    for(item in names(fv))
    {
      fv[[item]] <- sqrt(fv[[item]])*2*qnorm(0.975)
    }
    
    out$fciw <- fv
  }
  
  bciws <- which((target_rules=="eciw" | target_rules=="qciw") & !isFALSE(target_values))
  if(length(bciws)>0)
  {
    f_progress("Computing CI widths...")
    
    ciws <- calc_ciw_mc(N, sample, method=method)
    
    for(i in 1:length(bciws))
    {
      if(target_rules[bciws[i]]=="eciw") 
      {
        out$eciw[[target_metrics[bciws[i]]]] <- colMeans(ciws[[target_metrics[bciws[i]]]])
      }
      if(target_rules[bciws[i]]=="qciw")
      {
        out$qciw[[target_metrics[bciws[i]]]] <- apply(ciws[[target_metrics[bciws[i]]]], 2, quantile, target_values[[bciws[i]]][2])
      }
    }
    
    for(item in names(ciws))
    {
      if(is.na(match(item, target_metrics[bciws])))
      {
        ciws[[item]] <- NULL
      }
    }
    out$ciws <- ciws
  }
  
 
  b_assurance <- sum((target_rules=="assurance" & !isFALSE(target_values)))
  b_voi <- sum((target_rules=="voi" & !isFALSE(target_values)))
  if(b_assurance | b_voi)
  {
    if(is.null(threshold)) stop("NB-related stuff was requested but threshold is not specified") #Todo: move earlier to avoid this late error
# Step 6: Calculate se and sp
    if(is.na(match("se", colnames(sample))))
    {
      f_progress("Computing se/sp...")
      sample$sp <- sample$se <- NA
      for(i in 1:nrow(sample))
      {
        se_sp <- calc_se_sp(sample$dist_type[i],
                            c(sample$dist_parm1[i], sample$dist_parm2[i]),
                            sample$cal_int[i],
                            sample$cal_slp[i],
                            threshold,
                            sample$prev[i]
        )
        sample[i,c('se','sp')] <- se_sp
      }
    }
    
    
#Step 7: VoI
    f_progress("VoI and NB assuraance...")
    
    #require(evsiexval)
    #res <- evsiexval::EVSI_gf(sample[,c('prev','se','sp')], future_sample_sizes=N,  ignore_prior=TRUE, z=threshold)
    
    tnb1 <- sample$prev*sample$se - (1-sample$prev)*(1-sample$sp)*threshold/(1-threshold)
    tnb2 <- sample$prev - (1-sample$prev)*threshold/(1-threshold)
    tnbs <- cbind(0,tnb1,tnb2)
    maxnbs <- apply(tnbs, 1, max)
    maxenb <- max(colMeans(tnbs))
    evsi <- assurance <- rep(0, length(N))
    evpi <- mean(maxnbs)-maxenb
    assurance0 <- mean(apply(tnbs, 1, which.max)==which.max(colMeans(tnbs)))
    n_rep <- 1
    for(I in 1:n_rep)
    {
      for(i in 1:length(N))
      {
        n <- N[i]
        nd <- rbinom(rep(n,n_sim), n, sample$prev)
        ntp <- rbinom(rep(n,n_sim), nd, sample$se)
        nfp <- rbinom(rep(n,n_sim), n-nd, 1-sample$sp)
        
        nb1 <- ntp/n - nfp/n*threshold/(1-threshold)
        nb2 <- nd/n - (1-nd/n)*threshold/(1-threshold)
        
        winners <- apply(cbind(0,nb1,nb2), 1, which.max)
        winnertnbs <- tnbs[cbind(1:n_sim,winners)]
        evsi[i] <- evsi[i]+mean(winnertnbs)-maxenb
        assurance[i] <- assurance[i]+mean(winnertnbs==maxnbs)
      }
    }
    if(b_voi)
    {
      out$voi$evpi <- evpi #res$EVPI
      out$voi$evsi <- evsi/n_rep #res$EVSI
    }
    if(b_assurance)
    {
      out$assurance$assurance0 <- assurance0
      out$assurance$assurance <- assurance/n_rep #res$EVSIp
    }
  }
  
  out$sample <- sample
  out
}

