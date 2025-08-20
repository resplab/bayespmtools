
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Bayespmtools

<!-- badges: start -->
<!-- badges: end -->

The goal of Bayespmtools is to enable Bayesian sample size and precision
calculations for external validation of risk prediction models.

For the details of the methodology, please refer to the accompanying
paper: <https://arxiv.org/abs/2504.15923>

``` r
#Specify evidence:
evidence <- list(
  prev=list(type="beta", mean=0.427966984132821, sd=0.0295397309129426),
  cstat=list(type="logitnorm", mean=0.760628336908955, sd=0.00635806041351944),
  cal_mean=list(type="norm", mean=-0.00934717199436785, sd=0.124517605045825),
  cal_slp=list(type="norm", mean=0.995017759715243, sd = 0.0237278675967507))

#Specifying targets
#eciw=x indicates desired expected CI Width of x.
#qciw=c(a,b) indicates desired assurance CI Width of x at assurance level y.
targets <- list(eciw.cstat=0.1,
  eciw.cal_oe=0.22,
  eciw.cal_slp=0.30,
  qciw.cstat=c(0.1, 0.9),
  qciw.cal_oe=c(0.22, 0.9),
  qciw.cal_slp=c(0.30,0.9),
  assurance.nb=0.9)

library(bayespmtools)

#Main function call
res <- bpm_valsamp(evidence=evidence, #Evidence as a list
  dist_type="logitnorm", #Distribution type for calibrated risks
  method="sample", #Sample based or tw-level ("2s") method
  targets=targets, #Targets (as specified above)
  n_sim=100, #Number of Monte Carlo simulations
  threshold=0.2) #Risk threshold for NB VoI calculations
#> Processing evidence...
#> Generating Monte Carlo sample...
#> Imputing correlation...
#> Based on effective sample size: 281
#> Infering calibration intercept...
#> Computing CI sample size...
#> Computing se/sp...
#> VoI / NB assuraance...

print(res$N)
#>   eciw.cstat  eciw.cal_oe eciw.cal_slp   qciw.cstat  qciw.cal_oe qciw.cal_slp assurance.nb 
#>          356          405          910          405          524         1177          306
```
