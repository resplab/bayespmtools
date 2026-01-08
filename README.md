
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Bayespmtools

<!-- badges: start -->

<!-- badges: end -->

The goal of Bayespmtools is to enable Bayesian sample size and precision
calculations for external validation of risk prediction models.

For the details of the methodology, please refer to the accompanying
paper: <https://doi.org/10.1002/sim.70389>

``` r
#Specify evidence:
evidence <- list(
  prev ~ beta(116, 155),           # Outcome prevalence
  cstat ~ beta(3628, 1139),        # C-statistic
  cal_mean ~ norm(-0.009, 0.125),  # Mean calibration error
  cal_slp ~ norm(0.995, 0.024)     # Calibration slope
)

#Specifying targets
#eciw=x indicates desired expected CI Width of x.
#qciw=c(a,b) indicates desired assurance CI Width of b at assurance level a.
#voi.nb indicates targeting a given EVSI/EVPI ratio
targets <- list(
  eciw.cstat = 0.1,             # Expected CI width for c-statistic
  eciw.cal_oe = 0.22,           # Expected CI width for O/E ratio
  eciw.cal_slp = 0.30,          # Expected CI width for calibration slope
  qciw.cal_slp = c(0.90, 0.35),  # 90% assurance that CI width ≤ 0.35
  voi.nb = 0.90
)

library(bayespmtools)

#Main function call
res <- bpm_valsamp(
  evidence = evidence,
  targets = targets,
  n_sim = 1000,           # Number of Monte Carlo simulations
  threshold = 0.2         # Risk threshold for net benefit calculations
)
#> Processing evidence...
#> Generating Monte Carlo sample...
#> Imputing correlation, based on effective sample size: 272 ...
#> Infering calibration intercept...
#> Computing CI sample size...
#> Computing se/sp...
#> VoI / NB assuraance...

print(res$results)
#>   eciw.cstat  eciw.cal_oe eciw.cal_slp qciw.cal_slp       voi.nb 
#>          342          408         1067          869          717
```
