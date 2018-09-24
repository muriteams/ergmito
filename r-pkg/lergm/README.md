
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lergm: Estimation of Little ‘ERGMs’ using exact likelihood

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/lergm)](https://cran.r-project.org/package=lergm)

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("USCCANA/social-smarts/r-pkg/lergm")
```

## Example

An example from the manual

``` r
library(lergm)

# Generating a small graph
set.seed(12)
n <- 4
net <- sna::rgraph(n, tprob = .7)

model <- net ~ edges + mutual + balance

library(ergm)
#> Loading required package: statnet.common
#> 
#> Attaching package: 'statnet.common'
#> The following object is masked from 'package:base':
#> 
#>     order
#> Loading required package: network
#> network: Classes for Relational Data
#> Version 1.13.0.1 created on 2015-08-31.
#> copyright (c) 2005, Carter T. Butts, University of California-Irvine
#>                     Mark S. Handcock, University of California -- Los Angeles
#>                     David R. Hunter, Penn State University
#>                     Martina Morris, University of Washington
#>                     Skye Bender-deMoll, University of Washington
#>  For citation information, type citation("network").
#>  Type help("network-package") to get started.
#> 
#> ergm: version 3.8.0, created on 2017-08-18
#> Copyright (c) 2017, Mark S. Handcock, University of California -- Los Angeles
#>                     David R. Hunter, Penn State University
#>                     Carter T. Butts, University of California -- Irvine
#>                     Steven M. Goodreau, University of Washington
#>                     Pavel N. Krivitsky, University of Wollongong
#>                     Martina Morris, University of Washington
#>                     with contributions from
#>                     Li Wang
#>                     Kirk Li, University of Washington
#>                     Skye Bender-deMoll, University of Washington
#> Based on "statnet" project software (statnet.org).
#> For license and citation information see statnet.org/attribution
#> or type citation("ergm").
#> NOTE: Versions before 3.6.1 had a bug in the implementation of the
#> bd() constriant which distorted the sampled distribution somewhat.
#> In addition, Sampson's Monks datasets had mislabeled vertices. See
#> the NEWS and the documentation for more details.
#> 
#> Attaching package: 'ergm'
#> The following objects are masked from 'package:statnet.common':
#> 
#>     colMeans.mcmc.list, sweep.mcmc.list
ans_lergm <- lergm(model)
ans_ergm  <- ergm(model)
#> Observed statistic(s) balance are at their smallest attainable values. Their coefficients will be fixed at -Inf.
#> Starting maximum likelihood estimation via MCMLE:
#> Iteration 1 of at most 20:
#> Optimizing with step length 1.
#> The log-likelihood improved by 0.0006522.
#> Step length converged once. Increasing MCMC sample size.
#> Iteration 2 of at most 20:
#> Optimizing with step length 1.
#> The log-likelihood improved by 0.0005636.
#> Step length converged twice. Stopping.
#> Evaluating log-likelihood at the estimate. Using 20 bridges: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 .
#> This model was fit using MCMC.  To examine model diagnostics and check for degeneracy, use the mcmc.diagnostics() function.

# The lergm should have a larger value
ergm.exact(ans_lergm$coef, model)
#>           [,1]
#> [1,] -6.557266
ergm.exact(ans_ergm$coef, model)
#>      [,1]
#> [1,]  NaN

summary(ans_lergm)
#> 
#> Little ERGM estimates
#> 
#> ==========================
#> Summary of model fit
#> ==========================
#> 
#> Formula:   net ~ edges + mutual + balance
#> 
#> Iterations:  101 out of 20 
#> 
#> Monte Carlo MLE Results:
#>         Estimate Std. Error MCMC % p-value
#> edges    -0.3195     1.9166     29   0.871
#> mutual    2.3360     2.9648     29   0.451
#> balance  -8.2994    83.7136     29   0.923
#> 
#>      Null Deviance: 16.64  on 12  degrees of freedom
#>  Residual Deviance: 13.11  on  9  degrees of freedom
#>  
#> AIC: 31.11    BIC: 35.48    (Smaller is better.)
summary(ans_ergm)
#> 
#> ==========================
#> Summary of model fit
#> ==========================
#> 
#> Formula:   net ~ edges + mutual + balance
#> 
#> Iterations:  2 out of 20 
#> 
#> Monte Carlo MLE Results:
#>          Estimate Std. Error MCMC % p-value    
#> edges    0.002931   1.239047      0   0.998    
#> mutual  20.682872         NA     NA      NA    
#> balance      -Inf   0.000000      0  <1e-04 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>      Null Deviance: 16.64  on 12  degrees of freedom
#>  Residual Deviance:   NaN  on  9  degrees of freedom
#>  
#> AIC: NaN    BIC: NaN    (Smaller is better.) 
#> 
#>  Warning: The following terms have infinite coefficient estimates:
#>   balance
```
