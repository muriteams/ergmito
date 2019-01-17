
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ergmito: Estimation of Little ‘ERGMs’ using exact likelihood

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/ergmito)](https://cran.r-project.org/package=ergmito)
[![Travis build
status](https://travis-ci.org/muriteams/lergm.svg?branch=master)](https://travis-ci.org/muriteams/lergm)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/avadsgb4a3tg20fd?svg=true)](https://ci.appveyor.com/project/gvegayon/lergm)
[![Coverage
status](https://codecov.io/gh/muriteams/lergm/branch/master/graph/badge.svg)](https://codecov.io/github/muriteams/lergm?branch=master)

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("muriteams/lergm")
```

## Example

An example from the manual

## When `ergm` is not enough

``` r
library(ergmito)
library(sna)

# Generating a small graph
set.seed(12)
n <- 4
net <- sna::rgraph(n, tprob = .7)
gplot(net)
```

<img src="man/figures/README-net1-1.png" width="100%" />

``` r
model <- net ~ edges + mutual + balance

library(ergm)
ans_ergmito <- ergmito(model)
ans_ergm  <- ergm(model)

# The ergmito should have a larger value
ergm.exact(ans_ergmito$coef, model)
#>           [,1]
#> [1,] -6.556996
ergm.exact(ans_ergm$coef, model)
#>      [,1]
#> [1,]  NaN

summary(ans_ergmito)
#> 
#> ERGMito estimates
#>               Length Class          Mode   
#> call           2     -none-         call   
#> coef           3     -none-         numeric
#> iterations     1     -none-         numeric
#> loglikelihood  1     -none-         numeric
#> covar          9     -none-         numeric
#> coef.init      3     -none-         numeric
#> formulae       7     ergmito_loglik list   
#> network       16     -none-         numeric
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
#>         Estimate Std. Error MCMC % z value Pr(>|z|)    
#> edges    0.00116    1.21668      0   0.001    0.999    
#> mutual  20.68287         NA     NA      NA       NA    
#> balance     -Inf    0.00000      0    -Inf   <1e-04 ***
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

Checking convergence diagnostics

``` r
plot(ans_ergmito)
```

<img src="man/figures/README-convergence-diag-1.png" width="100%" />

## Do we get the same?

``` r
# Generating a small graph
set.seed(12123)
n   <- 4
net <- sna::rgraph(n, tprob = .3)
gplot(net)
```

<img src="man/figures/README-net2-1.png" width="100%" />

``` r
model <- net ~ edges + mutual

library(ergm)
ans_ergmito <- ergmito(model)
ans_ergm  <- ergm(model, control = control.ergm(
  MCMC.effectiveSize = 4000,
  seed = 444)
  )

# The ergmito should have a larger value
ergm.exact(ans_ergmito$coef, model) > ergm.exact(ans_ergm$coef, model)
#>      [,1]
#> [1,] TRUE

summary(ans_ergmito)
#> 
#> ERGMito estimates
#>               Length Class          Mode   
#> call           2     -none-         call   
#> coef           2     -none-         numeric
#> iterations     1     -none-         numeric
#> loglikelihood  1     -none-         numeric
#> covar          4     -none-         numeric
#> coef.init      2     -none-         numeric
#> formulae       7     ergmito_loglik list   
#> network       16     -none-         numeric
summary(ans_ergm)
#> 
#> ==========================
#> Summary of model fit
#> ==========================
#> 
#> Formula:   net ~ edges + mutual
#> 
#> Iterations:  2 out of 20 
#> 
#> Monte Carlo MLE Results:
#>        Estimate Std. Error MCMC % z value Pr(>|z|)
#> edges   -1.1003     0.9041      0  -1.217    0.224
#> mutual   1.1008     1.8194      0   0.605    0.545
#> 
#>      Null Deviance: 16.64  on 12  degrees of freedom
#>  Residual Deviance: 14.90  on 10  degrees of freedom
#>  
#> AIC: 18.9    BIC: 19.87    (Smaller is better.)
```

# Similarity indices

<https://cran.r-project.org/web/packages/proxy/proxy.pdf>

A Survey of Binary Similarity and Distance Measures Seung-Seok Choi,
Sung-Hyuk Cha, Charles C. Tappert Department of Computer Science, Pace
University New York, US

# Contributing

Please note that the ‘ergmito’ project is released with a [Contributor
Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.
