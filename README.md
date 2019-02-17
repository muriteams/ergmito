
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ergmito: Estimation of Little ‘ERGMs’ using exact likelihood

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/ergmito)](https://cran.r-project.org/package=ergmito)
[![Travis build
status](https://travis-ci.org/muriteams/ergmito.svg?branch=master)](https://travis-ci.org/muriteams/ergmito)
[![AppVeyor Build
status](https://ci.appveyor.com/api/projects/status/nl1irakr2g6y6w03?svg=true)](https://ci.appveyor.com/project/gvegayon/ergmito)
[![codecov](https://codecov.io/gh/muriteams/ergmito/branch/master/graph/badge.svg)](https://codecov.io/gh/muriteams/ergmito)

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("muriteams/ergmito")
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
#> [1,] -6.557076
ergm.exact(ans_ergm$coef, model)
#>      [,1]
#> [1,]  NaN

summary(ans_ergmito)
#> $coefs
#>           Estimate Std. Error     z value  Pr(>|z|)
#> edges   -0.3203616   1.355221 -0.23639058 0.8131296
#> mutual   2.3390394   2.096539  1.11566695 0.2645647
#> balance -9.3678948 100.942647 -0.09280413 0.9260592
#> 
#> $aic
#> [1] 19.11415
#> 
#> $bic
#> [1] 20.56887
#> 
#> $model
#> [1] "net ~ edges + mutual + balance"
#> 
#> attr(,"class")
#> [1] "ergmito_summary"
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
#>          Estimate Std. Error MCMC % z value Pr(>|z|)    
#> edges    0.006341   1.209717      0   0.005    0.996    
#> mutual  20.682872         NA     NA      NA       NA    
#> balance      -Inf   0.000000      0    -Inf   <1e-04 ***
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
#> $coefs
#>         Estimate Std. Error    z value  Pr(>|z|)
#> edges  -1.098612   0.912871 -1.2034693 0.2287948
#> mutual  1.098612   1.825742  0.6017347 0.5473507
#> 
#> $aic
#> [1] 18.90944
#> 
#> $bic
#> [1] 19.87925
#> 
#> $model
#> [1] "net ~ edges + mutual"
#> 
#> attr(,"class")
#> [1] "ergmito_summary"
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
