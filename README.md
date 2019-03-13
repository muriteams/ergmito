
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

This R package, which has been developed on top of the amazing work that
the [Statnet](https://github.com/statnet) team has done, implements
estimation and simulation methods for Exponential Random Graph Models of
small networks, in particular, less than 7 nodes. In the case of small
networks, the calculation of the likelihood of ERGMs becomes
computationally feasible, which allows us avoiding approximations and do
exact calculation, ultimately obtaining MLEs directly.

## Support

This material is based upon work support by, or in part by, the U.S.
Army Research Laboratory and the U.S. Army Research Office under grant
number W911NF-15-1-0577

Computation for the work described in this paper was supported by the
University of Southern California’s Center for High-Performance
Computing (hpcc.usc.edu).

## Installation

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
model <- net ~ edges + mutual + ctriad

library(ergm)
ans_ergmito <- ergmito(model)
ans_ergm  <- ergm(model)

# The ergmito should have a larger value when computing exact loglikelihood
ergm.exact(ans_ergmito$coef, model) >
  ergm.exact(ans_ergm$coef, model)
#>      [,1]
#> [1,] TRUE

summary(ans_ergmito)
#> $coefs
#>           Estimate Std. Error    z value  Pr(>|z|)
#> edges    0.5557799   1.772881  0.3134896 0.7539087
#> mutual   1.0722314   1.796543  0.5968304 0.5506206
#> ctriple -0.6000589   1.293807 -0.4637930 0.6427960
#> 
#> $aic
#> [1] 20.69027
#> 
#> $bic
#> [1] 22.14499
#> 
#> $model
#> [1] "net ~ edges + mutual + ctriad"
#> 
#> attr(,"class")
#> [1] "ergmito_summary"
summary(ans_ergm)
#> 
#> ==========================
#> Summary of model fit
#> ==========================
#> 
#> Formula:   net ~ edges + mutual + ctriad
#> 
#> Iterations:  2 out of 20 
#> 
#> Monte Carlo MLE Results:
#>         Estimate Std. Error MCMC % z value Pr(>|z|)
#> edges     0.5736     1.7862      0   0.321    0.748
#> mutual    1.0697     1.7993      0   0.594    0.552
#> ctriple  -0.5940     1.2914      0  -0.460    0.646
#> 
#>      Null Deviance: 16.64  on 12  degrees of freedom
#>  Residual Deviance: 14.66  on  9  degrees of freedom
#>  
#> AIC: 20.66    BIC: 22.11    (Smaller is better.)
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

## Estimating data with known parameters

The following example shows the estimation of a dataset that is included
in the package, `fivenets`. This set of five networks was generated
using the `new_rergmito` function which allows creating a function to
draw random ERGMs with a fixed set of parameters, in this case, `edges =
-2.0` and `nodematch("female") = 2.0`

``` r
data(fivenets)

model1 <- ergmito(fivenets ~ edges + nodematch("female"))

summary(model1) # This data has know parameters equal to -2.0 and 2.0
#> $coefs
#>                   Estimate Std. Error   z value    Pr(>|z|)
#> edges            -1.704748  0.5435573 -3.136280 0.001711055
#> nodematch.female  1.586965  0.6430475  2.467882 0.013591530
#> 
#> $aic
#> [1] 73.34109
#> 
#> $bic
#> [1] 77.52978
#> 
#> $model
#> [1] "fivenets ~ edges + nodematch(\"female\")"
#> 
#> attr(,"class")
#> [1] "ergmito_summary"
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
