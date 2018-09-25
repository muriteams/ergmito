
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

## When `ergm` is not enough

``` r
library(lergm)
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
ans_lergm <- lergm(model)
ans_ergm  <- ergm(model)

# The lergm should have a larger value
ergm.exact(ans_lergm$coef, model)
#>           [,1]
#> [1,] -6.557684
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
#> Iterations:  100 out of 20 
#> 
#> Monte Carlo MLE Results:
#>         Estimate Std. Error MCMC % z value Pr(>|z|)
#> edges     -0.317      1.917     29  -0.165    0.869
#> mutual     2.330      2.965     29   0.786    0.432
#> balance   -7.412     53.755     29  -0.138    0.890
#> 
#>      Null Deviance: 16.64  on  0  degrees of freedom
#>  Residual Deviance: 13.12  on -3  degrees of freedom
#>  
#> AIC: 31.12    BIC: 35.48    (Smaller is better.)
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
#> edges   -0.01198    1.20677      0   -0.01    0.992    
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
ans_lergm <- lergm(model)
ans_ergm  <- ergm(model, control = control.ergm(
  MCMC.effectiveSize = 4000,
  seed = 444)
  )

# The lergm should have a larger value
ergm.exact(ans_lergm$coef, model) > ergm.exact(ans_ergm$coef, model)
#>      [,1]
#> [1,] TRUE

summary(ans_lergm)
#> 
#> Little ERGM estimates
#> 
#> ==========================
#> Summary of model fit
#> ==========================
#> 
#> Formula:   net ~ edges + mutual
#> 
#> Iterations:  9 out of 20 
#> 
#> Monte Carlo MLE Results:
#>        Estimate Std. Error MCMC % z value Pr(>|z|)
#> edges    -1.099      1.291     29  -0.851    0.395
#> mutual    1.098      2.582     29   0.425    0.671
#> 
#>      Null Deviance: 16.64  on  0  degrees of freedom
#>  Residual Deviance: 14.91  on -2  degrees of freedom
#>  
#> AIC: 34.91    BIC: 39.76    (Smaller is better.)
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
