
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ergmito: Estimation of Little ‘ERGMs’ using exact likelihood <img src="man/figures/logo.png" align="right" width="180px"/>

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

## Citation

Vega Yon G (2018). *ergmito: Exponential Random Graph Models for Small
Networks*. R package version 0.1-0, \<URL:
<https://github.com/muriteams/ergmito>\>.

Vega Yon G, de la Haye K (2019). “Exponential Random Graph models for
Little Networks.” arXiv:1904.10406.

To see these entries in BibTeX format, use ‘print(<citation>,
bibtex=TRUE)’, ‘toBibtex(.)’, or set ‘options(citation.bibtex.max=999)’.

## Installation

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("muriteams/ergmito")
```

This requires compilation. Windows users can download the lates compiled
version from appveyor
[here](https://ci.appveyor.com/project/gvegayon/ergmito/build/artifacts).
The file to download is the one named `ergmito_[version number].zip`.
Once donwloaded, you can install typing the
following:

``` r
install.packages("[path to the zipfile]/ergmito_[version number].zip", repos = FALSE)
```

In the case of Mac users, and in particular, those with Mojave version,
they may need to install the following
<https://github.com/fxcoudert/gfortran-for-macOS/releases>

## Example

In the following example we simulate a small network with 4 vertices and
estimate the model parameters using `ergm` and `ergmito`. We start by
generating the graph

``` r
# Generating a small graph
library(ergmito)
library(ergm)
library(sna)

set.seed(12123)
n   <- 4
net <- rbernoulli(n, p = .3)
gplot(net)
```

<img src="man/figures/README-net2-1.png" width="80%" />

To estimate the model

``` r
model <- net ~ edges + istar(2)

# ERGMito (estimation via MLE)
ans_ergmito <- ergmito(model)

# ERGM (estimation via MC-MLE)
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
#> 
#> formula:  net ~ edges + istar(2) 
#> 
#>          Estimate Std. Error    z value  Pr(>|z|)
#> edges  -1.3774881   1.002596 -1.3739208 0.1694663
#> istar2  0.5065322   1.327007  0.3817103 0.7026763
summary(ans_ergm)
#> 
#> ==========================
#> Summary of model fit
#> ==========================
#> 
#> Formula:   net ~ edges + istar(2)
#> 
#> Iterations:  2 out of 20 
#> 
#> Monte Carlo MLE Results:
#>        Estimate Std. Error MCMC % z value Pr(>|z|)
#> edges   -1.4110     0.9995      0  -1.412    0.158
#> istar2   0.5607     1.3719      0   0.409    0.683
#> 
#>      Null Deviance: 16.64  on 12  degrees of freedom
#>  Residual Deviance: 13.35  on 10  degrees of freedom
#>  
#> AIC: 17.35    BIC: 18.32    (Smaller is better.)
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
#> 
#> ERGMito estimates
#> 
#> formula:  fivenets ~ edges + nodematch("female") 
#> 
#>                   Estimate Std. Error   z value    Pr(>|z|)
#> edges            -1.704748  0.5435573 -3.136280 0.001711055
#> nodematch.female  1.586965  0.6430475  2.467882 0.013591530
```

We can also compute GOF

``` r
fivenets_gof <- gof_ergmito(model1)
fivenets_gof
#> 
#> Goodness-of-fit for edges 
#> 
#>       obs min mean max lower upper lower prob. upper prob.
#> net 1   2   0  3.7  12     0     6      0.0081        0.96
#> net 2   7   0  3.7  12     0     6      0.0081        0.96
#> net 3   4   0  3.1  12     0     6      0.0206        0.99
#> net 4   5   0  5.6  12     2     8      0.0309        0.95
#> net 5   2   0  3.7  12     0     6      0.0081        0.96
#> 
#> 
#> Goodness-of-fit for nodematch.female 
#> 
#>       obs min mean max lower upper lower prob. upper prob.
#> net 1   2   0  2.8   6     0     5       0.022        0.99
#> net 2   5   0  2.8   6     0     5       0.022        0.99
#> net 3   3   0  1.9   4     0     3       0.079        0.95
#> net 4   5   0  5.6  12     2     8       0.031        0.95
#> net 5   1   0  2.8   6     0     5       0.022        0.99
#> 
#> Note: Exact confidence intervals where used. This implies that the requestes CI may differ from the one used (see ?gof_ergmito).
plot(fivenets_gof)
```

<img src="man/figures/README-fivenets-gof-1.png" width="80%" />

# Contributing

Please note that the ‘ergmito’ project is released with a [Contributor
Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.
