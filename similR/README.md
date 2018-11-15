
<!-- README.md is generated from README.Rmd. Please edit that file -->

# similR: Similarity and Distance metrics for binary matrices

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/similR)](https://cran.r-project.org/package=similR)

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("USCCANA/social-smarts/r-pkg/similR")
```

## Example

An example from the manual

``` r
library(similR)

data(powerset03)
 
# We can compute it over the entire set
head(similarity(powerset03, statistic="s14"))
#>      i j        s14
#> [1,] 1 2  0.6324555
#> [2,] 1 3 -0.2000000
#> [3,] 1 4  0.6324555
#> [4,] 1 5  0.4472136
#> [5,] 1 6 -0.3162278
#> [6,] 1 7 -0.2000000

# Or over two pairs
head(similarity(powerset03[[1]], powerset03[[2]], powerset03[[3]], statistic="s14"))
#>      i j        s14
#> [1,] 1 2  0.6324555
#> [2,] 1 3 -0.2000000
#> [3,] 2 3  0.6324555

# We can compute multiple distances at the same time
ans <- similarity(powerset03, statistic=c("hamming", "dennis", "jaccard"))
head(ans)
#>      i j   hamming     dennis   jaccard
#> [1,] 1 2 0.1666667  1.6329932 0.0000000
#> [2,] 1 3 0.3333333 -0.5773503 0.1666667
#> [3,] 1 4 0.1666667  1.6329932 0.0000000
#> [4,] 1 5 0.3333333  1.0000000 0.0000000
#> [5,] 1 6 0.5000000 -0.8164966 0.2000000
#> [6,] 1 7 0.3333333 -0.5773503 0.1666667
```
