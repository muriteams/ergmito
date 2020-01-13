## Test environments

* Ubuntu 18.04 R 3.6.2 (local).
* ubuntu 18.04, R 3.6.2 (on travis-ci)
* macOS Mojave 10.14.4, R 3.6.2 (on travis-ci)
* Windows release R 3.6.2 (64 and 32 on AppVeyor)
* Win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

* The note is regarding the size of the package. Since I'm using RcppArmadillo
  I don't have much control on the size of it. Yet, I'm constantly working on
  reducing it's size.
  
