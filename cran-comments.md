## Test environments

* Ubuntu 18.04 LTS, R 3.6.2 (locally).
* Ubuntu 18.04 LTS, R 3.6.2 with valgrind (locally).
* Ubuntu 16.04.6 LTS, R 3.6.2 (on travis-ci)
* Ubuntu 16.04.6 LTS, R Under development (unstable) (2020-01-13 r77662) (on travis-ci)
* macOS Mojave 10.14.4, R 3.6.2 (on travis-ci)
* Windows release R 3.6.2 (64 and 32 on AppVeyor)
* Windows, R version 3.6.2 (2019-12-12) (win-builder)
* Windows, R Under development (unstable) (2020-01-07 r77633) (win-builder)

## R CMD check results

0 errors | 0 warnings | 1 note (different depending on the OS)

* This is a new release.

* 1 note on windows: Possible mispelling unknown word ERGMs: It is the right
  spelling. In the literature, ERGM is used for single models while ERGMs for
  multiple.

* 1 note on Ubuntu: The note is regarding the size of the package. Since I'm
  using RcppArmadillo I don't have much control on the size of it. Yet, I'm
  constantly working on reducing it's size.
  
