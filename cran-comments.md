## Test environments

* Ubuntu 18.04 R 3.6.2 (local).
* Ubuntu 16.04.6 LTS, R 3.6.2 (on travis-ci)
* Ubuntu 16.04.6 LTS, R Under development (unstable) (2020-01-13 r77662)
* macOS Mojave 10.14.4, R 3.6.2 (on travis-ci)
* Windows release R 3.6.2 (64 and 32 on AppVeyor)
* Win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note (different depending on the OS)

* This is a new release.

* 1 note on windows: Possible mispelling unknown word ERGMs: It is the right
  spelling. In the literature, ERGM is used for single models while ERGMs for
  multiple.

* 1 note on Ubuntu: The note is regarding the size of the package. Since I'm
  using RcppArmadillo I don't have much control on the size of it. Yet, I'm
  constantly working on reducing it's size.
  
