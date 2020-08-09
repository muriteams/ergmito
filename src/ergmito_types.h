#include<RcppArmadillo.h>

using namespace Rcpp;

#ifndef ERGMITO_TYPES_H
#define ERGMITO_TYPES_H 1

typedef std::vector< int > vecint;
typedef std::vector< std::vector< int > > vecvecint;
#define AVOID_BIG_EXP 500

// Function to compute the normalizing constant
inline double kappa(
    const arma::colvec & params,
    const arma::rowvec & weights,
    const arma::mat    & statmat
) {
  
  return arma::as_scalar(weights * exp(statmat * params - AVOID_BIG_EXP));
  
}

#endif
