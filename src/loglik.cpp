#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

// Function to compute the normalizing constant
inline double kappa(
    const arma::colvec & params,
    const arma::rowvec & weights,
    const arma::mat    & statmat
) {
  
  return arma::as_scalar(weights * exp(statmat * params));
  
}

inline void exact_logliki(
    const arma::rowvec & x,
    const arma::colvec & params,
    const arma::rowvec & weights,
    const arma::mat    & statmat,
    arma::vec & ans,
    int i,
    bool as_prob = false
) {
  
  if (!as_prob) {
    ans.at(i) = arma::as_scalar(x * params) - 
      log(kappa(params, weights, statmat));
  } else {
    ans.at(i) = exp(arma::as_scalar(x * params))/ 
      kappa(params, weights, statmat);
  }
  
  return;
  
}

// [[Rcpp::export(name = "exact_loglik.")]]
arma::vec exact_loglik(
    const arma::mat & x,
    const arma::colvec & params,
    const std::vector< arma::rowvec > & weights,
    const std::vector< arma::mat > & statmat,
    bool as_prob = false
) {

  arma::vec ans(x.n_rows);
  int n = x.n_rows;
  
  for (int i = 0; i < n; ++i)
    exact_logliki(x.row(i), params, weights.at(i), statmat.at(i), ans, i, as_prob);
  
  return ans;
  
}
