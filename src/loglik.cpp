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

// Calculates the likelihood for a given network individually.
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
  
  // Checking the sizes
  if (weights.size() != statmat.size())
    stop("The weights and statmat lists must have the same length.");
  
  if (weights.size() > 1u) {
    
    for (int i = 0; i < n; ++i)
      exact_logliki(x.row(i), params, weights.at(i), statmat.at(i), ans, i, as_prob);
    
  } else {
    
    ans = x * params - log(kappa(params, weights.at(0), statmat.at(0)));
    
  }
  
  return ans;
  
}
