#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

inline void exact_logliki(
    const arma::colvec & x,
    const arma::rowvec & params,
    const arma::rowvec & weights,
    const arma::mat    & statmat,
    arma::vec & ans,
    int i
) {
  
  ans.at(i) = arma::as_scalar(params * x) - 
    log(arma::as_scalar(weights * exp(statmat * params.t())));
  
  return;
  
}

// [[Rcpp::export(name = "exact_loglik.")]]
arma::vec exact_loglik(
    const arma::mat & x,
    const arma::rowvec & params,
    const arma::rowvec & weights,
    const arma::mat & statmat
) {

  arma::vec ans(x.n_rows);
  int n = x.n_rows;

  for (int i = 0; i < n; ++i)
    exact_logliki(x.row(i).t(), params, weights, statmat, ans, i);
  
  return ans;
  
}
