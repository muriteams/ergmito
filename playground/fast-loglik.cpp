#include <RcppArmadillo.h>
// #include <omp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

inline void exact_logliki(
    const arma::colvec & stat0,
    const arma::rowvec & params,
    const arma::rowvec & weights,
    const arma::mat    & statsmat,
    arma::vec & ans,
    int i
) {
  
  ans.at(i) = arma::as_scalar(params * stat0) - 
    log(arma::as_scalar(weights * exp(statsmat * params.t())));
  
  return;
  
}

// [[Rcpp::export]]
arma::vec exact_loglik(
    arma::mat stat0,
    const arma::rowvec params,
    const arma::rowvec weights,
    const arma::mat statsmat,
    int ncores = 1
) {
  
  // // Setting number of threads
  // omp_set_num_threads(ncores);
  
  arma::vec ans(stat0.n_rows);
  int n = stat0.n_rows;

// #pragma omp parallel for shared(ans) firstprivate(params, weights, statsmat, n, stat0) \
//   default(none)
  for (int i = 0; i < n; ++i)
    exact_logliki(stat0.row(i).t(), params, weights, statsmat, ans, i);
  
  return ans;
  
}
