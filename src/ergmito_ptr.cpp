#include <RcppArmadillo.h>
#include "ergmito_types.h"
#include "likelihood.h"

using namespace Rcpp;

// [[Rcpp::export(rng = false, name = "new_ergmito_ptr_cpp")]]
SEXP new_ergmito_ptr(
    const NumericMatrix           & target_stats,
    const ListOf< NumericVector > & stats_weights,
    const ListOf< NumericMatrix > & stats_statmat,
    const NumericVector           & target_offset,
    const ListOf< NumericVector > & stats_offset
    ) {
  
  Rcpp::XPtr< ergmito_ptr > ptr(
    new ergmito_ptr(
      target_stats, stats_weights, stats_statmat,
      target_offset, stats_offset
    ),
    true
  );
  
  ptr.attr("class") = "ergmito_ptr";
  
  return wrap(ptr);
  
}

//' Vectorized version of log-likelihood function
//' 
//' @param x Matrix of statistic. `nnets * nstats`.
//' @param params Vector of coefficients.
//' @param weights A list of weights matrices (for `statmat`).
//' @param statmat A list of matrices with statistics for each row in `x`.
//' @noRd
// [[Rcpp::export(name = "exact_loglik_cpp", rng = false)]]
arma::vec exact_loglik(
    SEXP ptr,
    const arma::colvec & params,
    bool as_prob = false
) {
  
  Rcpp::XPtr< ergmito_ptr > p(ptr);
  return p->exact_loglik(params, as_prob);
  
}

//' Vectorized version of log-likelihood function
//' 
//' @param x Matrix of statistic. `nnets * nstats`.
//' @param params Vector of coefficients.
//' @param weights A list of weights matrices (for `statmat`).
//' @param statmat A list of matrices with statistics for each row in `x`.
//' @noRd
// [[Rcpp::export(name = "exact_gradient_cpp", rng = false)]]
arma::vec exact_gradient(
    SEXP ptr,
    const arma::colvec & params,
    bool as_prob = false
) {
  
  Rcpp::XPtr< ergmito_ptr > p(ptr);
  return p->exact_gradient(params, as_prob);
  
}

// [[Rcpp::export(rng = false)]]
List get_boundaries(SEXP ptr) {
  Rcpp::XPtr< ergmito_ptr > p(ptr);
  std::vector< arma:: mat > ans = p->get_boundaries();
  std::vector< std::vector< arma::uvec > > idx = p->get_dx_matches_target();
  return List::create(
    _["lower"] = wrap(ans[0u]),
    _["upper"] = wrap(ans[1u]),
    _["obs"]   = wrap(ans[2u]),
    _["match"] = wrap(idx)
  );
}


// Calculates the gradient for a given network individually.
inline arma::mat exact_hessiani(
    const arma::colvec & params,
    const arma::rowvec & stats_weights,
    const arma::mat    & stats_statmat,
    const arma::colvec & stats_offset
) {
  
  // Speeding up a bit calculations (this is already done)
  arma::colvec Z = exp(stats_statmat * params + stats_offset);
  double weighted_exp = arma::as_scalar(stats_weights * Z);
  double weighted_exp2 = pow(weighted_exp, 2.0);
  arma::rowvec WZ = stats_weights % Z.t();
  // Z.print("\nZ");
  // stats_weights.print("\nstats_weights");
  
  unsigned int K = params.size();
  arma::mat H(K, K);
  arma::rowvec WZS = WZ * stats_statmat;
  
  for (unsigned int k0 = 0u; k0 < K; ++k0) {
    for (unsigned int k1 = 0u; k1 <= k0; ++k1) {
      H(k0, k1) = - (arma::as_scalar(
        WZ * (stats_statmat.col(k0) % stats_statmat.col(k1)) * weighted_exp
      ) - WZS[k0] * WZS[k1]) / weighted_exp2;
      
      if (k0 != k1)
        H(k1, k0) = H(k0, k1);
    }
  }
  
  return H;
  
}


//' Vectorized version of gradient function
//' 
//' @param params Vector of coefficients.
//' @param weights A list of weights matrices (for `statmat`).
//' @param statmat A list of matrices with statistics for each row in `x`.
//' @noRd
// [[Rcpp::export(name = "exact_hessian_cpp", rng = false)]]
arma::mat exact_hessian(
    const arma::colvec                & params,
    const std::vector< arma::rowvec > & stats_weights,
    const std::vector< arma::mat >    & stats_statmat,
    const std::vector< arma::colvec > & stats_offset
) {
  
  // Checking the sizes
  if (stats_weights.size() != stats_statmat.size())
    stop("The weights and statmat lists must have the same length.");
  
  if (stats_weights.size() > 1u) {
    
    unsigned int n = stats_weights.size();
    unsigned int K = params.size();
    std::vector< arma::mat > ans(n);
    for (unsigned int i = 0u; i < n; ++i)
      ans[i].set_size(K, K);
    
    for (unsigned int i = 0u; i < n; ++i)
      ans[i] = exact_hessiani(
        params, stats_weights.at(i), stats_statmat.at(i),
        stats_offset.at(i)
      );
    
    arma::mat H(K, K);
    H = ans[0u];
    for (unsigned int i = 1u; i < n; ++i)
      H += ans[i];
    
    return H;
    
  } else {
    // In the case that all networks are from the same family, then this becomes
    // a trivial operation.
    return exact_hessiani(
      params, stats_weights.at(0u), stats_statmat.at(0u),
      stats_offset.at(0u)
    );
    
  }
  
}

