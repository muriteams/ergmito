#include <RcppArmadillo.h>

#ifndef SKIP_OMP
  #include <omp.h>
#endif

#include "ergmito_types.h"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

// Function to compute the normalizing constant
inline double kappa(
    const arma::colvec & params,
    const arma::rowvec & weights,
    const arma::mat    & statmat
) {
  
  return arma::as_scalar(weights * exp(statmat * params - AVOID_BIG_EXP));
  
}

// Calculates the likelihood for a given network individually.
inline void exact_logliki(
    const arma::rowvec & x,
    const arma::colvec & params,
    const arma::rowvec & stats_weights,
    const arma::mat    & stats_statmat,
    arma::vec & ans,
    int i,
    bool as_prob = false
) {
  
  if (!as_prob) {
    ans.at(i) = arma::as_scalar(x * params) - AVOID_BIG_EXP -
      log(kappa(params, stats_weights, stats_statmat));
  } else {
    ans.at(i) = exp(arma::as_scalar(x * params) - AVOID_BIG_EXP)/ 
      kappa(params, stats_weights, stats_statmat);
  }
  
  return;
  
}

//' Vectorized version of loglikelihood function
//' 
//' @param x Matrix of statistic. `nnets * nstats`.
//' @param params Vector of coefficients.
//' @param weights A list of weights matrices (for `statmat`).
//' @param statmat A list of matrices with statistics for each row in `x`.
//' @noRd
// [[Rcpp::export(name = "exact_loglik.", rng = false)]]
arma::vec exact_loglik(
    const arma::mat & x,
    const arma::colvec params,
    const std::vector< arma::rowvec > & stats_weights,
    const std::vector< arma::mat > & stats_statmat,
    bool as_prob = false,
    int ncores = 1
) {

  int n = x.n_rows;
  
#ifndef SKIP_OMP
  // Setting the cores
  omp_set_num_threads(ncores);
#endif
  
  // Checking the sizes
  if (stats_weights.size() != stats_statmat.size())
    stop("The weights and statmat lists must have the same length.");
  
  if (stats_weights.size() > 1u) {
    
    arma::vec ans(x.n_rows);
    
#pragma omp parallel for shared(x, stats_weights, stats_statmat, ans) \
    default(shared) firstprivate(params, as_prob, n)
    for (int i = 0; i < n; ++i)
      exact_logliki(
        x.row(i), params, stats_weights.at(i), stats_statmat.at(i), ans, i,
        as_prob
      );

    return ans;
  
  } else
    return x * params - AVOID_BIG_EXP -
      log(kappa(params, stats_weights.at(0), stats_statmat.at(0)));
  
}


// Calculates the gradient for a given network individually.
inline arma::colvec exact_gradienti(
    const arma::rowvec & x,
    const arma::colvec & params,
    const arma::rowvec & stats_weights,
    const arma::mat    & stats_statmat
) {

  // Speeding up a bit calculations (this is already done)
  arma::colvec exp_stat_params = exp(stats_statmat * params - AVOID_BIG_EXP);
  
  return x.t() - (
      stats_statmat.t() * (
          stats_weights.t() % exp_stat_params
      ))/arma::as_scalar(stats_weights * exp_stat_params);

}

//' Vectorized version of gradient function
//' 
//' @param x Matrix of statistic. `nnets * nstats`.
//' @param params Vector of coefficients.
//' @param weights A list of weights matrices (for `statmat`).
//' @param statmat A list of matrices with statistics for each row in `x`.
//' @noRd
// [[Rcpp::export(name = "exact_gradient.", rng = false)]]
arma::colvec exact_gradient(
    const arma::mat & x,
    const arma::colvec params,
    const std::vector< arma::rowvec > & stats_weights,
    const std::vector< arma::mat > & stats_statmat,
    int ncores
) {

  // Checking the sizes
  if (stats_weights.size() != stats_statmat.size())
    stop("The weights and statmat lists must have the same length.");

  // Setting the cores (not used right now)
#ifndef SKIP_OMP
  omp_set_num_threads(ncores);
#endif
  
  if (stats_weights.size() > 1u) {
    
    int n = x.n_rows;
    arma::mat ans(x.n_cols, n);
    ans.fill(0.0);
    
#pragma omp parallel for shared(x, stats_weights, stats_statmat, ans) \
    default(shared) firstprivate(params, n)
    for (int i = 0; i < n; ++i)
      ans.col(i) = exact_gradienti(
        x.row(i), params, stats_weights.at(i), stats_statmat.at(i)
      );
    
    return sum(ans, 1);

  } else {
    // In the case that all networks are from the same family, then this becomes
    // a trivial operation.
    return exact_gradienti(
      x.row(0), params, stats_weights.at(0), stats_statmat.at(0)
    );

  }

}

// Calculates the gradient for a given network individually.
inline arma::mat exact_hessiani(
    const arma::rowvec & x,
    const arma::colvec & params,
    const arma::rowvec & stats_weights,
    const arma::mat    & stats_statmat
) {
  
  // Speeding up a bit calculations (this is already done)
  arma::colvec Z = exp(stats_statmat * params);
  // Z.print("Z");
  
  unsigned int K = params.size();
  unsigned int n = stats_weights.size();
  arma::mat H(K, K);
  H.fill(0.0);
  
    // Calculating another product
  // Computing the hessian
  double weighted_exp = arma::as_scalar(stats_weights * Z);
  
  Rprintf("weighted_exp: %.4f\n", weighted_exp);
  
  for (unsigned int k0 = 0u; k0 < K; ++k0) {
    
    for (unsigned int k1 = 0u; k1 <= k0; ++k1) {
      H(k0, k1) = -
        (
            arma::as_scalar(stats_weights * (Z % stats_statmat.col(k0) % stats_statmat.col(k1))) *
              weighted_exp -
              arma::as_scalar(stats_weights * (Z % stats_statmat.col(k0))) * 
              arma::as_scalar(stats_weights * (Z % stats_statmat.col(k1)))
      ) / pow(weighted_exp, 2.0);
      
      if (k0 != k1)
        H(k1, k0) = H(k0, k1);
    }
  }
  
  return H;
  
}

//' Vectorized version of gradient function
//' 
//' @param x Matrix of statistic. `nnets * nstats`.
//' @param params Vector of coefficients.
//' @param weights A list of weights matrices (for `statmat`).
//' @param statmat A list of matrices with statistics for each row in `x`.
//' @noRd
// [[Rcpp::export(name = "exact_hessian.", rng = false)]]
arma::mat exact_hessian(
    const arma::mat & x,
    const arma::colvec params,
    const std::vector< arma::rowvec > & stats_weights,
    const std::vector< arma::mat > & stats_statmat,
    int ncores
) {
  
  // Checking the sizes
  if (stats_weights.size() != stats_statmat.size())
    stop("The weights and statmat lists must have the same length.");
  
  // Setting the cores (not used right now)
#ifndef SKIP_OMP
  omp_set_num_threads(ncores);
#endif
  
  if (stats_weights.size() > 1u) {
    
    int n = x.n_rows;
    unsigned int K = params.size();
    std::vector< arma::mat > ans(n);
    for (unsigned int i = 0u; i < n; ++i)
      ans[i].set_size(K, K);
    
#pragma omp parallel for shared(x, stats_weights, stats_statmat, ans) \
    default(shared) firstprivate(params, n)
      for (int i = 0; i < n; ++i)
        ans[i] = exact_hessiani(
          x.row(i), params, stats_weights.at(i), stats_statmat.at(i)
        );
    
    arma::mat H(K, K);
    for (unsigned int i = 0u; i < n; ++i)
      H += ans[i];
    
    return H;
    
  } else {
    // In the case that all networks are from the same family, then this becomes
    // a trivial operation.
    return exact_hessiani(
      x.row(0), params, stats_weights.at(0), stats_statmat.at(0)
    );
    
  }
  
}

// // [[Rcpp::export]]
// arma::vec exact_gradient(
//     const arma::mat & x,
//     const arma::colvec & params,
//     const std::vector< arma::rowvec > & weights,
//     const std::vector< arma::mat > & statmat,
// ) {
// 
//   arma::vec ans(params.size());
// 
//   ans += x - (weights * exp(statmat * params) % statmat)/kappa(params, weights, statmat) ;
//   // arma::as_scalar(weights * exp(statmat * params))
//   return ans;
// 
// }