#include <RcppArmadillo.h>
#include "ergmito_types.h"

#ifndef SKIP_OMP
#include <omp.h>
#endif


using namespace Rcpp;

/**
 * This class contains the needed objects to calculate likelihood. The idea is
 * that we don't need to pass matrices and vectors over and over again, so we
 * can simplify computations by just updating the parameter.
 * 
 * This was motivated by the fact that, during the optimization process, too
 * much memory is allocated and de-allocated. This is an effort to solve this
 * problem.
 */
class ergmito_model {
private:
  arma::vec res_loglik;
  arma::mat res_gradient;
  
public:
  const unsigned int n;
  const arma::mat                   & target_stats;
  const std::vector< arma::rowvec > & stats_weights;
  const std::vector< arma::mat >    & stats_statmat;
  
  // Initializing
  ergmito_model(
    const arma::mat                   & target_stats_,
    const std::vector< arma::rowvec > & stats_weights_,
    const std::vector< arma::mat >    & stats_statmat_
  ) :
    target_stats(target_stats_),
    stats_weights(stats_weights_),
    stats_statmat(stats_statmat_),
    n(target_stats_.n_rows)
    {
    
    // Checking dimmentions
    if (stats_weights.size() != stats_statmat.size())
      stop("Incorrect sizes. stats_weights and stats_statmat should have the same size");
    if (target_stats.n_rows != stats_weights.size())
      stop("Incorrect sizes. target_stats and stats_statmat should have the same size");
    
    res_loglik.resize(n);
    res_gradient.resize(target_stats.n_cols, n);
    
    return;
    
  };
  
  // Destructor function
  ~ergmito_model() {};
  
  // The loglikes
  arma::vec exact_loglik(const arma::colvec & params, bool as_prob = false, int ncores = 1);
  arma::vec exact_gradient(const arma::colvec & params, bool as_prob = false, int ncores = 1);
  
};

// Calculates the likelihood for a given network individually.
inline arma::vec ergmito_model::exact_loglik(
    const arma::colvec & params,
    bool as_prob,
    int ncores
) {
  
#ifndef SKIP_OMP
  // Setting the cores
  omp_set_num_threads(ncores);
#endif

#pragma omp parallel for shared(this->target_stats, this->stats_weights, this->stats_statmat, this->res) \
  default(shared) firstprivate(params, as_prob, n)
    for (int i = 0; i < n; ++i) {
      
      if (!as_prob) {
        this->res_loglik[i] = arma::as_scalar(this->target_stats.row(i) * params) - AVOID_BIG_EXP -
          log(kappa(params, this->stats_weights.at(i), this->stats_statmat.at(i)));
      } else {
        this->res_loglik[i] = exp(arma::as_scalar(this->target_stats.row(i) * params) - AVOID_BIG_EXP)/ 
          kappa(params, this->stats_weights.at(i), this->stats_statmat.at(i));
      }
      
    }
  
  return this->res_loglik;
  
}

//' Vectorized version of log-likelihood function
//' 
//' @param x Matrix of statistic. `nnets * nstats`.
//' @param params Vector of coefficients.
//' @param weights A list of weights matrices (for `statmat`).
//' @param statmat A list of matrices with statistics for each row in `x`.
//' @noRd
// [[Rcpp::export(name = "exact_loglik2.", rng = false)]]
arma::vec exact_loglik(
    SEXP ptr,
    const arma::colvec & params,
    bool as_prob = false,
    int ncores = 1
) {
  
  Rcpp::XPtr< ergmito_model > p(ptr);
  return p->exact_loglik(params, ncores, as_prob);
  
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
// [[Rcpp::export(name = "exact_gradient2.", rng = false)]]
arma::colvec ergmito_model::exact_gradient(
    const arma::colvec & params,
    bool as_prob,
    int ncores
) {
  
  // Setting the cores (not used right now)
#ifndef SKIP_OMP
  omp_set_num_threads(ncores);
#endif
  
  // this->res_gradient.fill(0.0);
  
#pragma omp parallel for shared(x, stats_weights, stats_statmat, ans) \
    default(shared) firstprivate(params, n)
    for (unsigned int i = 0u; i < n; ++i)
      res_gradient.col(i) = exact_gradienti(
        x.row(i), params, stats_weights.at(i), stats_statmat.at(i)
      );
  
  return sum(ans, 1);
  
  
}

//' Creates new pointer
//' @export
// [[Rcpp::export]]
SEXP new_ergmito_model(
    const arma::mat & target_stats,
    const std::vector< arma::rowvec > & stats_weights,
    const std::vector< arma::mat > & stats_statmat) {
  
  return wrap(
    Rcpp::XPtr< ergmito_model >(
        new ergmito_model(target_stats, stats_weights, stats_statmat),
        true
    )
  );
  
}


