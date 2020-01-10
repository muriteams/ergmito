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
  
  // These two objects contain outputs that will be used over and over again
  // during calls to the likelihood and gradient functions.
  arma::vec res_loglik;
  arma::mat res_gradient;
  
  /* Since the log-likelihood and the gradient function use both the normalizing
   * constant, we don't need to recompute it if the current set of parameters
   * is the same as the latest used.
   */
  arma::colvec current_parameters;
  arma::vec normalizing_constant;
  std::vector< arma::mat > exp_statmat_params;
  bool first_iter = true;
  
public:
  unsigned int n, k;
  arma::mat                   target_stats;
  std::vector< arma::rowvec > stats_weights;
  std::vector< arma::mat >    stats_statmat;
  
  // Initializing
  ergmito_model(
    arma::mat                   target_stats_,
    std::vector< arma::rowvec > stats_weights_,
    std::vector< arma::mat >    stats_statmat_
  );
  
  void update_normalizing_constant(const arma::colvec & params);
  
  // Destructor function
  ~ergmito_model() {};
  
  // The loglikes
  arma::vec exact_loglik(
      const arma::colvec & params,
      bool                 as_prob = false,
      int                  ncores  = 1
  );
  
  arma::vec exact_gradient(
      const arma::colvec & params,
      bool                 as_prob = false,
      int                  ncores  = 1
  );
  
  
};

// Constructor
inline ergmito_model::ergmito_model(
    arma::mat                   target_stats_,
    std::vector< arma::rowvec > stats_weights_,
    std::vector< arma::mat >    stats_statmat_
) : n(target_stats_.n_rows), k(target_stats_.n_cols)
{
  
  // Checking dimmentions
  if (stats_weights_.size() != stats_statmat_.size())
    stop("Incorrect sizes. stats_weights and stats_statmat should have the same size");
  if (target_stats_.n_rows != stats_weights_.size())
    stop("Incorrect sizes. target_stats and stats_statmat should have the same size");
  
  // Resizing needed outputs
  res_loglik.resize(n);
  stats_weights.resize(n);
  stats_statmat.resize(n);
  exp_statmat_params.resize(n);
  
  // Initializing 
  for (unsigned int i = 0; i < n; ++i) {
    
    // Checking the dimension, is it the same as the sufficient statistics?
    if (stats_weights_[i].size() != stats_statmat_[i].n_rows)
      stop(
        "length(stats_weights[[%i]]) != nrow(stats_statmat[[%i]]).",
        i + 1, i + 1
      );
    
    // Checking that the statistics have the same number of columns
    if (k != stats_statmat_[i].n_cols)
      stop("ncol(stats_statmat[[%i]]) != ncol(target_stats).", i + 1);
    
    /* Using advanced constructors from armadillo. In this case, I'm building
     * the objects s.t. it uses the raw memory pointer to the original data in
     * R. This is actually a bit dangerous but in our case can be very efficient
     */
    
    arma::rowvec tmpvec(
        stats_weights_[i].memptr(),
        stats_weights_[i].size(),
        false,
        true
    );
    
    arma::mat tmpmat(
        stats_statmat_[i].memptr(),
        stats_statmat_[i].n_rows,
        stats_statmat_[i].n_cols,
        false,
        true
    );
    
    // Moving the data to the desired location. This way we avoid duplicating
    // memory
    stats_weights[i] = std::move(tmpvec);
    stats_statmat[i] = std::move(tmpmat);
    exp_statmat_params[i].resize(stats_weights[i].size());
    
  }
  
  // Moving data
  arma::mat tmpmat(
      target_stats_.memptr(),
      target_stats_.n_rows,
      target_stats_.n_cols,
      false,
      true
  );
  target_stats.resize(tmpmat.n_rows, tmpmat.n_cols);
  target_stats = std::move(tmpmat);
  
  // Setting the sizes of commonly used containers
  res_gradient.resize(target_stats.n_cols, n);
  normalizing_constant.resize(n);
  current_parameters.resize(k);
  
  return;
  
};

inline void ergmito_model::update_normalizing_constant(const arma::colvec & params) {
  
  /* If the current set of parameters equals the latest computed, then we 
   * can skip computing some components of this, in particular, the
   * normalizing constant. 
   */
  if (this->first_iter | !arma::approx_equal(params, this->current_parameters, "absdiff", 1e-30)) {
    
    // Storing the current version of the parameters
    this->first_iter = false;
    std::copy(params.begin(), params.end(), this->current_parameters.begin());
    
    // Recalculating the normalizing constant and  exp(s(.) * theta)
    for (unsigned int i = 0; i < this->n; ++i) {
      
      this->exp_statmat_params[i] = exp(this->stats_statmat[i] * params - AVOID_BIG_EXP);
      
      this->normalizing_constant[i] = arma::as_scalar(
        this->stats_weights[i] * this->exp_statmat_params[i]
      );
    }
    
  }
  
  return;
}

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


  // Checking if we need to update the normalizing constant
  update_normalizing_constant(params);
  

#pragma omp parallel for shared(this->target_stats, this->stats_weights, this->stats_statmat, this->res) \
  default(shared) firstprivate(params, as_prob, n)
    for (unsigned int i = 0; i < this->n; ++i) {
      
      if (!as_prob) {
        
        this->res_loglik[i] =
          arma::as_scalar(this->target_stats.row(i) * params) -
          AVOID_BIG_EXP - log(normalizing_constant[i]);
        
      } else {
        
        this->res_loglik[i] =
          exp(arma::as_scalar(this->target_stats.row(i) * params) - AVOID_BIG_EXP)/ 
          normalizing_constant[i];
        
      }
      
    }
  
  return this->res_loglik;
  
}

//' Vectorized version of gradient function
//' 
//' @param x Matrix of statistic. `nnets * nstats`.
//' @param params Vector of coefficients.
//' @param weights A list of weights matrices (for `statmat`).
//' @param statmat A list of matrices with statistics for each row in `x`.
//' @noRd
inline arma::colvec ergmito_model::exact_gradient(
    const arma::colvec & params,
    bool as_prob,
    int ncores
) {
  
  // Setting the cores (not used right now)
#ifndef SKIP_OMP
  omp_set_num_threads(ncores);
#endif
  
  // Checking if we need to update the normalizing constant
  update_normalizing_constant(params);
  
#pragma omp parallel for shared(x, this->stats_weights, this->stats_statmat) \
    default(shared) firstprivate(params, n)
    for (unsigned int i = 0u; i < n; ++i) {
      
      // Speeding up a bit calculations (this is already done)
      // arma::colvec exp_stat_params = exp(
      //   this->stats_statmat.at(i) * params - AVOID_BIG_EXP
      // );
      
      res_gradient.col(i) = this->target_stats.row(i).t() - (
          this->stats_statmat[i].t() * (
              this->stats_weights[i].t() % this->exp_statmat_params[i]
          )) / this->normalizing_constant[i];
      
    }
  
  return sum(res_gradient, 1);
  
  
}

//' Creates new pointer
//' @export
// [[Rcpp::export]]
SEXP new_ergmito_model(
    const arma::mat & target_stats,
    const std::vector< arma::rowvec > & stats_weights,
    const std::vector< arma::mat > & stats_statmat) {
  
  Rcpp::XPtr< ergmito_model > ptr(
      new ergmito_model(target_stats, stats_weights, stats_statmat),
      true
  );
  
  ptr.attr("class") = "ergmito_model";
  
  return wrap(ptr);
  
}

//' Vectorized version of log-likelihood function
//' 
//' @param x Matrix of statistic. `nnets * nstats`.
//' @param params Vector of coefficients.
//' @param weights A list of weights matrices (for `statmat`).
//' @param statmat A list of matrices with statistics for each row in `x`.
//' @noRd
// [[Rcpp::export(name = "exact_loglik2.", rng = false)]]
arma::vec exact_loglik2(
    SEXP ptr,
    const arma::colvec & params,
    bool as_prob = false,
    int ncores = 1
) {
  
  Rcpp::XPtr< ergmito_model > p(ptr);
  return p->exact_loglik(params, as_prob, ncores);
  
}

//' Vectorized version of log-likelihood function
//' 
//' @param x Matrix of statistic. `nnets * nstats`.
//' @param params Vector of coefficients.
//' @param weights A list of weights matrices (for `statmat`).
//' @param statmat A list of matrices with statistics for each row in `x`.
//' @noRd
// [[Rcpp::export(name = "exact_gradient2.", rng = false)]]
arma::vec exact_gradient2(
    SEXP ptr,
    const arma::colvec & params,
    bool as_prob = false,
    int ncores = 1
) {
  
  Rcpp::XPtr< ergmito_model > p(ptr);
  return p->exact_gradient(params, as_prob, ncores);
  
}
