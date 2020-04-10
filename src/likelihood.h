#include <RcppArmadillo.h>
#include "ergmito_types.h"

#ifndef ERGMITO_LIKELIHOOD_H
#define ERGMITO_LIKELIHOOD_H 1

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
class ergmito_ptr {
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
  
  /* This bool marks the data in case that the networks in the model have
   * the same set of sufficient statistics. This is important since in such
   * a case, we don't need to check for matching lenghts. 
   */
  bool same_stats = false;

public:
  unsigned int n, k;
  arma::mat                     target_stats;
  std::vector< arma::rowvec * > stats_weights;
  std::vector< arma::mat * >    stats_statmat;
  
  // Offset terms
  arma::colvec target_offset;
  std::vector< arma::colvec * > stats_offset;
  
  // Initializing
  ergmito_ptr(
    const NumericMatrix           & target_stats_,
    const ListOf< NumericVector > & stats_weights_,
    const ListOf< NumericMatrix > & stats_statmat_,
    const NumericVector           & target_offset_,
    const ListOf< NumericVector > & stats_offset_
  );
  
  void update_normalizing_constant(const arma::colvec & params);
  
  // Destructor function
  ~ergmito_ptr() {
    
    // Cleaning up
    for (auto i = stats_weights.begin(); i != stats_weights.end(); ++i)
      delete *i;
    
    for (auto i = stats_statmat.begin(); i != stats_statmat.end(); ++i)
      delete *i;
    
    for (auto i = stats_offset.begin(); i != stats_offset.end(); ++i)
      delete *i;
    
  };
  
  // The loglikes
  arma::vec exact_loglik(
      const arma::colvec & params,
      bool                 as_prob = false
  );
  
  arma::vec exact_gradient(
      const arma::colvec & params,
      bool                 as_prob = false
  );
  
  
};

// Constructor
inline ergmito_ptr::ergmito_ptr(
    const NumericMatrix & target_stats_,
    const ListOf< NumericVector > & stats_weights_,
    const ListOf< NumericMatrix > & stats_statmat_,
    const NumericVector           & target_offset_,
    const ListOf< NumericVector > & stats_offset_
) :
  n(target_stats_.nrow()),
  k(target_stats_.ncol()),
  target_stats((double *) &target_stats_[0], target_stats_.nrow(), target_stats_.ncol(), false, true),
  target_offset((double *) &target_offset_[0], target_offset_.size(), false, true)
{
  
  // Do the networks share the same vector of weights?
  if (stats_weights_.size() == 1u)
    same_stats = true;
  
  // Checking dimmentions (only check against target_stats if different than 1)
  if (stats_weights_.size() != stats_statmat_.size())
    stop("Incorrect sizes. stats_weights and stats_statmat should have the same size");
  if (!same_stats && (target_stats_.nrow() != stats_weights_.size()))
    stop("Incorrect sizes. target_stats and stats_statmat should have the same size");
  
  // Checking offset
  if (target_stats.n_rows != target_offset.size())
    stop("The length of target_offset should be equal to the number of rows in target_stats.");
  
  if (stats_offset_.size() != stats_statmat_.size())
    stop("The length of stats_offset should be the same as stats_statmat.");
  
  // Resizing needed outputs
  res_loglik.set_size(n);
  res_gradient.set_size(k, n);
  current_parameters.set_size(k);
  
  // The size of these objects is conditional on whether the networks share
  // statistics.
  normalizing_constant.set_size(same_stats ? 1u: n);
  stats_weights.reserve(same_stats ? 1u: n);
  stats_statmat.reserve(same_stats ? 1u: n);
  stats_offset.reserve(same_stats ? 1u: n);
  
  exp_statmat_params.resize(same_stats ? 1u: n);
  
  // Initializing 
  for (unsigned int i = 0u; i < n; ++i) {
    
    // Checking the dimension, is it the same as the sufficient statistics?
    if (stats_weights_[i].size() != stats_statmat_[i].nrow())
      stop(
        "length(stats_weights[[%i]]) != nrow(stats_statmat[[%i]]).",
        i + 1, i + 1
      );
    
    // Checking that the statistics have the same number of columns
    if (k != (unsigned int) stats_statmat_[i].ncol())
      stop("ncol(stats_statmat[[%i]]) != ncol(target_stats).", i + 1);
    
    /* Using advanced constructors from armadillo. In this case, I'm building
     * the objects s.t. it uses the raw memory pointer to the original data in
     * R. This is actually a bit dangerous but in our case can be very efficient
     */
    
    stats_weights.push_back(
      new arma::rowvec(
        (double *) &stats_weights_[i][0],
        stats_weights_[i].size(),
        false,
        true
    ));
    
    stats_statmat.push_back(
      new arma::mat(
        (double *) &stats_statmat_[i][0],
        stats_statmat_[i].nrow(),
        stats_statmat_[i].ncol(),
        false,
        true
    ));
    
    stats_offset.push_back(
      new arma::colvec(
        (double *) &stats_offset_[i][0],
        stats_offset_[i].size(),
        false,
        true
    ));
    
    // Preparing the exp statmat
    exp_statmat_params[i].set_size(stats_weights[i]->size());
    
    // This should only be done once if the same stats
    if (same_stats)
      break;
    
  }
  
  return;
  
}

inline void ergmito_ptr::update_normalizing_constant(
    const arma::colvec & params
  ) {
  
  /* If the current set of parameters equals the latest computed, then we 
   * can skip computing some components of this, in particular, the
   * normalizing constant. 
   */
  bool needs_update =
    this->first_iter ? true :
    arma::any(abs(params - this->current_parameters) > 1e-20);
  if (needs_update) {
    
    // Storing the current version of the parameters
    this->first_iter = false;
    std::copy(params.begin(), params.end(), this->current_parameters.begin());
    
    // Recalculating the normalizing constant and  exp(s(.) * theta)
    for (unsigned int i = 0; i < this->n; ++i) {
      
      this->exp_statmat_params[i] = exp(
        (*this->stats_statmat[i]) * params - AVOID_BIG_EXP +
          (*this->stats_offset[i])
      );
      
      this->normalizing_constant[i] = arma::as_scalar(
        (*this->stats_weights[i]) * this->exp_statmat_params[i]
      );
      
      // No need to recompute this if all share the same sufficient stats space.
      if (same_stats)
        break;
    }
    
  }
  
  return;
}

// Calculates the likelihood for a given network individually.
inline arma::vec ergmito_ptr::exact_loglik(
    const arma::colvec & params,
    bool as_prob
) {

  // Checking if we need to update the normalizing constant
  update_normalizing_constant(params);
  
  for (unsigned int i = 0; i < this->n; ++i) {
    
    unsigned int j = this->same_stats ? 0u : i;
    
    if (!as_prob) {
      
      this->res_loglik[i] =
        arma::as_scalar(this->target_stats.row(i) * params) -
        AVOID_BIG_EXP +
        this->target_offset[i] -
        log(normalizing_constant[j]);
      
    } else {
      
      this->res_loglik[i] =
        exp(
          arma::as_scalar(this->target_stats.row(i) * params) -
            AVOID_BIG_EXP + 
            this->target_offset[i]
      )/normalizing_constant[j];
      
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
inline arma::colvec ergmito_ptr::exact_gradient(
    const arma::colvec & params,
    bool as_prob
) {
  
  // Checking if we need to update the normalizing constant
  update_normalizing_constant(params);
  
  for (unsigned int i = 0u; i < n; ++i) {
    
    unsigned int j = this->same_stats ? 0u : i;
    
    res_gradient.col(i) = this->target_stats.row(i).t() - (
        this->stats_statmat[j]->t() * ( 
            this->stats_weights[j]->t() % this->exp_statmat_params[j]
        )) / this->normalizing_constant[j];
    
  }
  
  return sum(res_gradient, 1);
  
  
}

#endif

