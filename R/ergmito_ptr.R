#' Creates a new `ergmito_ptr`
#' 
#' After calculating the support of the sufficient statistics, the second
#' most computationally expensive task is computing log-likelihoods, 
#' Gradients, and Hessian matrices of ERGMs. This function creates a pointer to an
#' underlying class that is optimized to improve memory allocation and 
#' save computation time when possible.
#' 
#' @details This function is for internal used only. Non-advanced users
#' are not encouraged to use it. See [ergmito_formulae] and [exact_loglik]
#' for user friendly wrappers of this function.
#' @section Recycling computations:
#' 
#' Some components of the likelihood, its gradient, and hessian can be 
#' pre-computed and recycled when needed. For example, it is usually the
#' case that in optimization gradients are computed using a current state
#' of the model's parameter, which implies that the normalizing constant
#' and some other matrix products will be the same between the log-likelihood
#' and the gradient. Because of this, the underlying class `ergmito_ptr`
#' will only re-calculate these shared components if the parameter used
#' changes as well. This saves a significant amount of computation time.
#' 
#' @section Scope of the class methods:
#' 
#' To save space, the class creates pointers to the matrices of sufficient
#' statistics that the model uses. This means that once these objects are
#' deleted the log-likelihood and the gradient functions become invalid
#' from the computational point of view. 
#' 
#' @param target_stats,stats_weights,stats_statmat,target_offset,stats_offset see [exact_loglik].
#' 
new_ergmito_ptr <- function(
  target_stats,
  stats_weights,
  stats_statmat,
  target_offset,
  stats_offset
) {
  
  stopncall <- function(...) stop(..., call. = FALSE)

  # Checking target stats
  if (!inherits(target_stats, "matrix"))
    stopncall("-target_stats- must be an object of class matrix.")
  
  if (!is.numeric(target_stats))
    stopncall("-target_stats- should be numeric.")
  
  # Weights should be numeric vectors
  if (any(!sapply(stats_weights, is.numeric))) 
    stopncall("All the vectors in -stats_weights- should be numeric.")
  
  if (any(!sapply(stats_weights, is.vector))) 
    stopncall("All the elements in -stats_weights- should be vectors.")
  
  # Statmat should be numeric matrices
  if (any(!sapply(stats_statmat, inherits, what = "matrix"))) 
    stopncall("All the elements in -stats_statmat- should be matrix.")
  
  if (any(!sapply(stats_statmat, is.numeric))) 
    stopncall("All the elements in -stats_statmat- should be numeric.")
  
  # Stats offset are numeric vectors
  if (any(!sapply(stats_offset, is.numeric))) 
    stopncall("All the vectors in -stats_offset- should be numeric.")
  
  if (any(!sapply(stats_offset, is.vector))) 
    stopncall("All the elements in -stats_offset- should be vectors.")
  
  new_ergmito_ptr.(
    target_stats  = target_stats,
    stats_weights = stats_weights,
    stats_statmat = stats_statmat,
    target_offset = target_offset,
    stats_offset  = stats_offset
  )
  
}