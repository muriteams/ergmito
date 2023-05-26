#' Vectorized calculation of ERGM exact log-likelihood
#' 
#' This function can be compared to [ergm::ergm.exact] with the statistics not
#' centered at `x`, the vector of observed statistics.
#' 
#' @param x Matrix. Observed statistics
#' @param params Numeric vector. Parameter values of the model.
#' @template stats
#' @param target_offset Numeric vector of length `nrow(target_stats)`.
#' @param stats_offset List of numeric vectors, each of length equal to the
#' lengths of vectors in `stats_weights` (see details).
#' @param ... Arguments passed to the default methods.
#' @param as_prob Logical scalar. When `TRUE`, the function returns probabilities
#' instead of log-likelihoods.
#' 
#' @section Sufficient statistics:
#' 
#' One of the most important components of `ergmito` is calculating the full
#' support of the model's sufficient statistics. Right now, the package uses
#' the function [ergm::ergm.allstats] which returns a list of two objects:
#'
#' - `weights`: An integer vector of counts.
#' - `statmat`: A numeric matrix with the rows as unique vectors of sufficient
#'   statistics.
#' 
#' Since `ergmito` can vectorize operations, in order to specify weights and
#' statistics matrices for each network in the model, the user needs to pass
#' two lists `stats_weights` and `stats_statmat`. While both lists have to
#' have the same length (since its elements are matched), this needs not to
#' be the case with the networks, as the user can specify a single set of
#' weights and statistics that will be recycled (smartly).
#' 
#' In the case of offset terms, these can be passed directly via `target_offset`
#' and `stats_offset`. The first is a numeric vector of length equal to the
#' number of networks in the model. The later is a list of vectors that is
#' matched against `stats_weights`, so each of it's elements must be a
#' numeric vector of the same length that in the list of weights. By default
#' the offset terms are set to equal zero.
#' 
#' 
#' @examples 
#' data(fivenets)
#' ans <- ergmito(fivenets ~ edges + nodematch("female"))
#' 
#' # This computes the likelihood for all the networks independently
#' with(ans$formulae, {
#'   exact_loglik(
#'     x      = target_stats,
#'     params = coef(ans),
#'     stats_weights = stats_weights,
#'     stats_statmat = stats_statmat
#'   )
#' })
#' 
#' # This should be close to zero
#' with(ans$formulae, {
#'   exact_gradient(
#'     x      = target_stats,
#'     params = coef(ans),
#'     stats_weights = stats_weights,
#'     stats_statmat = stats_statmat
#'   )
#' })
#' 
#' # Finally, the hessian
#' with(ans$formulae, {
#'   exact_hessian(
#'     params = coef(ans),
#'     stats_weights = stats_weights,
#'     stats_statmat = stats_statmat
#'   )
#' })
#' 
#' @export
exact_loglik <- function(
  x,
  params,
  ...,
  as_prob = FALSE
  )
  UseMethod("exact_loglik")

#' @export
# @rdname exact_loglik
exact_loglik.ergmito_ptr <- function(x, params, ..., as_prob) {
  exact_loglik_cpp(x, params = params, as_prob = as_prob)
}

#' @export
#' @rdname exact_loglik
exact_loglik.default <- function(
  x,
  params,
  stats_weights,
  stats_statmat,
  target_offset = double(nrow(x)),
  stats_offset  = lapply(stats_weights, function(i) double(length(i))),
  ...,
  as_prob = FALSE
) {
  
  # Need to calculate it using chunks of size 200, otherwise it doesn't work(?)
  n      <- nrow(x)
  chunks <- make_chunks(n, 4e5)

  if (!is.list(stats_weights))
    stats_weights <- list(stats_weights)
  
  if (!is.list(stats_statmat))
    stats_statmat <- list(stats_statmat)
  
  if (!is.list(stats_offset))
    stats_offset <- list(stats_offset)
  
  # Computing in chunks
  ans <- vector("double", n)
  if (length(stats_weights) > 1L) {
    
    for (s in seq_along(chunks$from)) {
      
      i <- chunks$from[s]
      j <- chunks$to[s]
      
      ergmito_ptr <- new_ergmito_ptr(
        target_stats  = x[i:j, , drop = FALSE],
        stats_weights = stats_weights[i:j],
        stats_statmat = stats_statmat[i:j],
        target_offset = x[i:j],
        stats_offset  = stats_offset[i:j]
      )
      
      ans[i:j] <- exact_loglik_cpp(ergmito_ptr, params, as_prob = as_prob)
      
    }
  } else {
    
    # In this case, this doesn't change
    for (s in seq_along(chunks$from)) {
      
      i <- chunks$from[s]
      j <- chunks$to[s]
      
      # Creating the model pointer
      ergmito_ptr <- new_ergmito_ptr(
        target_stats  = x[i:j, , drop = FALSE],
        stats_weights = stats_weights,
        stats_statmat = stats_statmat,
        target_offset = target_offset[i:j],
        stats_offset  = stats_offset
      )
      
      ans[i:j] <- exact_loglik_cpp(ergmito_ptr, params, as_prob = as_prob)
      
    }
  }
  
  ans
  
}

# This function uis just used for testing
exact_loglik2 <- function(params, stat0, stats) {
  
  sum(params * stat0) - log(stats$weights %*% exp(stats$statmat %*% params))
  
}

#' @rdname exact_loglik
#' @export
exact_gradient <- function(x, params, ...) UseMethod("exact_gradient")

#' @export
# @rdname exact_loglik
exact_gradient.ergmito_ptr <- function(x, params, ...) {
  exact_gradient_cpp(x, params = params)
}

#' @export
#' @rdname exact_loglik
exact_gradient.default <- function(
  x,
  params,
  stats_weights,
  stats_statmat,
  target_offset = double(nrow(x)),
  stats_offset  = lapply(stats_weights, function(i) double(length(i))),
  ...
) {
  
  # Need to calculate it using chunks of size 200, otherwise it doesn't work(?)
  chunks <- make_chunks(nrow(x), 4e5)
  
  n <- nrow(x)
  
  if (!is.list(stats_weights))
  stats_weights <- list(stats_weights)

  if (!is.list(stats_statmat))
    stats_statmat <- list(stats_statmat)
  
  if (!is.list(stats_offset))
    stats_offset <- list(stats_offset)
  
  # Computing in chunks
  ans <- matrix(0, nrow = length(params), ncol=1L)
  for (s in seq_along(chunks$from)) {
    
    i <- chunks$from[s]
    j <- chunks$to[s]
    
    # Creating the model pointer
    ergmito_ptr <- new_ergmito_ptr(
      target_stats  = x[i:j, , drop = FALSE],
      stats_weights = stats_weights[i:j],
      stats_statmat = stats_statmat[i:j],
      target_offset = target_offset[i:j],
      stats_offset  = stats_offset[i:j]
    )
    
    ans <- ans + exact_gradient_cpp(ergmito_ptr, params)
    
  }
  
  ans
  
}

#' @rdname exact_loglik
#' @export
exact_hessian <- function(
  params,
  stats_weights,
  stats_statmat,
  stats_offset = lapply(stats_weights, function(i) double(length(i)))
) {
  
  # Need to calculate it using chunks of size 200, otherwise it doesn't work(?)
  n <- length(stats_weights)
  chunks <- make_chunks(n, 4e5)
  
  # Checking the weights and stats mat
  if (n == 1) {
    # If only one observation
    
    if (!is.list(stats_weights))
      stats_weights <- list(stats_weights)
    
    if (!is.list(stats_statmat))
      stats_statmat <- list(stats_statmat)
    
    if (!is.list(stats_offset))
      stats_offset <- list(stats_offset)
    
  } else if (n > 1) {
    # If more than 1, then perhaps we need to recycle the values
    
    if (!is.list(stats_weights)) {
      stats_weights <- list(stats_weights)
    } else if (length(stats_weights) != n) {
      stop("length(stats_weights) != nrow(x). When class(stats_weights) == 'list', the number",
           " of elements should match the number of rows in statistics (x).", 
           call. = FALSE)
    }
    
    if (!is.list(stats_statmat)) {
      stats_statmat <- list(stats_statmat)
    } else if (length(stats_statmat) != n) {
      stop("length(statmat) != nrow(x). When class(statmat) == 'list', the number",
           " of elements should match the number of rows in statistics (x).", 
           call. = FALSE)
    }
    
    if (!is.list(stats_offset)) {
      stats_offset <- list(stats_offset)
    } else if (length(stats_offset) != n) {
      stop("length(stats_offset) != nrow(x). When class(stats_offset) == 'list', the number",
           " of elements should match the number of rows in statistics (x).", 
           call. = FALSE)
    }
    
  } else 
    stop("nrow(x) == 0. There are no observed statistics.", call. = FALSE)
  
  # Computing in chunks
  ans <- matrix(0, nrow = length(params), ncol = length(params))
  for (s in seq_along(chunks$from)) {
    
    i <- chunks$from[s]
    j <- chunks$to[s]
    
    ans <- ans + exact_hessian_cpp(
      params,
      stats_weights = stats_weights[i:j],
      stats_statmat = stats_statmat[i:j],
      stats_offset  = stats_offset[i:j]
    )
    
  }
  
  ans
  
}

