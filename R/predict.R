#' Prediction method for `ergmito` objects
#' 
#' Takes an [ergmito] object and makes prediction at tie level. See details for
#' information regarding its implementation.
#' 
#' @param object An object of class [ergmito].
#' @param newdata New set of networks (or network) on which to make the prediction.
#' @param ... Passed to [new_rergmito], the workhorse of this method.
#' 
#' @details After fitting a model with a small network (or a set of them), we
#' can use the parameter estimates to calculate the likelihood of observing any
#' given tie in the network, this is, the marginal probabilites at the tie level.
#' 
#' In particular, the function takes the full set of networks on the support of
#' the model and adds them up weighted by the probability of observing them as
#' predicted by the ERGM, formally:
#' 
#' \deqn{%
#' \hat{A} = \sum_i \mathbf{Pr}(A = a_i)\times a_i
#' }{%
#' \hat A = \sum_i Pr(A = a[i]) x a[i]
#' }
#' 
#' Where \eqn{\hat{A}}{\hat A} is the predicted adjacency matrix, and
#' \eqn{a_i}{a[i]} is the i-th network in the support of the model. This calculation
#' is done for each individual network used in the model.
#' 
#' @return A list of adjacency matrix of length `nnets(object)` or, if
#' specified `nnets(newdata)`.
#' 
#' @export
#' @examples 
#' data(fivenets)
#' 
#' # bernoulli graph
#' fit <- ergmito(fivenets ~ edges) 
#' 
#' # all ties have the same likelihood
#' # which is roughly equal to:
#' # mean(nedges(fivenets)/(nvertex(fivenets)*(nvertex(fivenets) - 1)))
#' predict(fit) 
#' 
#' 
#' # If we take into account vertex attributes, now the story is different!
#' fit <- ergmito(fivenets ~ edges + nodematch("female")) 
#' 
#' # Not all ties have the same likelihood, since it depends on homophily!
#' predict(fit) 
predict.ergmito <- function(object, newdata = NULL, ...) {
  
  # The best way of implementing this is with the sampler, which already 
  # computes the probability of each sample to be observed.
  
  # Step 1: Get the network from the model
  model  <- stats::formula(object)
  model. <- stats::update.formula(model, networks[[i]] ~ .)
  environment(model.) <- environment()
  
  if (!is.null(newdata)) {
    if (inherits(newdata, "network"))
      networks <- list(newdata)
    else if (inherits(newdata, "matrix")) {
      
      # Square
      if (dim(newdata)[1] != dim(newdata)[2])
        stop("If `newdata` is specified and is of class 'matrix', it should be ",
             "a square matrix.", call. = FALSE)
      
      networks <- list(newdata)
        
    } else if (is.list(newdata) & all(sapply(newdata, inherits, what = "network"))) {
      networks <- newdata
    } else if (is.list(newdata) & all(sapply(newdata, inherits, what = "matrix"))) {
      
      test <- sapply(newdata, dim)
      if (!all(test[1,] == test[2,]))
        stop(
          "When a list of matrices are passed to `newdata`, all must be squared.",
          call. = FALSE)
      
      networks <- newdata
      
    } else
      stop("`newdata` should be either a list of (or single) networks or ",
           "square matrices.", call. = FALSE)
      
  } else 
    networks <- eval(model[[2]], environment(model))
  
  # Step 2: Loop through the networks to generate the predictions:
  predictions <- vector("list", nnets(networks))
  for (i in seq_along(networks)) {
    
    # Checking if we have computed this previously
    for (j in 1L:i) {
      
      if (i == 1L | i == j)
        break
      
      # If compared previously, then assign that value
      equal_networks <-
        same_dist(
          networks[[i]], networks[[j]],
          object$formulae$used_attrs$attr
          # subset(object$formulae$vertex_attrs, type == "vertex")
          )
      if (equal_networks) {
        
        # Put the prediction in the right order of things
        if (length(object$formulae$vertex_attrs)) {
          ord <- attr(equal_networks, "map10")
          predictions[[i]] <- predictions[[j]][ord,][,ord]
        } else {
          predictions[[i]] <- predictions[[j]]
        }
        
        break
      }
    }
    
    # Did we found any matching network, i.e., previously computed quantitites?
    if (!is.null(predictions[[i]]))
      next
    
    # Generating sample, and later on, adding up matrices
    sampler <- new_rergmito(
      model = model.,
      theta = coef(object),
      # sizes = nvertex(networks[[i]]),
      ...
      )
    
    mats  <- sampler$networks
    probs <- sampler$probabilities
    
    predictions[[i]] <- mats[[1L]]*probs[1L]
    for (j in 2L:length(probs))
      predictions[[i]] <- predictions[[i]] + mats[[j]]*probs[j]

  }
  
  predictions
  
}

