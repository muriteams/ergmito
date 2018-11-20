#' Little ERGM sampler
#' 
#' Using 
#' 
#' @param model A formula.
#' @param theta Named vector. Model parameters.
#' @param x An object of class `lergm_sampler`.
#' @param sizes Integer vector. Values between 2 to 5 (6 becomes too intensive).
#' @param mc.cores Integer. Passed to [parallel::mclapply]
#' @param ... Further arguments passed to [ergm::ergm.allstats].
#' 
#' @export
#' @importFrom parallel mclapply
new_rlergm <- function(model, theta = NULL, sizes = 2:4, mc.cores = 2L,...) {
  
  # Getting the estimates
  if (!length(theta))
    theta <- coef(lergm(model))
  
  # Obtaining the network(s) object
  net <- eval(model[[2]])
  
  # Checking stats0
  dots <- list(...)
  if (nnets(net) > 1L) {
    
    if (length(dots$zeroobs) && dots$zeroobs)
      warning(
        "The option `zeroobs` was set to FALSE. `zeroobs = TRUE` invalidates the pooled ERGM.",
        call. = FALSE)
    
    dots$zeroobs <- FALSE
  } 
  # dots$zeroobs <- FALSE
  
  # Generating powersets
  ans   <- new.env()
  ans$networks <- lapply(sizes, powerset)
  names(ans$networks) <- sizes
  
  # Computing probabilities
  model. <- stats::update.formula(model, pset ~ .)
  ans$prob <- parallel::mclapply(ans$networks, function(psets) {
    
    # Getting the corresponding powerset
    pr <- vector("double", length(psets))
    environment(model.) <- environment()
    
    S <- NULL
    for (i in seq_along(psets)) {
      
      # Updating the model
      pset   <- psets[[i]]
      
      # In the first iteration we need to compute the statsmat
      if (i == 1L)
        S <- do.call(ergm::ergm.allstats, c(list(formula = model.), dots))
      
      pr[i] <- exact_loglik(
        params = theta,
        stat0  = summary(model.),
        stats  = S
        )
      
    }
    
    pr
    
  }, mc.cores = mc.cores)
  
  # Turning into probability
  ans$prob <- lapply(ans$prob, exp)
  names(ans$prob) <- sizes
  
  # Sampling function
  ans$sample <- function(n, s) {
    
    # All should be able to be sampled
    test <- which(!(s %in% sizes))
    if (length(test))
      stop("Some values of `s` are not included in the sampling function.",
           call. = FALSE)
    
    s <- as.character(s)
    
    sample(ans$networks[[s]], n, replace = TRUE, prob = ans$prob[[s]])
      
  }
  
  # Call
  ans$call <- match.call()
  
  structure(
    ans,
    class = "lergm_sampler"
  )
  
}

#' @export
#' @rdname new_rlergm
print.lergm_sampler <- function(x, ...) {
  
  cat("Little ERGM simulator\n")
  cat("Call   :", deparse(x$call), "\n")
  cat("sample :", deparse(x$sample)[1], "\n")
  
  invisible(x)
  
}
