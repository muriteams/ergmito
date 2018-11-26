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
  
  environment(model) <- parent.frame()
  
  # Getting the estimates
  if (!length(theta))
    theta <- coef(lergm(model, zeroobs = FALSE))
  
  # Obtaining the network(s) object
  net   <- eval(model[[2]], envir = environment(model))
  terms <- attr(stats::terms(model), "term.labels")
  
  # Checking stats0
  dots <- list(...)
  if (nnets(net) > 1L) {
    
    if (length(dots$zeroobs) && dots$zeroobs)
      warning(
        "The option `zeroobs` was set to FALSE. `zeroobs = TRUE` invalidates the pooled ERGM.",
        call. = FALSE)
    
  } 
  
  dots$zeroobs <- FALSE
  
  # Generating powersets
  ans          <- new.env()
  ans$networks <- lapply(sizes, powerset)
  ans$theta    <- theta
  names(ans$networks) <- sizes
  
  # Computing probabilities
  model. <- stats::update.formula(model, pset ~ .)
  ans$counts <- parallel::mclapply(ans$networks, function(psets) {
    
    # Getting the corresponding powerset
    stats <- matrix(NA, nrow = length(psets), ncol = length(terms), dimnames = list(NULL, terms))
    environment(model.) <- environment()
    
    # Counts
    pset <- psets[[1]]
    S    <- do.call(ergm::ergm.allstats, c(list(formula = model.), dots))
    
    
    if (all(terms %in% c("edges", "mutual")) & Sys.getenv("LERGM_TEST") == "") {
      
      stats[,terms] <- count_stats(psets, terms)
      
    } else {
      
      for (i in seq_along(psets)) {
        
        # Updating the model
        pset   <- psets[[i]]
        
        # In the first iteration we need to compute the statsmat
        stats[i, terms] <- summary(model.)
        
      }
      
    }
    
    c(list(stats = stats), S)
    
  }, mc.cores = mc.cores)
  
  
  # Computing probabilities
  ans$prob <- lapply(sizes, function(s) vector("double", 2^(s*(s-1))))
  names(ans$prob) <- as.character(sizes)
  
  # Function to compute probabilities
  ans$calc_prob <- function(theta = NULL) {
    
    if (!length(theta))
      theta <- ans$theta
    
    for (i in seq_along(sizes))
      ans$prob[[i]] <- exp(exact_loglik(
        x       = ans$counts[[i]]$stats,
        params  = theta,
        weights = ans$counts[[i]]$weights,
        statmat = ans$counts[[i]]$statmat
      ))
      
    invisible()
  }
  
  # Calling the prob function
  ans$calc_prob()
  
  # Sampling function
  ans$sample <- function(n, s, theta = NULL) {
    
    # If no new set of parameters is used, then 
    if (length(theta)) {
      ans$calc_prob(theta)
      on.exit(ans$calc_prob())
    } 
    
    # All should be able to be sampled
    test <- which(!(s %in% sizes))
    if (length(test))
      stop("Some values of `s` are not included in the sampling function.",
           call. = FALSE)
    
    s <- as.character(s)
    
    idx <- seq_along(ans$networks[[s]])
    idx <- sample(idx, n, replace = TRUE, prob = ans$prob[[s]])
    ans$networks[[s]][idx]
      
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

#' @export
#' @importFrom utils ls.str
str.lergm_sampler <- function(object, ...) {
  
  utils::ls.str(object, ...)
  
}
