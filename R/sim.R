#' ERGMito sampler
#' 
#' Using 
#' 
#' @param model A formula.
#' @param theta Named vector. Model parameters.
#' @param x An object of class `ergmito_sampler`.
#' @param sizes Integer vector. Values between 2 to 5 (6 becomes too intensive).
#' @param mc.cores Integer. Passed to [parallel::mclapply]
#' @param ... Further arguments passed to [ergm::ergm.allstats].
#' 
#' @return An environment with the following objects:
#' 
#' - `calc_prob`
#' - `call`
#' - `counts`
#' - `networks`
#' - `prob`
#' - `sample` A function to draw samples. `n` specifies the number of samples to
#'   draw, `s` the size of the networks, and `theta` the parameter to use to
#'   calculate the likelihoods.
#' - `theta`
#' 
#' 
#' 
#' @export
#' @importFrom parallel mclapply
new_rergmito <- function(model, theta = NULL, sizes = NULL, mc.cores = 2L,...) {
  
  environment(model) <- parent.frame()
  
  # What are the sizes
  if (!length(sizes)) {
    sizes <- nvertex(model)
    sizes <- sort(unique(sizes))
  }
  
  # Getting the estimates
  if (!length(theta))
    theta <- coef(ergmito(model, zeroobs = FALSE))
  
  # Obtaining the network(s) object
  net   <- eval(model[[2]], envir = environment(model))
  terms <- attr(stats::terms(model), "term.labels")
  
  # Checking stats0
  dots <- list(...)
  
  # Generating powersets
  ans          <- new.env()
  ans$networks <- lapply(sizes, powerset)
  ans$theta    <- theta
  names(ans$networks) <- sizes
  
  # Are we addinig attributes?
  if (nnets(net) == 1 && network::is.network(net)) {
    attrs <- list()
    
    vattrs <- network::list.vertex.attributes(net)
    if (length(vattrs)) {
      
      # Getting the attributes
      attrs$vertex.attr      <- lapply(vattrs, network::get.vertex.attribute, x = net)
      attrs$vertex.attrnames <- vattrs
      
      # Adding attributes to the networks
      if (length(ans$neworks) > 1L) {
        warning(
          "When `length(size) > 1`, attributes from the networks in `x` cannont",
          " be added (don't know what goes with what). We will skip adding the",
          " observed attributes to the family of networks. This could result",
          " in an error if the model includes nodal attributes.", call. = FALSE)
      } else {
        
        for (i in seq_along(ans$networks[[1]]))
          ans$networks[[1L]][[i]] <- network::network(
            x                = ans$networks[[1L]][[i]],
            vertex.attr      = attrs$vertex.attr,
            vertex.attrnames = attrs$vertex.attrnames
          )
      }
      
    } else
      attrs <- NULL
  }
  
  
  
  # Computing probabilities
  model.     <- stats::update.formula(model, psets[[i]] ~ .)
  ans$counts <- parallel::mclapply(ans$networks, function(psets) {
    
    # Getting the corresponding powerset
    stats <- matrix(NA, nrow = length(psets), ncol = length(terms), dimnames = list(NULL, terms))
    environment(model.) <- environment()
    
    # Counts
    i    <- 1L
    S    <- do.call(ergm::ergm.allstats, c(list(formula = model.), dots))
    
    if (all(terms %in% c("edges", "mutual")) & Sys.getenv("ergmito_TEST") == "") {
      
      stats[,terms] <- count_stats(psets, terms)
      
    } else {
      
      # In the first iteration we need to compute the statsmat
      for (i in seq_along(psets)) 
        stats[i, terms] <- summary(model.)
      
    }
    
    c(list(stats = stats), S)
    
  }, mc.cores = mc.cores)
  
  
  # Computing probabilities
  ans$prob <- lapply(sizes, function(s) vector("double", 2^(s*(s-1))))
  names(ans$prob) <- as.character(sizes)
  
  # Function to compute probabilities
  ans$calc_prob <- function(theta = NULL, s = NULL) {
    
    if (!length(s))
      s <- names(ans$prob)
    
    if (!length(theta))
      theta <- ans$theta
    
    for (i in s)
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
    
    s <- as.character(s)
    # All should be able to be sampled
    test <- which(!(s %in% as.character(sizes)))
    if (length(test))
      stop("Some values of `s` are not included in the sampling function.",
           call. = FALSE)
    
    # If no new set of parameters is used, then 
    if (length(theta)) {
      oldp <- ans$prob[s]
      ans$calc_prob(theta, s)
      on.exit(ans$prob[s] <- oldp)
    } 
    
    ans$networks[[s]][
      sample.int(
        n       = length(ans$prob[[s]]),
        size    = n,
        replace = TRUE,
        prob    = ans$prob[[s]],
        useHash = FALSE
        )
      ]
      
  }
  
  # Call
  ans$call <- match.call()
  ans$network0 <- net
  ans$sizes    <- nvertex(net)
  
  structure(
    ans,
    class = "ergmito_sampler"
  )
  
}

#' @export
#' @rdname new_rergmito
print.ergmito_sampler <- function(x, ...) {
  
  cat("ERGMito simulator\n")
  cat("Call   :", deparse(x$call), "\n")
  cat("sample :", deparse(x$sample)[1], "\n")
  
  invisible(x)
  
}

#' @export
#' @importFrom utils ls.str
str.ergmito_sampler <- function(object, ...) {
  
  utils::ls.str(object, ...)
  
}
