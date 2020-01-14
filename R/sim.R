replicate_vertex_attr <- function(x, attrname, value) {
  
  for (i in seq_along(x)) 
    for (a in seq_along(attrname))
    network::set.vertex.attribute(x[[i]], attrname = attrname[[a]], value=value[[a]])
  
  x
}

#' ERGMito sampler
#' 
#' Create a sampler object that allows you simulating streams of small networks
#' fast.
#' 
#' @param model A formula.
#' @param theta Named vector. Model parameters.
#' @param x An object of class `ergmito_sampler`.
#' @param sizes Integer vector. Values between 2 to 5 (6 becomes too intensive).
#' @param ncores Integer. Number of processors to use.
#' @param cl An object of class [cluster][parallel::makeCluster].
#' @param force Logical. When `FALSE` (default) will try to use `ergmito`'s stat
#' count functions (see [count_stats]). This means that if one of the requested
#' statistics in not available in `ergmito`, then we will use `ergm` to compute
#' them, which is significantly slower (see details).
#' @param ... Further arguments passed to [ergm::ergm.allstats].
#' 
#' @details 
#' While the \CRANpkg{ergm} package is very efficient, it was not built to do some
#' of the computations required in the ergmito package. This translates in having
#' some of the functions of the package (ergm) with poor speed performance. This
#' led us to "reinvent the wheel" in some cases to speed things up, this includes
#' calculating observed statistics in a list of networks.
#' 
#' @return An environment with the following objects:
#' 
#' - `calc_prob` A function to calculate each graph's probability under the
#'   specified model.
#' - `call` A language object with the call.
#' - `counts` A list with 3 elements: `stats` the sufficient statistics of each network,
#'   `weights` and `statmat` the overall matrices of sufficient statistics used to
#'   compute the likelihood.
#' - `network0` The baseline network used to either fit the model or obtain
#'   attributes.
#' - `networks` A list with the actual sample space of networks.
#' - `prob` A numeric
#' - `sample` A function to draw samples. `n` specifies the number of samples to
#'   draw, `s` the size of the networks, and `theta` the parameter to use to
#'   calculate the likelihoods.
#' - `theta` Named numeric vector with the current values of the model parameters.
#'  
#' 
#' @export
#' @importFrom parallel clusterEvalQ makeCluster makeForkCluster clusterExport
#' @examples 
#' 
#' # We can generate a sampler from a random graph
#' set.seed(7131)
#' ans <- new_rergmito(rbernoulli(4) ~ edges)
#' 
#' # Five samples
#' ans$sample(5, s = 4)
#' 
#' # or we can use some nodal data:
#' data(fivenets)
#' ans <- new_rergmito(
#'   fivenets[[3]] ~ edges + nodematch("female"),
#'   theta = c(-1, 1)
#' )
#' 
#' # Five samples
#' ans$sample(5, s = 4) # All these networks have a "female" vertex attr
#' 
new_rergmito <- function(
  model, 
  theta  = NULL,
  sizes  = NULL,
  cl     = NULL,
  ncores = 1L,
  force  = FALSE,
  ...
  ) {
  
  # environment(model) <- parent.frame()
  # What are the sizes
  if (!length(sizes)) {
    sizes <- nvertex(model)
    sizes <- sort(unique(sizes))
  }
  
  if ((any(sizes) > 7) & !force)
    stop(
      "For simulating networks with more than 7 vertices, you need to set ",
      "force = TRUE.", call. = FALSE
    )
  
  # Checking whether we are dealing with undirected networks
  are_directed <- is_directed(model)
  
  if (!(sum(are_directed) %in% c(nnets(model), 0)))
    stop(
      "All networks must be either directed or undirected. We don't support mixed sets",
      call. = FALSE
      )
  undirected <- !are_directed[1]
  
  # Getting the estimates
  if (!length(theta)) {
    environment(model) <- parent.frame()
    theta <- coef(ergmito(model, zeroobs = FALSE))
  }
  
  # Obtaining the network(s) object
  net   <- eval(model[[2]], envir = environment(model))
  terms <- attr(stats::terms(model), "term.labels")
  
  # Checking stats0
  dots <- list(...)
  dots$zeroobs <- FALSE
  
  # Generating powersets
  ans          <- new.env()
  ans$networks <- lapply(sizes, powerset, directed = !undirected, force = force)
  
  ans$theta    <- theta
  names(ans$networks) <- sizes
  
  # Analyzing formula
  ergm_model       <- analyze_formula(model)
  ergm_model_attrs <- which(sapply(ergm_model$attrnames, length) > 0)
  
  # Capturing attributes (if any)
  sampler_w_attributes <- FALSE
  if (length(ergm_model_attrs) && nnets(net) != 1L) {
    
    stop(
      "Nodal attributes cannot be used when nnets() > 1L and sizes != nvertex().",
      call. = FALSE
      )
    
  } else if (length(ergm_model_attrs) && nnets(net) == 1L) {
    
    sampler_w_attributes <- TRUE
    
    for (a in ergm_model_attrs) {
      ergm_model$attrs[[a]] <-
        if (is.null(ergm_model$attrnames[[a]]))
          numeric(0)
      else
        network::get.vertex.attribute(net, attrname = ergm_model$attrnames[[a]])
    }
  }
  
  if ((Sys.getenv("ERGMITO_TEST") == "") && (!undirected | all(ergm_model$names %in% AVAILABLE_STATS()))) {
    # THE ERGMITO WAY ----------------------------------------------------------
    
    
    # We will use this updated version
    model. <- stats::update.formula(model, net_i ~ .)
    environment(model.) <- environment()
    for (s in names(ans$networks)) {
      
      # Generating all statistics (both have to match!)
      if (nnets(net) == 1 && nvertex(net) == as.integer(s)) {
        net_i <- net
      } else {
        net_i <- rbernoulli(as.integer(s))
      } 

      # Computing all stats
      S <- do.call(ergm::ergm.allstats, c(list(formula = model.), dots))
      
      # Doing parlapply over the networks. Notice that we split the chunks using
      # splitIndices instead of applying the function directly since count_stats
      # receives lists (so multiple networks at the same time).
      if (!is.null(cl) | ncores > 1L) {
        
        # Creating cluster object, if needed
        cl <- parallel::makeCluster(ncores)
        on.exit(tryCatch(parallel::stopCluster(cl), error = function(e) e))
        parallel::clusterEvalQ(cl, library(ergmito))
        parallel::clusterExport(cl, c("ergm_model", "ans"), envir = environment())
        
        ans$counts[[s]] <- parallel::parLapply(
          cl  = cl,
          X   = parallel::splitIndices(length(ans$networks[[s]]), ncores), 
          fun = function(idx) {
            
            # Making room
            smat <- matrix(ncol = length(ergm_model$names), nrow = length(idx))
            for (j in seq_along(ergm_model$names)) 
              smat[,j] <- count_stats(
                X     = ans$networks[[s]][idx],
                terms = ergm_model$names[j],
                # All individuals have the same data
                attrs = replicate(length(idx), ergm_model$attrs[[j]], simplify = FALSE)
              )
            
            # Returning a single matrix with network statistics
            smat
            
          })
        
      } else {
        
        ans$counts[[s]] <- lapply(
          X   = parallel::splitIndices(length(ans$networks[[s]]), ncores),
          FUN = function(idx) {
            
            # Making room
            smat <- matrix(ncol = length(ergm_model$names), nrow = length(idx))
            for (j in seq_along(ergm_model$names)) 
              smat[,j] <- count_stats(
                X     = ans$networks[[s]][idx],
                terms = ergm_model$names[j],
                # All individuals have the same data
                attrs = replicate(length(idx), ergm_model$attrs[[j]], simplify = FALSE)
              )
            
            # Returning a single matrix with network statistics
            smat
            
          })
        
      }
        
        
        
      }
      
      # Adding the results
      ans$counts[[s]] <- do.call(rbind, ans$counts[[s]])
      ans$counts[[s]] <- c(list(stats = ans$counts[[s]]), S)
      
    
  } else {
    # THE ERGM WAY -------------------------------------------------------------
    
    if (force)
      warning_ergmito(
        "To generate this sampler we need to use statnet's ergm functions since",
        " not all the requested statistics are available in ergmito. This will ",
        "take longer...",
        call. = FALSE, immediate. = TRUE
      )
    else
      stop(
        "To generate this sampler we need to use statnet's ergm functions since",
        " not all the requested statistics are available in ergmito. If you",
        " would like to procede, use the option `force = TRUE`", call. = FALSE
      )
    
    # We have to switch this flag so that the other methods (member functions)
    # don't try to turn things into network objects (these will already be)
    # that class of objects
    sampler_w_attributes <- FALSE
    
    # Are we addinig attributes?
    if (nnets(net) == 1 && network::is.network(net)) {
      attrs <- list()
      
      vattrs <- setdiff(network::list.vertex.attributes(net), c("na", "vertex.names"))
      if (length(vattrs)) {
        
        # Getting the attributes
        attrs$vertex.attr      <- lapply(vattrs, network::get.vertex.attribute, x = net)
        attrs$vertex.attrnames <- vattrs
        
        # Adding attributes to the networks
        if (length(ans$neworks) > 1L) {
          warning_ergmito(
            "When `length(size) > 1`, attributes from the networks in `x` cannont",
            " be added (don't know what goes with what). We will skip adding the",
            " observed attributes to the family of networks. This could result",
            " in an error if the model includes nodal attributes.", call. = FALSE)
        } else {
          
          ans$networks[[1L]] <- matrix_to_network(ans$networks[[1L]], directed = !undirected)
          ans$networks[[1L]] <- replicate_vertex_attr(
            x        = ans$networks[[1L]],
            attrname = ergm_model$attrnames,
            value    = ergm_model$attrs
            )
        }
        
      } else
        attrs <- NULL
    }
    
    for (s in names(ans$networks)) {
      
      # All networks in ans$networks should be of class network
      ans$networks[[s]] <- lapply(ans$networks[[s]], network::network, directed = !undirected)
      
      # Computing probabilities
      if (nnets(net) == 1) {
        
        if (nvertex(net) != as.integer(s))
          net_i <- rbernoulli(as.integer(s))
        else 
          net_i <- net
        
      } else if (nnets(net) > 1) {
        
        if (length(ergm_model_attrs))
          stop("Nodal attributes cannot be used when nnets() > 1L and sizes != nvertex().",
               call. = FALSE)
        
        net_i <- rbernoulli(as.integer(s))
        
      }
        
      # Need to make sure we are setting this to be a directed/undirected graph.
      net_i <- network::network(net_i, directed = !undirected)
      
      model. <- stats::update.formula(model, net_i ~ .)
      environment(model.) <- environment()
      i <- 1L
      S <- do.call(ergm::ergm.allstats, c(list(formula = model.), dots))
      
      # Updating the model (again)
      model. <- stats::update.formula(model, ans$networks[[s]][[i]] ~ .)
      
      # Creating cluster object, if needed
      on.exit(tryCatch(parallel::stopCluster(cl), error = function(e) e))
      cl <- parallel::makeCluster(ncores)
      parallel::clusterEvalQ(cl, library(ergmito))
      parallel::clusterExport(cl, c("ans", "model.", "terms"), envir = environment())
      
      ncores <- length(cl)
      
      if (ncores > 1L) {
        
        ans$counts[[s]] <- parallel::parLapply(
          cl  = cl,
          X   = seq_along(ans$networks[[s]]),
          fun = function(i) {
          
          # Updating the environment
          environment(model.) <- environment()
          
          # In the first iteration we need to compute the statsmat
          summary(model.)
          
        })
        
      } else {
        
        ans$counts[[s]] <- lapply(
          X   = seq_along(ans$networks[[s]]),
          FUN = function(i) {
            
            # Updating the environment
            environment(model.) <- environment()
            
            # In the first iteration we need to compute the statsmat
            summary(model.)

          })
        
      }
        
      ans$counts[[s]] <- do.call(rbind, ans$counts[[s]])
      ans$counts[[s]] <- c(list(stats = ans$counts[[s]]), S)
      
    }

  }
  
  # Computing probabilities
  ans$prob <- lapply(sizes, function(s) {
    vector(
      "double",
      2 ^ ( s*(s-1)/ifelse(undirected, 2, 1) )
      )
  })
  names(ans$prob) <- as.character(sizes)
  
  # Function to compute probabilities
  ans$calc_prob <- function(theta = NULL, s = NULL) {
    
    if (!length(s))
      s <- names(ans$prob)
    
    if (!length(theta))
      theta <- ans$theta
    
    for (i in s) {
      ans$prob[[i]] <- exp(exact_loglik(
        x             = ans$counts[[i]]$stats,
        params        = theta,
        stats.weights = ans$counts[[i]]$weights,
        stats.statmat = ans$counts[[i]]$statmat
      ))
      
      # Should be within epsilon
      test <- which(!is.finite(ans$prob[[i]]))
      if (any(test) | abs(sum(ans$prob[[i]]) - 1) > .00001)
        stop(
          "The sum of each network's probability does not equal to one.",
          " This may be due to model parameters that are too far away from 0.",
          " The current value of theta is: [",
          paste(sprintf("%.4f", theta), collapse = ", "),
          "].",
          " You should try with other parameters.",
          call. = FALSE
          )
      
      # Normalizing to be 1
      ans$prob[[i]] <- ans$prob[[i]]/sum(ans$prob[[i]])
      
    }
    
    invisible()
  }
  
  # Calling the prob function
  ans$calc_prob()
  
  # A getter function ----------------------------------------------------------
  ans$get_networks <- function(idx, s) {
    
    s <- as.character(s)
    test <- which(!(s %in% as.character(sizes)))
    if (length(test))
      stop("Some values of `s` are not included in the sampling function:",
           paste0(s[test], collapse = ", "), ".",
           call. = FALSE)
    
    nets <- ans$networks[[s]][idx]
    
    if (sampler_w_attributes) {
      nets <- matrix_to_network(nets)
      replicate_vertex_attr(
        nets,
        attrname = ergm_model$attrnames,
        value    = ergm_model$attrs
      ) 
    } else
      nets
    
  }
  
  # Sampling functions ---------------------------------------------------------
  ans$sample <- function(n, s, theta = NULL, as_indexes = FALSE) {
    
    s <- as.character(s)
    # All should be able to be sampled
    test <- which(!(s %in% as.character(sizes)))
    if (length(test))
      stop("Some values of `s` are not included in the sampling function.",
           paste0(s[test], collapse = ", "), ".",
           call. = FALSE)
    
    # If no new set of parameters is used, then 
    if (length(theta)) {
      oldp <- ans$prob[s]
      ans$calc_prob(theta, s)
      on.exit(ans$prob[s] <- oldp)
    } 
    
    idx <- sample.int(
      n       = length(ans$prob[[s]]),
      size    = n,
      replace = TRUE,
      prob    = ans$prob[[s]],
      useHash = FALSE
    )
    
    if (!as_indexes) 
      ans$get_networks(idx, s)
    else 
      idx
    
  }
  
  # Call
  ans$call      <- match.call()
  ans$network0  <- net
  ans$sizes     <- sizes
  
  structure(
    ans,
    class = "ergmito_sampler"
  )
  
}

#' @export
#' @rdname new_rergmito
#' @param i,j `i` is an integer vector indicating the indexes of the networks to
#' draw, while `j` the corresponding sizes. These need not to be of the same size.
#' @details The indexing method, `[.ergmito_sampler`, allows extracting networks
#' directly by passing indexes. `i` indicates the index of the networks to draw,
#' which go from 1 through `2^(n*(n-1))`, and `j` indicates the requested
#' size. 
#' @return The indexing method `[.ergmito_sampler` returns a named list of length
#' `length(j)`.
`[.ergmito_sampler` <- function(x, i, j, ...) {
  
  # Checking sizes
  j <- as.character(j)
  test <- which(!(j %in% as.character(x$sizes)))
  if (length(test))
    stop(
      "Some values of `j` (requested sizes) are not included in the sampling function: ",
      paste(j[test], collapse = ", "), ".", call. = FALSE
      )
  
  # Sampling networks
  ans <- structure(vector("list", length(j)), names = j)
  for (k in j)
    ans[[j]] <- x$get_networks(i, k)
  
  return(ans)
  
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
