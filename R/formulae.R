
#' Processing formulas in `ergmito`
#' 
#' The ERGMitos R package allows estimating pooled ERGMs by aggregating
#' independent networks together. 
#' 
#' @param model A formula. The left-hand-side can be either a small network, or
#' a list of networks. 
#' @param gattr_model A formula. Model especification for graph attributes. This
#' is useful when using multiple networks.
#' @param stats A list.
#' @param obs_stats Observed statistics. If multiple networks, then a list, otherwise
#' a named vector (see [ergm::summary_formula]). 
#' @param ... Further arguments passed to [ergm::ergm.allstats].
#' @param env Environment in which `model` should be evaluated.
#' @return A list of class `ergmito_loglik`.
#' 
#' - `loglik` A function. The log-likelihood function.
#' - `grad` A function. The gradient of the model.
#' - `stats` If the number of networks is greater than 1, then a list of objects
#' as returned by [ergm::ergm.allstats] (one per network). Otherwise, a single
#' list returned by the same function.
#' - `model` A formula. The model passed.
#' - `npars` Integer. Number of parameters.
#' - `nnets` Integer. Number of networks to estiamte.
#' 
#' @export
ergmito_formulae <- function(
  model,
  gattr_model = NULL,
  stats       = NULL,
  obs_stats   = NULL,
  env         = parent.frame(),
  ...
  ) {
  
  # Collecting extra options
  dots <- list(...)
  if (!length(dots$zeroobs))
    dots$zeroobs <- FALSE

  # Capturing model
  if (!inherits(model, "formula"))
    model <- eval(model, envir = env)
  
  # What is the first component
  LHS <- eval(model[[2]], envir = env)
  
  # Checking the appropriate types
  if (inherits(LHS, "list")) {
    
    # Are all either matrices or networks?
    test <- which(!(sapply(LHS, class) %in% c("matrix", "network")))
    if (length(test))
      stop("One of the components is not a matrix `", deparse(model[[2]]),
         "` is of class ", class(LHS), ".", call. = FALSE)
    
  } else if (inherits(LHS, "matrix") | inherits(LHS, "network")) {
    
    # Nesting into a list. Later we will loop over the elements, so we can
    # apply the same logic regardless of if it is one or more networks
    LHS <- list(LHS)
    
  } else
    stop("LHS of the formula should be either a list of networks or a single ",
         "network.", call. = FALSE)
  
  # We will evaluate the formula in the current environment
  model. <- model
  model. <- stats::update.formula(model., LHS[[i]] ~ .)
  
  formulaeenv <- environment()
  environment(model.) <- formulaeenv
  
  # Checking observed stats and stats
  if (!length(obs_stats))
    obs_stats <- vector("list", nnets(LHS))
  
  if (!length(stats))
    stats <- vector("list", nnets(LHS))
  
  # Need to improve the speed of this!
  for (i in 1L:nnets(LHS)) {
    # Calculating gattrs_model
    
    if (length(gattr_model))
      g <- gmodel(gattr_model, LHS[[i]])[1,]
    else
      g <- NULL
    
    # Calculating observed statistics
    if (!length(stats[[i]]))
      stats[[i]] <- c(summary(model.), g)
    
    # Should it be normalized to 0?
    if (length(dots$zeroobs) && dots$zeroobs)
      stats[[i]][] <- rep(0, length(stats[[i]][]))
    
    # Calculating statistics and weights
    if (!length(obs_stats[[i]]))
      obs_stats[[i]] <- do.call(ergm::ergm.allstats, c(list(formula = model.), dots))
    
    # Adding graph parameters to the statmat
    if (length(g)) {
      obs_stats[[i]]$statmat <- cbind(
        obs_stats[[i]]$statmat,
        matrix(
          stats[[i]][names(g)], nrow = nrow(obs_stats[[i]]$statmat), ncol=length(g),
          byrow = TRUE, dimnames = list(NULL, names(g))
          )
        )
    }
    
  }
  
  # Coercing objects
  stats   <- do.call(rbind, stats)
  weights <- lapply(obs_stats, "[[", "weights")
  statmat <- lapply(obs_stats, "[[", "statmat")
  
  structure(list(
    loglik = function(params, stats = NULL) {
      
      # Are we including anything 
      ans <- if (!length(stats)) {
        exact_loglik(
          params = params, x = formulaeenv$stats,
          weights = formulaeenv$weights, statmat = formulaeenv$statmat
        )
      } else {
        exact_loglik(
          params = params, x = formulaeenv$stats,
          weights = stats$weights, statmat = stats$statmat
        )
      }
      
      max(sum(ans), -.Machine$double.xmax)
      
    },
    grad  = function(params, stats = NULL) {
      
      # # Are we including anything 
      # if (!length(stats))  
      #   stats <- originenv$stats
      # 
      # exact_loglik_gr(params, stats0, stats)
      
    },
    stats     = stats,
    obs_stats = obs_stats,
    model     = stats::as.formula(model, env = env),
    npars     = ncol(stats),
    nnets     = nnets(LHS)
  ), class="ergmito_loglik")
  
  
  
}

gmodel <- function(model, net) {
  
  netattrs <- network::list.network.attributes(net)
  ans <- lapply(netattrs, network::get.network.attribute, x = net)
  names(ans) <- netattrs
    
  stats::model.matrix(
    stats::update.formula(model, ~ 0 + .),
    as.data.frame(ans)
    )
  
  
}

#' @export
print.ergmito_loglik <- function(x, ...) {
  
  cat("ergmito log-likelihood function\n")
  cat("Number of networks: ", x$nnets, "\n")
  cat("Model: ", deparse(x$model), "\n")
  cat("Available elements by using the $ operator:\n")
  cat(sprintf("loglik: %s", deparse(x$loglik)[1]))
  cat(sprintf("grad  : %s", deparse(x$grad)[1]))
  
  invisible(x)
}

#' Vectorized calculation of ERGM exact loglikelihood
#' 
#' This function can be compared to [ergm::ergm.exact] with the statistics
#' centered at `x`, the observed statistic.
#' 
#' @param x Matrix. Observed statistics
#' @param params Numeric vector. Parameter values of the model.
#' @param weights,statmat Vector and Matrix as returned by [ergm::ergm.allstats].
#' @export
exact_loglik <- function(x, params, weights, statmat) {
  
  # Need to calculate it using chunks of size 200, otherwise it doesn't work(?)
  chunks <- make_chunks(nrow(x), 4e5)
  
  n <- nrow(x)
  
  # Checking the weights and stats mat
  if (n == 1) {
    # If only one observation
    
    if (!is.list(weights))
      weights <- list(weights)
    
    if (!is.list(statmat))
      statmat <- list(statmat)
    
  } else if (n > 1) {
    # If more than 1, then perhaps we need to recycle the values
    
    if (!is.list(weights)) {
      weights <- list(weights)
    } else if (length(weights) != n) {
      stop("length(weights) != nrow(x). When class(weights) == 'list', the number",
           " of elements should match the number of rows in statistics (x).", 
           call. = FALSE)
    }
    
    if (!is.list(statmat)) {
      statmat <- list(statmat)
    } else if (length(statmat) != n) {
      stop("length(statmat) != nrow(x). When class(statmat) == 'list', the number",
           " of elements should match the number of rows in statistics (x).", 
           call. = FALSE)
    }
    
  } else 
    stop("nrow(x) == 0. There are no observed statistics.", call. = FALSE)
  

  # Computing in chunks
  ans <- vector("double", n)
  for (s in seq_along(chunks$from)) {
    
    i <- chunks$from[s]
    j <- chunks$to[s]
    
    ans[i:j] <- exact_loglik.(x[i:j, ,drop=FALSE], params, weights, statmat)
    
  }
  
  ans
  
}

# This function uis just used for testing
exact_loglik2 <- function(params, stat0, stats) {
  
  sum(params * stat0) - log(stats$weights %*% exp(stats$statmat %*% params))
  
}

exact_loglik_gr <- function(params, stat0, stats) {
  
  exp_sum  <- exp(stats$statmat %*% params)
  
  stat0 - 
    1/log(stats$weights %*% exp_sum)[1]*(t(stats$statmat) %*% (exp_sum*stats$weights))
  
  
}

