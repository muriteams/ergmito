
#' Processing formulas in `lergm`
#' 
#' The little ERGMs R package allows estimating pulled ERGMs by aggregating
#' independent networks together. 
#' 
#' @param model A formula. The left-hand-side can be either a small network, or
#' a list of networks. 
#' @param stats A list.
#' @param obs_stats Observed statistics. If multiple networks, then a list, otherwise
#' a named vector (see [ergm::summary_formula]). 
#' @param ... Further arguments passed to [ergm::ergm.allstats].
#' @param env Environment in which `model` should be evaluated.
#' @return A list of class `lergm_loglik`.
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
lergm_formulae <- function(
  model,
  stats     = NULL,
  obs_stats = NULL,
  env       = parent.frame(),
  ...
  ) {
  
  # Capturing model
  if (!inherits(model, "formula"))
    model <- eval(model, envir = env)

  # What is the first component
  LHS <- eval(model[[2]], envir = env)

  if (inherits(LHS, "list")) {
    
    # Checking stats0
    dots <- list(...)
    if (length(dots$zeroobs) && dots$zeroobs) {
      warning(
        "The option `zeroobs` was set to FALSE. `zeroobs = TRUE` invalidates the pooled ERGM.",
        call. = FALSE)
    } 
    
    dots$zeroobs <- FALSE
    
    # Checking length of stats
    if (!length(stats))
      stats <- vector("list", length(LHS))
    
    if (!length(obs_stats))
      obs_stats <- vector("list", length(LHS))
    
    # Creating one model per model
    f. <- vector("list", length(LHS))
    for (i in seq_along(f.)) {
      
      # Rewriting the model
      model. <- model
      model.[[2]] <- substitute(
        NET[[i.]],
        list(i. = i, NET = model[[2]])
      )
      
      # Getting the functions
      f.[[i]] <- lergm_formulae(
        model     = model.,
        stats     = stats[[i]],
        obs_stats = obs_stats[[i]],
        env       = env,
        ...
        )
    }
    
    # Additive loglike function
    structure(
      list(
        # Joint log-likelihood
        loglik = function(params, stats = NULL) {
          ll <- 0
          for (i in seq_along(f.))
            ll <- ll + f.[[i]]$loglik(params, stats[[i]])
          ll
          },
        # Joint gradient
        grad = function(params, stats = NULL) {
          gr <- vector("numeric", length(params))
          for (i in seq_along(f.))
            gr <- gr + f.[[i]]$grad(params, stats[[i]])
          gr
          },
        obs_stats = lapply(f., "[[", "obs_stats"),
        stats     = lapply(f., "[[", "stats"),
        model     = model,
        npars     = f.[[1]]$npars,
        nnets     = length(f.),
        zeroobs   = dots$zeroobs
        ),
      class="lergm_loglik"
      )
    
    
  } else if (inherits(LHS, "matrix") | inherits(LHS, "network")) {
    
    # Collecting options
    dots <- list(...)
    
    # Network size
    n <- nrow(LHS)
    
    # Baseline statistics, i.e. the zero out. If this is not applied, then
    # the ERGM coefficients should be interpreted in a different fashion. By
    # default in the ERGM package statistics are centered with respect to the
    # observed network. This is critical when estimated te pooled version of 
    # ERGMs.
    #
    # We MUST NOT center the statistics when estimating the pooled version.
    if (!length(dots$zeroobs))
      dots$zeroobs <- TRUE
    
    environment(model) <- env
    
    if (!length(obs_stats))
      obs_stats <- summary(model)
    
    stats0 <- obs_stats
    if (dots$zeroobs)
      stats0[] <- rep(0, length(obs_stats))
    stats0 <- matrix(stats0, nrow=1)
      
    
    # Calculating statistics and weights
    if (!length(stats))
      stats <- ergm::ergm.allstats(formula = model, ...)
    
    # This environment will like in the same place as where the loglikelihood
    # function was created, so it can be access after it.
    originenv <- new.env(parent = emptyenv())
    originenv$stats   <- stats
    
    # Returning the likelihood function
    structure(list(
      loglik = function(params, stats = NULL) {
        
        # Are we including anything 
        if (!length(stats))  
          stats <- originenv$stats
        
        # Computing the log-likelihood
        exact_loglik(params = params, x = stats0, weights = stats$weights, statmat = stats$statmat)
        
      },
      grad  = function(params, stats = NULL) {
        
        # Are we including anything 
        if (!length(stats))  
          stats <- originenv$stats
        
        exact_loglik_gr(params, stats0, stats)
        
      },
      stats = stats,
      obs_stats = obs_stats,
      model = stats::as.formula(model, env = env),
      npars = ncol(originenv$stats$statmat),
      nnets = 1L
    ), class="lergm_loglik")

  } else 
    stop("One of the components is not a matrix `", deparse(model[[2]]),
         "` is of class ", class(LHS), ".", call. = FALSE)
  
  
}

#' @export
print.lergm_loglik <- function(x, ...) {
  
  cat("lergm log-likelihood function\n")
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
  chunks <- make_chunks(nrow(x), 2e5)
  
  ans <- vector("double", nrow(x))
  for (s in seq_along(length(chunks$from))) {
    
    i <- chunks$from[s]
    j <- chunks$to[s]
    
    ans[c((i+1):j)] <- exact_loglik.(x[(i+1):j, ,drop=FALSE], params, weights, statmat)
    
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

