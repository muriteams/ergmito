
#' Processing formulas in `ergmito`
#' 
#' Analyze formula objects returning the matrices of weights and sufficient 
#' statistics to be used in the model together with the log-likelihood and
#' gradient functions for joint models. 
#' 
#' @param model A formula. The left-hand-side can be either a small network, or
#' a list of networks. 
#' @param gattr_model A formula. Model especification for graph attributes. This
#' is useful when using multiple networks.
#' @param stats.weights,stats.statmat Lists of sufficient statistics and their
#' respective weights.
#' @param target.stats Observed statistics. If multiple networks, then a list, otherwise
#' a named vector (see [ergm::summary_formula]). 
#' @param ... Further arguments passed to [ergm::ergm.allstats].
#' @param env Environment in which `model` should be evaluated.
#' @return A list of class `ergmito_loglik`.
#' 
#' - `loglik` A function. The log-likelihood function.
#' - `grad` A function. The gradient of the model.
#' - `stats.weights`,`stats.statmat` two list of objects as returned by
#' [ergm::ergm.allstats].
#' - `model` A formula. The model passed.
#' - `npars` Integer. Number of parameters.
#' - `nnets` Integer. Number of networks to estiamte.
#' - `vertex.attr` Character vector. Vertex attributes used in the model.
#' - `term.names` Names of the terms used in the model.
#' 
#' @aliases ergmito_loglik
#' @export
ergmito_formulae <- function(
  model,
  gattr_model   = NULL,
  target.stats  = NULL,
  stats.weights = NULL,
  stats.statmat = NULL,
  env           = parent.frame(),
  ...
  ) {
  
  # Collecting extra options
  dots <- list(...)
  if (length(dots$zeroobs) && dots$zeroobs) {
    zeroobs <- TRUE
  } else {
    zeroobs <- FALSE
  }
  dots$zeroobs <- FALSE

  # Capturing model
  if (!inherits(model, "formula"))
    model <- eval(model, envir = env)
  
  # What is the first component
  LHS <- eval(model[[2]], envir = env)
  
  # Checking whether this model has attributes on it or not
  vattrs <- attr(model_has_attrs(model), "anames")
  
  # Checking the appropriate types
  if (inherits(LHS, "list")) {
    
    # Are all either matrices or networks?
    test <- which(!(
      sapply(LHS, inherits, what = "matrix") | sapply(LHS, inherits, what = "network")
      ))
    if (length(test))
      stop(
        "One of the components is not a matrix `", deparse(model[[2]]),
        "` is of class ",
        paste(sapply(LHS, class)[test], collapse = ", "),
        ".",
        call. = FALSE)
    
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
  
  environment(model.) <- environment()
  
  # Checking observed stats and stats
  if (!length(stats.weights) | !length(stats.statmat)) {
    stats.weights <- vector("list", nnets(LHS))
    stats.statmat <- stats.weights
  }
  
  # Has the user passed target statistics?
  if (!length(target.stats))
    target.stats <- vector("list", nnets(LHS))
  else # This is painful, the current way I have this written needs me to pass
       # the target.stats as a list instead of a matrix, which is not the best.
       # For now it should be OK.
    target.stats <- lapply(1:nrow(target.stats), function(i) target.stats[i,])
  
  # Checking types
  test <- sapply(LHS, inherits, what = "network")
  if (length(which(test))) {
    test[test] <- test[test] & !sapply(LHS[test], network::is.directed)
    
    test <- which(test)
    if (length(test))
      stop(
        "One or more networks in the model are undirected. Currently, ",
        "undirected networks are not supported by ergmito (but will, in the ",
        "future). The following networks are marked as undirected: ",
        paste(test, collapse = ", "), ".", call. = FALSE
        )
    
  }
  
  # Need to improve the speed of this!
  for (i in 1L:nnets(LHS)) {
    
    # Calculating statistics and weights
    if (!length(stats.statmat[[i]])) {
      
      # Smart thing to do, have I calculated this earlier? ---------------------
      # 1. For now we do it if the model is attr-less model
      matching_net <- NULL
      if (i > 1L) {
        
        # Can we find a match in the previous set?
        for (j in 1L:(i-1L)) {
          
          # Minimum (and only for now): Have the same size
          if ( same_dist(LHS[[i]], LHS[[j]], vattrs) ) {
            matching_net <- j
            break
          }
          
        }
        
      }
      
      if (!is.null(matching_net)) {
        stats.statmat[[i]] <- stats.statmat[[matching_net]]
        stats.weights[[i]] <- stats.weights[[matching_net]]
      } else {
        allstats_i <- do.call(
          ergm::ergm.allstats, 
          # We correct for zero obs later
          c(list(formula = model.), dots)
          )
        stats.statmat[[i]] <- allstats_i$statmat
        stats.weights[[i]] <- allstats_i$weights
      }
      
      
    }
    
  }
  
  if (all(sapply(target.stats, length) == 0)) {
    
    # Should we use summary.formula?
    model_analysis <- analyze_formula(model)
    if (all(model_analysis$names %in% AVAILABLE_STATS())) {
      
      target.stats <- count_stats(model)
      
      # Still need to get it once b/c of naming
      i <- 1
      colnames(target.stats) <- names(
        summary(model.)
      )
      
    } else {
      for (i in seq_len(nnets(LHS))) {
        # Calculating observed statistics
        if (!length(target.stats[[i]]))
          target.stats[[i]] <- c(summary(model.))

      }
      
      # Coercing objects
      target.stats   <- do.call(rbind, target.stats)
    }
    
  }
  
  g <- vector("list", nnets(LHS))
  for (i in seq_len(nnets(LHS))) {
    
    # Calculating gattrs_model
    if (length(gattr_model))
      g[[i]] <- gmodel(gattr_model, LHS[[i]])
    
    # Adding graph level attributes
    # Adding graph parameters to the statmat
    if (length(g[[i]])) {
      
      stats.statmat[[i]] <- cbind(
        stats.statmat[[i]],
        matrix(
          data     = g[[i]],
          nrow     = nrow(stats.statmat[[i]]),
          ncol     = length(g[[i]]),
          byrow    = TRUE,
          dimnames = list(NULL, names(g[[i]]))
        )
      )
    }
    
  }
  
  if (is.list(target.stats))
    target.stats <- do.call(rbind, target.stats)
  
  if (all(sapply(g, length) != 0))
    target.stats <- cbind(target.stats, do.call(rbind, g))
  
  # Should it be normalized to 0?
  if (zeroobs)
    for (i in seq_len(nnets(LHS))) {
      stats.statmat[[i]] <- stats.statmat[[i]] - 
        matrix(
          data  = target.stats[i, ],
          nrow  = nrow(stats.statmat[[i]]),
          ncol  = ncol(target.stats),
          byrow = TRUE
          )
      
      target.stats[i, ] <- 0
      
    }
  
  
  structure(list(
    loglik = function(params, stats.weights, stats.statmat, target.stats, ncores = 1L) {
      
      ans <- sum(exact_loglik(
        params        = params,
        x             = target.stats,
        stats.weights = stats.weights,
        stats.statmat = stats.statmat,
        ncores        = ncores
      ))
      
      # If awfully undefined
      if (!is.finite(ans))
        return(-.Machine$double.xmax * 1e-100)
      else
        return(ans)
      # max(ans, -.Machine$double.xmax/1e100)
      
    },
    grad  = function(params, stats.weights, stats.statmat, target.stats, ncores = 1L) {

      ans <- exact_gradient(
        params        = params,
        x             = target.stats,
        stats.weights = stats.weights,
        stats.statmat = stats.statmat,
        ncores        = ncores
      )

      test <- which(!is.finite(ans))
      if (length(test))
        ans[test] <- sign(ans[test]) * .Machine$double.xmax / 1e200
      
      ans
      
    },
    # hess  = function(params, stats.weights, stats.statmat, target.stats, ncores = 1L) {
    #   
    #   ans <- exact_hessian(
    #     params        = params,
    #     x             = target.stats,
    #     stats.weights = stats.weights,
    #     stats.statmat = stats.statmat,
    #     ncores        = ncores
    #   )
    #   
    #   test <- which(!is.finite(ans))
    #   if (length(test))
    #     ans[test] <- sign(ans[test]) * .Machine$double.xmax / 1e200
    #   
    #   ans
    #   
    # },
    target.stats  = target.stats,
    stats.weights = stats.weights,
    stats.statmat = stats.statmat,
    model         = stats::as.formula(model, env = env),
    npars         = ncol(target.stats),
    nnets         = nnets(LHS),
    vertex.attrs  = vattrs,
    term.names    = colnames(target.stats)
  ), class="ergmito_loglik")
  
  
  
}

#' Test whether the model terms list attributes
#' 
#' It simply looks for the regex pattern [(].*\".+\".*[])] in the formula
#' terms.
#' @param x A [stats::formula].
#' @noRd
model_has_attrs <- function(x) {
  
  if (!inherits(x, "formula"))
    stop("`x` must be a formula.", call. = FALSE)
  
  trms <- stats::terms(x)
  
  if (any(grepl("[(].*\".+\".*[)]", attr(trms, "term.labels")))) {
    
    # Listing attribute names
    pat    <- "(?<=\")(.+)(?=\")|(?<=\')(.+)(?=\')"
    anames <- NULL
    for (a in attr(trms, "term.labels")) {
      m      <- regexpr(pat, a, perl = TRUE)
      anames <- c(anames, regmatches(a, m))
    }
    
    return(structure(TRUE, anames = unique(anames)))
    
  } else
    return(structure(FALSE, anames = NULL))
  
}

gmodel <- function(model, net) {
  
  netattrs <- network::list.network.attributes(net)
  ans <- lapply(netattrs, network::get.network.attribute, x = net)
  names(ans) <- netattrs
    
  ans <- stats::model.matrix(
    stats::update.formula(model, ~ 0 + .),
    as.data.frame(ans)
  )
  
  structure(
    ans[1, , drop=TRUE],
    names = colnames(ans)
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
#' @param stats.weights,stats.statmat Vector and Matrix as returned by [ergm::ergm.allstats].
#' @param ncores Integer scalar. Number of cores to use with OpenMP (if available).
#' @export
exact_loglik <- function(x, params, stats.weights, stats.statmat, ncores = 1L) {
  
  # Need to calculate it using chunks of size 200, otherwise it doesn't work(?)
  chunks <- make_chunks(nrow(x), 4e5)
  
  n <- nrow(x)
  
  # Checking the weights and stats mat
  if (n == 1) {
    # If only one observation
    
    if (!is.list(stats.weights))
      stats.weights <- list(stats.weights)
    
    if (!is.list(stats.statmat))
      stats.statmat <- list(stats.statmat)
    
  } else if (n > 1) {
    # If more than 1, then perhaps we need to recycle the values
    
    if (!is.list(stats.weights)) {
      stats.weights <- list(stats.weights)
    } else if (length(stats.weights) != n) {
      stop("length(stats.weights) != nrow(x). When class(stats.weights) == 'list', the number",
           " of elements should match the number of rows in statistics (x).", 
           call. = FALSE)
    }
    
    if (!is.list(stats.statmat)) {
      stats.statmat <- list(stats.statmat)
    } else if (length(stats.statmat) != n) {
      stop("length(stats.statmat) != nrow(x). When class(stats.statmat) == 'list', the number",
           " of elements should match the number of rows in statistics (x).", 
           call. = FALSE)
    }
    
  } else 
    stop("nrow(x) == 0. There are no observed statistics.", call. = FALSE)
  

  # Computing in chunks
  ans <- vector("double", n)
  if (length(stats.weights) > 1L) {
    for (s in seq_along(chunks$from)) {
      
      i <- chunks$from[s]
      j <- chunks$to[s]
      
      
      ans[i:j] <- exact_loglik.(
        x[i:j, , drop = FALSE],
        params,
        stats_weights = stats.weights[i:j],
        stats_statmat = stats.statmat[i:j],
        ncores        = ncores
        )
      
    }
  } else {
    for (s in seq_along(chunks$from)) {
      
      i <- chunks$from[s]
      j <- chunks$to[s]
      
      ans[i:j] <- exact_loglik.(
        x[i:j, ,drop=FALSE],
        params,
        stats_weights = stats.weights,
        stats_statmat = stats.statmat,
        ncores        = ncores
        )
      
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
exact_gradient <- function(x, params, stats.weights, stats.statmat, ncores = 1L) {
  
  # Need to calculate it using chunks of size 200, otherwise it doesn't work(?)
  chunks <- make_chunks(nrow(x), 4e5)
  
  n <- nrow(x)
  
  # Checking the weights and stats mat
  if (n == 1) {
    # If only one observation
    
    if (!is.list(stats.weights))
      stats.weights <- list(stats.weights)
    
    if (!is.list(stats.statmat))
      stats.statmat <- list(stats.statmat)
    
  } else if (n > 1) {
    # If more than 1, then perhaps we need to recycle the values
    
    if (!is.list(stats.weights)) {
     stats.weights <- list(stats.weights)
    } else if (length(stats.weights) != n) {
      stop("length(stats.weights) != nrow(x). When class(stats.weights) == 'list', the number",
           " of elements should match the number of rows in statistics (x).", 
           call. = FALSE)
    }
    
    if (!is.list(stats.statmat)) {
      stats.statmat <- list(stats.statmat)
    } else if (length(stats.statmat) != n) {
      stop("length(statmat) != nrow(x). When class(statmat) == 'list', the number",
           " of elements should match the number of rows in statistics (x).", 
           call. = FALSE)
    }
    
  } else 
    stop("nrow(x) == 0. There are no observed statistics.", call. = FALSE)
  
  
  # Computing in chunks
  ans <- matrix(0, nrow = length(params), ncol=1L)
  for (s in seq_along(chunks$from)) {
    
    i <- chunks$from[s]
    j <- chunks$to[s]
    
    ans <- ans + exact_gradient.(
      x[i:j, ,drop=FALSE],
      params,
      stats_weights = stats.weights[i:j],
      stats_statmat = stats.statmat[i:j],
      ncores        = ncores
      )
    
  }
  
  ans
  
}

#' @rdname exact_loglik
#' @export
exact_hessian <- function(x, params, stats.weights, stats.statmat, ncores = 1L) {
  
  # Need to calculate it using chunks of size 200, otherwise it doesn't work(?)
  chunks <- make_chunks(nrow(x), 4e5)
  
  n <- nrow(x)
  
  # Checking the weights and stats mat
  if (n == 1) {
    # If only one observation
    
    if (!is.list(stats.weights))
      stats.weights <- list(stats.weights)
    
    if (!is.list(stats.statmat))
      stats.statmat <- list(stats.statmat)
    
  } else if (n > 1) {
    # If more than 1, then perhaps we need to recycle the values
    
    if (!is.list(stats.weights)) {
      stats.weights <- list(stats.weights)
    } else if (length(stats.weights) != n) {
      stop("length(stats.weights) != nrow(x). When class(stats.weights) == 'list', the number",
           " of elements should match the number of rows in statistics (x).", 
           call. = FALSE)
    }
    
    if (!is.list(stats.statmat)) {
      stats.statmat <- list(stats.statmat)
    } else if (length(stats.statmat) != n) {
      stop("length(statmat) != nrow(x). When class(statmat) == 'list', the number",
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
    
    ans <- ans + exact_hessian.(
      x[i:j, ,drop=FALSE],
      params,
      stats_weights = stats.weights[i:j],
      stats_statmat = stats.statmat[i:j],
      ncores        = ncores
    )
    
  }
  
  ans
  
}

# inline arma::colvec exact_gradienti(
#   const arma::rowvec & x,
#   const arma::colvec & params,
#   const arma::rowvec & stats_weights,
#   const arma::mat    & stats_statmat
# ) {
#   
#   return x.t() - (stats_statmat.t() * (stats_weights.t() % exp(stats_statmat * params)))/
#     kappa(params, stats_weights, stats_statmat);
#   
# }

# 
# exact_loglik_gr <- function(params, stat0, stats) {
#   
#   exp_sum  <- exp(stats$statmat %*% params)
#   
#   RHS <- 1/log(stats$weights %*% exp_sum) *
#     (t(stats$statmat) %*% (exp_sum*as.vector(stats$weights)))
#   RHS <- matrix(RHS, nrow=nrow(stat0), ncol=ncol(stat0), byrow = TRUE)
#   
#   stat0 - RHS
#   
# }

