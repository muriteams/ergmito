#' Processing formulas in `ergmito`
#' 
#' Analyze formula objects returning the matrices of weights and sufficient 
#' statistics to be used in the model together with the log-likelihood and
#' gradient functions for joint models. 
#' 
#' @param model A formula. The left-hand-side can be either a small network, or
#' a list of networks. 
#' @param gattr_model A formula. Model specification for graph attributes. This
#' is useful when using multiple networks.
#' @param stats_weights,stats_statmat Lists of sufficient statistics and their
#' respective weights.
#' @param target_stats Observed statistics. If multiple networks, then a list, otherwise
#' a named vector (see [ergm::summary_formula]). 
#' @param ... Further arguments passed to [ergm::ergm.allstats].
#' @param env Environment in which `model` should be evaluated.
#' @return A list of class `ergmito_loglik`.
#' 
#' - `loglik` A function. The log-likelihood function.
#' - `grad` A function. The gradient of the model.
#' - `stats_weights`,`stats_statmat` two list of objects as returned by
#' [ergm::ergm.allstats].
#' - `model` A formula. The model passed.
#' - `npars` Integer. Number of parameters.
#' - `nnets` Integer. Number of networks in the model.
#' - `vertex.attr` Character vector. Vertex attributes used in the model.
#' - `term.names` Names of the terms used in the model.
#' 
#' @aliases ergmito_loglik
#' @examples 
#' data(fivenets)
#' model <- ergmito_formulae(fivenets ~ edges + nodematch("female"))
#' print(model)
#' model$loglik(c(-2, 2))
#' @export
ergmito_formulae <- function(
  model,
  model_update  = NULL,
  target_stats  = NULL,
  stats_weights = NULL,
  stats_statmat = NULL,
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
  if (!length(stats_weights) | !length(stats_statmat)) {
    stats_weights <- vector("list", nnets(LHS))
    stats_statmat <- stats_weights
  }
  
  # Has the user passed target statistics?
  if (!length(target_stats))
    target_stats <- vector("list", nnets(LHS))
  else # This is painful, the current way I have this written needs me to pass
       # the target_stats as a list instead of a matrix, which is not the best.
       # For now it should be OK.
    target_stats <- lapply(1:nrow(target_stats), function(i) target_stats[i,])
  
  # Checking types
  test     <- sapply(LHS, inherits, what = "network")
  directed <- TRUE
  if (length(which(test))) {
    test[test] <- test[test] & !sapply(LHS[test], network::is.directed)
    
    test <- which(test)
    
    if ((length(test) != 0) & (length(test) < length(LHS))) {
      
      stop(
        "All networks should be of the same type. Right now, networks ",
        paste(test, collapse = ", "), " are undirected while networks ",
        paste(setdiff(1:length(LHS), test), collapse = ", "), " are directed.",
        call. = FALSE
        )
      
    } else if (length(test)) {
      
      warning_ergmito(
        "ergmito does not fully support undirected graphs (yet). We will ",
        "continue with the estimation process, but simulation has limited supported ",
        "for now.", call. = FALSE
        )
      
      directed <- FALSE
      
    }
    
  }
  
  # Need to improve the speed of this!
  for (i in 1L:nnets(LHS)) {
    
    # Calculating statistics and weights
    if (!length(stats_statmat[[i]])) {
      
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
        stats_statmat[[i]] <- stats_statmat[[matching_net]]
        stats_weights[[i]] <- stats_weights[[matching_net]]
      } else {
        allstats_i <- do.call(
          ergm::ergm.allstats, 
          # We correct for zero obs later
          c(list(formula = model.), dots)
          )
        stats_statmat[[i]] <- allstats_i$statmat
        stats_weights[[i]] <- allstats_i$weights
      }
      
      
    }
    
  }
  
  if (all(sapply(target_stats, length) == 0)) {
    
    # Should we use summary.formula?
    model_analysis <- analyze_formula(model)
    
    # Checking gw terms
    if (any(grepl("^d?gw", model_analysis$names)))
      stop(
      "Currently, geometrically weighted terms are not supported in ergmito.",
      " For more information, see https://github.com/muriteams/ergmito/issues/17.",
      call. = FALSE
      )
    
    if (directed && all(model_analysis$names %in% AVAILABLE_STATS())) {
      
      target_stats <- count_stats(model)
      
      # Still need to get it once b/c of naming
      i <- 1
      colnames(target_stats) <- names(
        summary(model.)
      )
      
    } else {
      for (i in seq_len(nnets(LHS))) {
        # Calculating observed statistics
        if (!length(target_stats[[i]]))
          target_stats[[i]] <- c(summary(model.))

      }
      
      # Coercing objects
      target_stats   <- do.call(rbind, target_stats)
    }
    
  }
  
  if (is.list(target_stats))
    target_stats <- do.call(rbind, target_stats)
  
  # Should it be normalized to 0?
  if (zeroobs)
    for (i in seq_len(nnets(LHS))) {
      stats_statmat[[i]] <- stats_statmat[[i]] - 
        matrix(
          data  = target_stats[i, ],
          nrow  = nrow(stats_statmat[[i]]),
          ncol  = ncol(target_stats),
          byrow = TRUE
        )
      
      target_stats[i, ] <- 0
      
    }
  
  # Updating the model, if needed
  if (length(model_update)) {
    
    # Updating the model (must remove the network to apply the formula)
    model_final      <- as.formula(
      sprintf("~ %s", paste(colnames(target_stats), collapse = " + "))
      )
    model_final      <- stats::update.formula(model_final, model_update)
    
    # Updating the target_stats first
    g_attrs <- graph_attributes_as_df(LHS)
    target_stats <- stats::model.frame(
      formula = model_final,
      data    = cbind(
        as.data.frame(target_stats),
        g_attrs
        )
    )
    # Taking it back to a matrix
    target_stats <- as.matrix(target_stats)
    
    # Doing the same for the rest of the networks
    for (i in seq_len(nnets(LHS))) {
      
      # Merging the data
      stats_statmat[[i]] <- cbind(
        as.data.frame(stats_statmat[[i]]),
        do.call(rbind, replicate(
          nrow(stats_statmat[[i]]),
          g_attrs[i, , drop = FALSE],
          simplify = FALSE
          ))
      )
      
      # Updating the model
      stats_statmat[[i]] <- stats::model.frame(
        formula = model_final,
        data    = stats_statmat[[i]]
        )
      
      # And taking it back to a matrix
      stats_statmat[[i]] <- as.matrix(stats_statmat[[i]])
      
    }
    
    # Reupdating the model (so users look it as is)
    model_final <- stats::update.formula(model, model_update)
    
  } else {
    model_final <- model
  }
  
  # Initializing the pointer
  ergmito_ptr   <- new_ergmito_ptr(
    target_stats  = target_stats,
    stats_weights = stats_weights,
    stats_statmat = stats_statmat
    )
  
  
  # Building joint likelihood function
  loglik <- function(params, ...) {
    
    ans <- sum(exact_loglik(ergmito_ptr, params = params, ...))
    
    # If awfully undefined
    if (!is.finite(ans))
      return(-.Machine$double.xmax * 1e-100)
    else
      return(ans)
    
  }
  
  # Building joint gradient
  grad <- function(params, ...) {
    
    ans <- exact_gradient(ergmito_ptr, params = params, ...)
    
    test <- which(!is.finite(ans))
    if (length(test))
      ans[test] <- sign(ans[test]) * .Machine$double.xmax / 1e200
    
    ans
    
  }
  
  hess <- function(params, ...) {
    
    exact_hessian(
      x = target_stats,
      params = params,
      stats_weights = stats_weights,
      stats_statmat = stats_statmat
      )
    
  }
  
  
  structure(
    list(
      loglik = loglik,
      grad  = grad,
      hess  = hess,
      target_stats  = target_stats,
      stats_weights = stats_weights,
      stats_statmat = stats_statmat,
      model         = stats::as.formula(model, env = env),
      model_update  = model_update,
      model_final   = model_final,
      npars         = ncol(target_stats),
      nnets         = nnets(LHS),
      vertex.attrs  = vattrs,
      term.names    = colnames(target_stats)
    ),
    class="ergmito_loglik"
  )
  
  
  
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

#' This function extracts networks attributes and returns a data
#' frame
#' @noRd
graph_attributes_as_df <- function(net) {
  
  if (inherits(net, "list") && nnets(net) > 1L) {
    net <- matrix_to_network(net)
    return(do.call(rbind, lapply(net, graph_attributes_as_df)))
  }
  
  if (!network::is.network(net))
    net <- matrix_to_network(net)
  
  # Listing attributes, extracting and naming
  netattrs   <- network::list.network.attributes(net)
  ans        <- lapply(netattrs, network::get.network.attribute, x = net)
  names(ans) <- netattrs
  
  # Returning as data.frame
  as.data.frame(ans)
}

gmodel <- function(model, net) {
  
  netattrs   <- network::list.network.attributes(net)
  ans        <- lapply(netattrs, network::get.network.attribute, x = net)
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

