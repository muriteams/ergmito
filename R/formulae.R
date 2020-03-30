#' Offset processer
#' @param x An object of class formula
#' @details
#' It will identify if the model has offset terms or not. In the case of
#' having offsets, the resulting object has an attribute indicating which
#' of the terms are offset terms.
#' @return A list of formulas.
#' @noRd
parse_offset <- function(x) {
  
  terms_x <- terms(x)
  
  if (!length(attr(terms_x, "term.labels")))
    stop(
      "Invalid model:\n",
      deparse(x), 
      "\nA offset-only model cannot be fitted.", call. = FALSE)
  
  # Are there any offset terms?
  offsets_x <- attr(terms_x, "offset")
  if (!length(offsets_x))
    return(list(x))
  
  # Updating offset terms
  variables <- rownames(attr(terms_x, "factors"))
  variables_free <- gsub("^offset[(]|[)]$", "", variables[offsets_x])
  
  # We need to memorize this later
  offsets_x <- structure(offsets_x, names = variables[offsets_x])
  variables[offsets_x] <- variables_free

  if (attr(terms_x, "response") > 0) 
    ans <- sprintf("%s ~ %s", variables[1L], variables[-1L])
  else 
    ans <- sprintf(" ~ %s", variables)
  
  ans <- lapply(ans, formula, env = environment(x))
  attr(ans, "ergmito_offset") <- offsets_x

  return(ans)
    
}

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
#' - `vertex_attr` Character vector. Vertex attributes used in the model.
#' - `term_names` Names of the terms used in the model.
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
  target_offset = NULL,
  stats_offset  = NULL,
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
    target_stats <- lapply(1:nrow(target_stats), function(i) target_stats[i, , drop=FALSE])
  
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
  
  # Calculating the support of the sufficient statistics -----------------------
  for (i in 1L:nnets(LHS)) {
    
    # Calculating statistics and weights
    if (!length(stats_statmat[[i]])) {
      
      # Smart thing to do, have I calculated this earlier? ---------------------
      # 1. For now we do it if the model is attr-less model

      # Can we find a match in the previous set?
      matching_net <- NULL
      for (j in 1L:(i-1L)) {
        
        if (i == 1)
          break
        
        # Minimum (and only for now): Have the same size
        if ( same_dist(LHS[[i]], LHS[[j]], vattrs) ) {
          
          matching_net <- j
          break
          
        }
        
      }
        
      # If we calculated this earlier, then copy, else, need to recalc
      if (!is.null(matching_net)) {
        
        stats_statmat[[i]] <- stats_statmat[[matching_net]]
        stats_weights[[i]] <- stats_weights[[matching_net]]
        next
        
      } 
      
      # Computing if we need
      allstats_i <- do.call(
        ergm::ergm.allstats, 
        # We correct for zero obs later
        c(list(formula = model.), dots)
        )
      stats_statmat[[i]] <- allstats_i$statmat
      stats_weights[[i]] <- allstats_i$weights
      
    }
    
  }
  
  # Computing target statistics ------------------------------------------------
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
    
    # Can we compute it directly with ergmito? If not, we default to
    # ergm's summary function
    if (directed && all(model_analysis$names %in% AVAILABLE_STATS())) {
      
      target_stats <- count_stats(model)
      
      # Still need to get it once b/c of naming
      i <- 1
      colnames(target_stats) <- names(summary(model.))
      
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
  
  # Should it be normalized to 0? ----------------------------------------------
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
  
  # Updating the model, if needed, since this could affect offset, we need to
  # preassing information
  if (is.null(target_offset))
    target_offset <- double(nrow(target_stats))
  
  if (is.null(stats_offset))
    stats_offset <- lapply(stats_weights, function(i) double(length(i)))
  
  # Are we updating the model? ------------------------------------------------
  offset_terms <- NULL
  if (length(model_update)) {
    
    # Updating the model (must remove the network to apply the formula)
    model_final      <- as.formula(
      sprintf("~ %s", paste(colnames(target_stats), collapse = " + "))
      )
    model_final      <- stats::update.formula(model_final, model_update)
    
    # Parsing offset terms. This removes the offset() around the term
    # names and then it adds an attribute that tags what are the offset
    # variables so it is easier to extract them from the model.
    model_final <- parse_offset(model_final)
    
    # Updating the target_stats first
    g_attrs <- graph_attributes_as_df(LHS)
    target_stats <- cbind(as.data.frame(target_stats), g_attrs)
    target_stats <- lapply(model_final, stats::model.frame, data = target_stats)
    
    # If we have offset terms, we need to separate them from the stats mat.
    offset_terms <- attr(model_final, "ergmito_offset")
    if (length(offset_terms)) {
      
      target_offset <- do.call(cbind, target_stats[offset_terms])
      target_offset <- rowSums(as.matrix(target_offset))
      target_stats  <- target_stats[-offset_terms]
      
    }
    
    # As a matrix
    target_stats  <- as.matrix(do.call(cbind, target_stats))
    
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
      stats_statmat[[i]] <- lapply(
        model_final,
        stats::model.frame,
        data = stats_statmat[[i]]
        )
        
      
      # Parsing offset terms, if any
      if (length(offset_terms)) {
        
        stats_offset[[i]]  <- do.call(cbind, stats_statmat[[i]][offset_terms])
        stats_offset[[i]]  <- rowSums(as.matrix(stats_offset[[i]]))
        stats_statmat[[i]] <- stats_statmat[[i]][-offset_terms]
        
      }
        
      # And taking it back to a matrix
      stats_statmat[[i]] <- as.matrix(do.call(cbind, stats_statmat[[i]]))
        
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
    stats_statmat = stats_statmat,
    target_offset = target_offset,
    stats_offset  = stats_offset
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
  
  # And the joint hessian
  hess <- function(params, ...) {
    
    exact_hessian(
      params        = params,
      stats_weights = stats_weights,
      stats_statmat = stats_statmat,
      stats_offset  = stats_offset
      )
    
  }
  
  
  structure(
    list(
      loglik        = loglik,
      grad          = grad,
      hess          = hess,
      target_stats  = target_stats,
      stats_weights = stats_weights,
      stats_statmat = stats_statmat,
      target_offset = target_offset,
      stats_offset  = stats_offset,
      model         = stats::as.formula(model, env = env),
      model_update  = model_update,
      model_final   = model_final,
      npars         = ncol(target_stats),
      nnets         = nnets(LHS),
      vertex_attrs  = vattrs,
      term_names    = colnames(target_stats)
    ),
    class = "ergmito_loglik"
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

