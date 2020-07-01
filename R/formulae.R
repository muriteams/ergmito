#' Processing formulas in `ergmito`
#' 
#' Analyze formula objects returning the matrices of weights and sufficient 
#' statistics to be used in the model together with the log-likelihood and
#' gradient functions for joint models. 
#' 
#' @param model A formula. The left-hand-side can be either a small network, or
#' a list of networks. 
#' @param stats_weights,stats_statmat Lists of sufficient statistics and their
#' respective weights.
#' @param target_stats Observed statistics. If multiple networks, then a list, otherwise
#' a named vector (see [ergm::summary_formula]). 
#' @param model_update A formula. If specified, the after computing the
#' sufficient statistics (observed and support), the model is updated using
#' [stats::model.frame()]. This includes processing offset terms.
#' @template offset
#' @param ... Further arguments passed to [ergm::ergm.allstats].
#' @param env Environment in which `model` should be evaluated.
#' @details 
#' One of the main advantages of been able to compute exact likelihoods is that
#' we can build arbitrarily complex models in the same way that we would do in
#' the context of Generalized Linear Models, this is, adding offset terms,
#' interaction effects, or transformations of statistics without much effort.
#' 
#' In particular, if the user passes a formula via `model_update`, the 
#' cannonical additive ERGM can be modified to include other terms, for example,
#' if we wanted to add an interaction effect of the `nodematch("age")` with
#' network size, we can simply type
#' 
#' ```
#' model_update = ~ . + I(nodematch.age * n)
#' ```
#' 
#' The [I()] function allows operating over variables in the model, in this case,
#' we took the `nodematch.age` variable (which is the name that [ergm::ergm()] 
#' assigns to it after computing the sufficient statistics) and multiplied it by
#' `n`, which is the network size (this variable is included by default).
#' 
#' By default, the ergm package calculates up to 2^16 unique values for the
#' vector of sufficient statistics. This results in issues if the user tries to
#' fit a model with too heterogenous networks or sets of attributes. To deal 
#' with this it suffices with adding the option `maxNumChangeStatVectors` in
#' the ergmito call, e.g.:
#' 
#' ```
#' # Networks of size 5 have up to 2^20 unique sets of sufficient statistics
#' ergmito(..., maxNumChangeStatVectors = 2^20)
#' ```
#' 
#' See more in ?[ergm::ergm.allstats].
#' 
#' @return A list of class `ergmito_loglik`.
#' 
#' - `loglik` A function. The log-likelihood function.
#' - `grad` A function. The gradient of the model.
#' - `stats_weights`,`stats_statmat` two list of objects as returned by
#' [ergm::ergm.allstats].
#' - `target_offset`,`stats_offset` A vector of offset terms and a list of
#' vectors of offset terms, one for the target stats and the other for the
#' support of the sufficient statistics (defaults to 0).
#' - `model` A formula. The model passed.
#' - `npars` Integer. Number of parameters.
#' - `nnets` Integer. Number of networks in the model.
#' - `vertex_attr` Character vector. Vertex attributes used in the model.
#' - `term_names` Names of the terms used in the model.
#' 
#' @aliases ergmito_loglik
#' @examples 
#' data(fivenets)
#' model0 <- ergmito_formulae(fivenets ~ edges + nodematch("female"))
#' print(model0)
#' model0$loglik(c(-2, 2))
#' 
#' # Model with interaction effects and an offset term
#' model1 <- ergmito_formulae(
#'   fivenets ~ edges + nodematch("female"),
#'   model_update = ~ . + offset(edges) + I(edges * nodematch.female)
#' )
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
  
  # Capturing model
  if (!inherits(model, "formula"))
    model <- eval(model, envir = env)
  
  # What is the first component
  LHS <- eval(model[[2]], envir = env)
  
  # Analyzing the formula
  model_analysis <- analyze_formula(model, LHS)

  # Checking if statmat weights and offsets are passed, then we assume the user
  # knows what is doing, so we should skip all the checks
  test <- c(
    target_stats  = is.null(target_stats),
    stats_weights = is.null(stats_weights),
    stats_statmat = is.null(stats_statmat),
    target_offset = is.null(target_offset),
    stats_offset  = is.null(stats_offset)
  )
  
  if (all(test)) {
    
    # Collecting extra options
    dots <- list(...)
    if (length(dots$zeroobs) && dots$zeroobs) {
      zeroobs <- TRUE
    } else {
      zeroobs <- FALSE
    }
    dots$zeroobs <- FALSE
    
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
          if ( same_dist(LHS[[i]], LHS[[j]], model_analysis$all_attrs$attr) ) {
            
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
        allstats_i <- tryCatch(do.call(
          ergm::ergm.allstats, 
          # We correct for zero obs later
          c(list(formula = model.), dots)
        ), error = function(e) e
        )
        
        # We need to cach for the error that shows in the function.
        if (inherits(allstats_i, "error") | is.null(allstats_i)) {
          
          if (!is.null(allstats_i) && grepl("initialization.+not found", allstats_i$message)) {
            stop(
              "The term you are trying to use was not found. The following is ",
              "the full error message returned by the ergm package:",
              allstats_i$message, call. = FALSE
              )
          }
          
          msg <- if (is.null(allstats_i)) "use force = TRUE."
          else allstats_i$message
          
          stop(
            "The function ergm::ergm.allstats returned with an error. Most of ",
            "time this error means that the model you are trying to comput is ",
            "too large. If you are sure you want to continue, add the option ",
            "maxNumChangeStatVectors = 2^20 if you are using directed graphs of",
            "size 5, or try setting force = TRUE. For more info see ",
            "help(\"ergm.allstats\", \"ergm\"). Here is ",
            "the error reported by the function:\n",
            paste0(msg, collapse = "\n"),
            call. = FALSE
            )
        }
        
        
        stats_statmat[[i]] <- allstats_i$statmat
        stats_weights[[i]] <- allstats_i$weights
        
      }
      
    }
    
    # Computing target statistics ------------------------------------------------
    
    # Checking gw terms
    if (any(grepl("^d?gw", model_analysis$term_names)))
      stop(
        "Currently, geometrically weighted terms are not supported in ergmito.",
        " For more information, see https://github.com/muriteams/ergmito/issues/17.",
        call. = FALSE
      )
    
    # Can we compute it directly with ergmito? If not, we default to
    # ergm's summary function
    if (directed && all(model_analysis$term_names %in% AVAILABLE_STATS())) {
      
      target_stats <- count_stats(model)
      
      # Still need to get it once b/c of naming
      i <- 1
      colnames(target_stats) <- names(ergm::summary_formula(model.))
      
    } else {
      
      for (i in seq_len(nnets(LHS))) {
        # Calculating observed statistics
        if (!length(target_stats[[i]]))
          target_stats[[i]] <- c(ergm::summary_formula(model.))
        
      }
      
      # Coercing objects
      target_stats   <- do.call(rbind, target_stats)
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
    
    # Are we updating the model? ------------------------------------------------
    model_update <- model_update_parser(model_update)
    
    graph_attrs <- graph_attributes_as_df(LHS)
    model_frame <- model_frame_ergmito(
      formula        = model,
      formula_update = model_update, 
      data           = target_stats,
      g_attrs.       = graph_attrs
    )
    
    target_stats  <- model_frame$stats
    target_offset <- model_frame$offsets
    model_final   <- model_frame$model
    
    
    # Doing the same for the rest of the networks
    stats_offset <- vector("list", length(stats_weights))
    
    for (i in seq_len(nnets(LHS))) {
      
      # Merging the data
      model_frame <- model_frame_ergmito(
        formula        = model,
        formula_update = model_update, 
        data           = stats_statmat[[i]],
        g_attrs.       = graph_attrs[i, , drop = FALSE]
      )
      
      stats_statmat[[i]] <- model_frame$stats
      stats_offset[[i]]  <- model_frame$offsets
      
    }
    
  } else if (sum(!test) != length(test)) { # Here is not all, but some
    
    stop(
      "When passing target_stats, statmat_stats, statmat_weights, target_offset, ",
      "statmat_offset either all of none is specified. Right now only : -",
      paste(names(test[!test]), collapse = "-, -"), "- are specified.",
      call. = FALSE
      )
    
  } else {
    # We only need to build the model_final object
    model_final <- if (!is.null(model_update))
      stats::update.formula(model, model_update)
    else
      model
  }
  
  # Checking offsets, if some offsets go to -Inf then it means that the values
  # shouldn't be included in the model, i.e., we are truncating the sample space
  # this should apply for all that's done forward, so we need to subset the values
  # overall
  
  # First, the space of sufficient statistics
  for (i in seq_along(stats_statmat)) {
    
    inf_test    <- is.infinite(stats_offset[[i]])
    
    if (any(inf_test & (sign(stats_offset[[i]]) > 0)))
      stop(
        "Some offset terms in the support have +Inf, thus the log-likelihood ",
        "function is not well defined.",
        call. = FALSE
        )
    
    inf_test <- which(inf_test)
    
    # If all well defined, then go next
    if (!length(inf_test))
      next
    else if (length(inf_test) == nrow(stats_statmat[[i]]))
      stop(
        "All the offset terms on the support go to -Inf, thus the log-likelihood ",
        "is not well defined.",
        call. = FALSE
        )
    
    # Subsetting
    stats_statmat[[i]] <- stats_statmat[[i]][-inf_test, , drop = FALSE]
    stats_weights[[i]] <- stats_weights[[i]][-inf_test]
    stats_offset[[i]]  <- stats_offset[[i]][-inf_test]
    
  }
  
  # Checking the offset of the target stats
  excluded <- is.infinite(target_offset)
  if (any(excluded & (sign(target_offset) > 0)))
    stop(
      "Some offset terms have +Inf, thus the log-likelihood function is not ",
      "well defined.",
      call. = FALSE
    )
  
  excluded <- which(excluded)
  if (length(excluded) == length(target_offset)) {
    warning_ergmito(
      "All of the observed offset terms have -Inf, thus the log-likelihood function ",
      "describes a bernoulli graph.",
      call. = FALSE
    )
    target_offset <- 0
    target_stats  <- matrix(
      0, nrow = 1, ncol = ncol(target_stats),
      dimnames = list(NULL, colnames(target_stats))
      )
    
  } else if (length(excluded) > 0) {
    
    # We will need to further remove observations from the suffstats support
    target_offset <- target_offset[-excluded]
    target_stats  <- target_stats[-excluded, , drop = FALSE]

    stats_statmat <- stats_statmat[-excluded]
    stats_weights <- stats_weights[-excluded]
    stats_offset  <- stats_offset[-excluded]
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
  loglik <- function(params, ..., as_prob = FALSE, total = TRUE) {
    
    if (total) {
    
      ans <- sum(exact_loglik(ergmito_ptr, params = params, ..., as_prob = as_prob))
      
      # If awfully undefined
      if (!as_prob && !is.finite(ans))
        return(-.Machine$double.xmax * 1e-100)
      else 
        return(ans)
      
    } else {
      
      exact_loglik(ergmito_ptr, params = params, ..., as_prob = as_prob)
      
    }
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
      nnets         = nnets(LHS) - length(excluded),
      used_attrs    = model_analysis$all_attrs,
      term_fun      = model_analysis$term_names,
      term_names    = colnames(target_stats),
      term_attrs    = model_analysis$term_attrs,
      excluded      = excluded
    ),
    class = "ergmito_loglik"
  )
  
  
  
}

#' @export
print.ergmito_loglik <- function(x, ...) {
  
  cat("ergmito log-likelihood function\n")
  cat("Number of networks: ", x$nnets, "\n")
  cat("Model: ", deparse(x$model_final), "\n")
  cat("Available elements by using the $ operator:\n")
  cat(sprintf("loglik: %s", deparse(x$loglik)[1]))
  cat(sprintf("grad  : %s", deparse(x$grad)[1]))
  
  invisible(x)
}

