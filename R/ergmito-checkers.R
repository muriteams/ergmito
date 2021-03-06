#' @rdname check_convergence
#' @param target_stats,stats_statmat See [ergmito_formulae].
#' @param threshold Numeric scalar. Confidence range for flagging an observed
#' statistic as potentially near the boundary.
#' @param warn logical scalar.
check_support <- function(
  target_stats,
  stats_statmat,
  threshold = .8,
  warn      = TRUE
  ) {
  
  res <- structure(
    logical(ncol(target_stats)),
    names = colnames(target_stats)
  )
  
  # In the case of failing to converge, this tells what should be
  # the best guess for the sign of Inf.
  possible_sign <- double(ncol(target_stats))
  
  for (k in 1L:ncol(target_stats)) {
    
    # Looking for not in the interior at the k-th parameter
    stat_range <- lapply(stats_statmat, "[", i=, j = k, drop = TRUE)
    stat_range <- lapply(stat_range, range)
    stat_range <- do.call(rbind, stat_range)
    
    res[k] <- mean(
      (target_stats[, k] - stat_range[,1L])/
        (stat_range[, 2L] - stat_range[, 1L] + 1e-20)
      )
    
    # If on average is less than .5, then is negative, otherwise
    # is positive.
    possible_sign[k] <- ifelse(res[k] < .5, -1, 1)
    
  }
  
  attr(res, "threshold") <- threshold
  test <- which((res >= threshold) | (res <= threshold))
  if (length(test)) {
    
    if (warn)
      warning_ergmito("The observed statistics (target.statistics) are near or at the ",
              "boundary of its support, i.e. the Maximum Likelihood Estimates may",
              "not exist or be hard to be estimated. In particular,", 
              " the statistic(s) \"", paste(names(res)[test], collapse="\", \""), 
              "\".", call. = FALSE, immediate. = TRUE)
    
    attr(res, "interior") <- FALSE
    attr(res, "which")    <- test
    attr(res, "sign")     <- possible_sign
    
  } else {
    attr(res, "interior") <- TRUE
    attr(res, "which")    <- NULL
    attr(res, "sign")     <- possible_sign
  }
  
  res
  
  
}

#' This is used to generate notes
#' @noRd
CONVERGENCE_DICTIONARY <- list(
  `00` = NULL,
  `01` = "optim converged, but the Hessian is not p.s.d.",
  `10` = "optim did not converged, but the estimates look OK.",
  `11` = "optim did not converged, and the Hessian is not p.s.d.",
  `20` = "A subset of the parameters estimates was replaced with +/-Inf.",
  `21` = paste(
    "A subset of the parameters estimates was replaced with +/-Inf, ",
    "and the Hessian matrix is not p.s.d."
  ),
  `30` = "All parameters went to +/-Inf suggesting that the MLE may not exists.",
  `31` = "All parameters went to +/-Inf suggesting that the MLE may not exists. Also, the Hessian is not p.s.d."
)

map_convergence_message <- function(x) {
  CONVERGENCE_DICTIONARY[[sprintf("%02d", x)]]
}

#' Check the convergence of ergmito estimates
#' 
#' This is an internal function used to check the convergence of the optim function.
#' 
#' @param optim_output A list output from the [stats::optim] function.
#' @param model An object of class [ergmito_loglik].
#' @param support As returned by `check_support`.
#' @param crit Numeric scalar. Level at which a parameter estimate
#' will be questioned.
#' @return A list with the following components:
#' 
#' - `par` Updated set of parameters
#' - `vcov` Updated variance-covariance matrix
#' - `valid` Vector of integers with the parameters that are marked as OK.
#' - `status` Return code of the analysis. See details.
#' - `note` A note describing the status.
#' 
#' @section Return codes: 
#' 
#' The function makes an analysis of the outcome of the model and makes the corresponding
#' adjustments when required. In particular, we check:
#' 
#' 1. Whether the optimization algorithm converged or not
#' 
#' 2. If the obtained estimates maximize the function. If this is not the case,
#'    the function checks whether the MLE may not exist. This usually happens
#'    when the log-likelihood function can improve by making increments to parameters
#'    that are already tagged as large. If the ll improves, then the value is
#'    replaced with `Inf` (+- depending on the sign of the parameter).
#'    
#' 3. If the Hessian is semi-positive-definite, i.e. if it is invertible. If it 
#'    is not, it usually means that the function did not converged, in which 
#'    case we will use [MASS::ginv] instead.
#'    
#' The return codes are composed of two numbers, the first number gives information
#' regarding of the parameter estimates, while the second number give information
#' about the variance-covariance matrix.
#' 
#' Column 1:
#' 
#' - 0: Converged and estimates at the max.
#' - 1: It did not converged, but I see no issue in the max.
#' - 2: One or more estimates went to +/-Inf
#' - 3: All went to hell. All estimates went to +/-Inf
#' 
#' Column 2: 
#' 
#' - 0: Hessian is p.s.d.
#' - 1: Hessian is not not p.s.d.
#' 
#' Possible codes and corresponding messages:
#' 
#' - 00 All OK (no message).
#' - 01 \Sexpr{ergmito:::map_convergence_message(01)}. 
#' - 10 \Sexpr{ergmito:::map_convergence_message(10)}. 
#' - 11 \Sexpr{ergmito:::map_convergence_message(11)}. 
#' - 20 \Sexpr{ergmito:::map_convergence_message(20)}. 
#' - 21 \Sexpr{ergmito:::map_convergence_message(21)}. 
#' - 30 \Sexpr{ergmito:::map_convergence_message(30)}. 
#' 
#' @keywords Internal
check_convergence <- function(
  optim_output,
  model,
  support,
  crit = 5.0
  ) {
  
  # Baseline check... we are passing the right type.
  if (!inherits(model, "ergmito_loglik"))
    stop(
      "`model` has to be an object of class \"ergmito_loglik\". This is an ",
      "internal function, are you sure you want to use this directly?.",
      call. = TRUE
      )

  
  # Checking values of the convergence
  to_check <- c(
    which(abs(optim_output$par) > crit),
    attr(support, "which")
    )
  to_check <- sort(unique(to_check))
  
  k <- length(optim_output$par)
  estimates <- list(
    par    = structure(optim_output$par, names = model$term_names),
    vcov   = matrix(
      0.0, nrow = k, ncol = k,
      dimnames = with(model, list(term_names, term_names))
      ),
    valid  = 1L:k,
    status = ifelse(optim_output$convergence == 0L, 0L, 1L),
    note   = NULL,
    ll     = optim_output$value
    )
  
  # We will update this later
  # estimates$vcov[] <- optim_output$hessian
  estimates$vcov[] <- model$hess(optim_output$par)
  
  # Step 1: Checking parameter estimates ---------------------------------------
  if (length(to_check)) {
    
    # Should we replace with Inf?
    newpars  <- optim_output$par
    modified <- NULL 
    for (i in to_check) {
      
      # Near to 0/1 should always be set as infinite
      if ((support[i] <= 1e-10) | (support[i] >= (1 - 1e-10))) {
        
        newpars[i] <- attr(support, "sign")[i] *Inf
        modified   <- c(modified, i)
        
      }
      
    }
    
    # Updating parameters, if needed
    estimates$par[] <- newpars
    estimates$valid <- setdiff(estimates$valid, modified)
    
    if (length(modified)) {
      # Updating the hessian matrix. We cannot use infite values for this step
      # since optimHess will return with an error. That's why we just use a
      # very large value instead
      newpars <- estimates$par
      newpars[!is.finite(newpars)] <- sign(newpars[!is.finite(newpars)]) * 100
      
      estimates$vcov[] <- model$hess(newpars)
      estimates$vcov[!is.finite(estimates$par),] <- 0
      estimates$vcov[,!is.finite(estimates$par)] <- 0
      
      # The observed likelihood will change as well, it may be the case that it
      # becomes undefined b/c of the fact that 0 * Inf = NaN, yet the right
      # value is well defined
      estimates$ll <- model$loglik(newpars)
    }
    
    # Are we in hell?
    if (!length(estimates$valid)) {
      
      warning_ergmito(
        "All parameters went to +-Inf. This suggests the MLE may not exist.",
        call. = FALSE
      )
      
      estimates$status <- 30L
      
    } else if (length(modified)) {
      
      estimates$status <- 20L
      
    }
    
    
  }
  
  # Step 2: Checking variance cov-matrix ---------------------------------------
  
  # Trying to compute the variance co-variance matrix, on the right ones.
  vcov. <- tryCatch(solve(-estimates$vcov), error = function(e) e)
  
  if (inherits(vcov., "error")) {
    
    # Trying to estimate using the Generalized Inverse
    vcov. <- tryCatch(MASS::ginv(-estimates$vcov), error = function(e) e)
    if (inherits(vcov., "error"))
      estimates$vcov[] <- NA
    else
      estimates$vcov[] <- vcov.
    
    # Wasn't able to fully estimate it...
    estimates$status <- estimates$status + 1L
    
  } else {
    estimates$vcov[] <- vcov.
  }
  
  # Returning, asis
  estimates$note <- map_convergence_message(estimates$status)
  return(estimates)

  
}


#' Possible problems:
#' 
#' - The sufficient statistics lay too proximate to the boundaries of its 
#'   support. For example a triadic term has been included in a model in which
#'   the observed graph has either zero or the max number of possible triads.
#'   
#'   If this happens, it a lot of times it could be that the MLE does not exists.
#'  
#' - A multi-modal distribution: This makes the maximization problem hard in 
#'   terms of finding a global maximum.
#'   
#' Flags:
#' 
#' - Not semi-positive definite Hessian: The global has not been found. Most of
#'   the cases it could be related to a model with a set of sufficient statistics
#'   near the boundary
#'   
#' - Parameters estimate too big. Anything outside of 5 is a candidate for a 
#'   problem in the function. 
#'   
#' - Maximum number of iterations reached: This is related to the previous
#'   situation. If the sufficient statistics are near the boundary of the
#'   support, then in many cases the theoretical optimum lies at +/-Infinite.
#'   
#'   
#'   
#' 
#' @noRd 
