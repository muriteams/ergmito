
check_degeneracy <- function(target.stats, stats, threshold = .8, warn = TRUE) {
  
  # Retrieving the matrix of sufficient statistics
  statsmats <- lapply(stats, "[[", "statmat")
  
  res <- structure(
    vector("logical", ncol(target.stats)),
    names = colnames(target.stats)
    )
  
  for (k in 1L:ncol(target.stats)) {
    
    # Retrieving the space range
    stats_range <- stats
    
    # Looking for degeneracy at the k-th parameter
    stat_range <- lapply(statsmats, "[", i=, drop = TRUE)
    stat_range <- lapply(stat_range, range)
    stat_range <- do.call(rbind, stat_range)
    
    res[k] <- mean((target.stats[, k] == stat_range[, 1L]) | 
      (target.stats[, k] == stat_range[, 2L]))
    
  }
  
  attr(res, "threshold") <- threshold
  test <- which(res >= threshold)
  if (length(test)) {
    
    if (warn)
    warning("The observed statistics (target.statistics) are near or at the",
            "boundary of its support, i.e. the Maximum Likelihood Estimates may",
            "not exist or be hard to be estimated. In particular,", 
            " the statistics \"", paste(names(res)[test], collapse="\", \""), 
            "\".", call. = FALSE, immediate. = TRUE)
    
    attr(res, "degenerate") <- TRUE
    attr(res, "which")      <- test
  } else {
    attr(res, "degenerate") <- FALSE
    attr(res, "which")      <- NULL
  }
  
  res
  
  
}

#' Estimation of ERGMs using Exact likelihood functions
#' 
#' As a difference from [ergm::ergm][ergm], `ergmito` uses the exact log-likelihood
#' function for fitting the model. This implies that all the `2^(n*(n-1))` 
#' graphs are generated for computing the normalizing constant of the ERGM
#' model. This implies that models with up to 5 nodes are relatively simple
#' to fit, since more than that can become infeasible.
#' 
#' @param x,object An object of class `ergmito`
#' @param model Model to estimate. See [ergm::ergm]. The only difference with
#' `ergm` is that the LHS can be a list of networks.
#' @param gattr_model A formula. Model especification for graph attributes. This
#' is useful when using multiple networks.
#' @param optim.args List. Passed to [stats::optim].
#' @param target.stats A matrix of target statistics (see [ergm::ergm]).
#' @param stats List as returned by [ergm::ergm.allstats]. When this is provided,
#' the function does not call `ergm.allstats`, which can be useful in simulations.
#' @param init Either a numeric vector, or a matrix with `ntries` rows. Sets
#' the starting parameters for the optimization routine.
#' @param use.grad Logical. When `TRUE` passes the gradient function to `optim`.
#' This is intended for testing only (internal use).
#' @param ntries Integer. In case of no convergence, number of times that the 
#' `optim` function should be re-run.
#' @param ... Further arguments passed to the method. In the case of `ergmito`,
#' `...` are passed to [ergmito_formulae].
#' 
#' @seealso The function [plot.ergmito] for post-estimation diagnostics.
#' 
#' @return An list of class `ergmito`:
#' 
#' - `call`          The program call.
#' - `coef`          Named vector. Parameter estimates.
#' - `iterations`    Integer. Number of times the loglikelihood was evaluated
#'   (see [stats::optim]).
#' - `loglikelihood` Numeric. Final value of the objective function.
#' - `covar`         Square matrix of size `length(coef)`. Variance-covariance matrix
#' - `coef.init`     Named vector of length `length(coef)`. Initial set of parameters
#'   used in the optimization.
#' - `formulae`      An object of class [ergmito_loglik][ergmito_formulae].
#' - `network`       Networks passed via `model`.
#' 
#' @export
#' @section MLE:
#' 
#' Maximum Likelihood Estimates are obtained using the [stats::optim] function.
#' The default method for maximization is `BFGS` using both the loglikelihood
#' function and its corresponding gradient.
#' 
#' Anecdotically, while the optimization process should be realitvely straight
#' forward, in some cases a near-degenerate set of observed statistics may
#' yield suboptimal solutions. Because of this reason, the parameter `ntries`
#' allows the user to specify how many runs of `optim` should be done before
#' returning. In general, optimization via `optim` is very fast, which is why
#' this part of the function doesn't impact the actual wall time significantly
#' (the most part of the time that the function takes to run is caused by the
#' computation of the stats natrices using the [ergm::ergm.allstats] function.)
#' 
#' Another important factor to consider is the existance of the MLE estiates.
#' As shown in Handcock (2003), if the observed statitcs are near the border
#' if the support function (e.g. too many edges or almost none), then, even if
#' the MLE estimates exists, the optimization function may not be able to reach
#' the optima. Moreover, if the target (observed) statistics live in the boundary,
#' then the MLE estimates do not exists. In general, this should not be an issue
#' in the context of the pooled model, as the variability of observed statistics
#' should be enough to avoid those situations.
#' 
#' The function `ergmito` will try to identify possible cases of model degeneracy,
#' and if identified, then try to re estimate the model parameters using larger
#' values than the ones obtained, if the log-likelihood is greater, then it is 
#' assumed that the model is degenerate and the corresponding values will be
#' replaced with either `+Inf` or  `-Inf`.
#' 
#' 
#' @examples 
#' 
#' # Generating a small graph
#' set.seed(12)
#' n <- 4
#' net <- rbernoulli(n, p = .7)
#' 
#' model <- net ~ edges + mutual
#' 
#' library(ergm)
#' ans_ergmito <- ergmito(model)
#' ans_ergm  <- ergm(model)
#' 
#' # The ergmito should have a larger value
#' ergm.exact(ans_ergmito$coef, model)
#' ergm.exact(ans_ergm$coef, model)
#' 
#' summary(ans_ergmito)
#' summary(ans_ergm)
#' 
#' # Example 2: Estimating an ERGMito using data with know DGP parameters -----
#' data(fivenets) 
#' 
#' model1 <- ergmito(fivenets ~ edges + nodematch("female"))
#' summary(model1) # This data has know parameters equal to -2.0 and 2.0
#' 
#' 
#' @importFrom stats optim terms rnorm
#' @importFrom MASS ginv
ergmito <- function(
  model,
  gattr_model  = NULL,
  stats        = NULL,
  optim.args   = list(),
  init         = NULL,
  ntries       = 5L,
  use.grad     = TRUE,
  target.stats = NULL,
  ...
  ) {
  
  # Generating the objective function
  ergmitoenv <- environment(model)
  formulae   <- ergmito_formulae(
    model,
    gattr_model  = gattr_model, 
    target.stats = target.stats,
    stats        = stats,
    env          = ergmitoenv,
    ...
    )

  # Verifying existance of MLE
  degeneracy <- check_degeneracy(formulae$target.stats, formulae$stats)
  
  npars  <- formulae$npars
  
  # Checking initial parameters
  if (!length(init))
    init <- matrix(stats::runif(npars * ntries, -1.0, 1.0), ncol = npars)
  else if (length(init) == npars && ntries != 1) {
    warning(
      "The set of initial parameters will be extended to match `ntries`.",
      " The extended initial points will be drawn from a N(0, 1).",
      call. =  FALSE)
    
    init <- rbind(init, matrix(stats::runif(npars * (ntries - 1L), -1.0, 1.0), ncol = npars))
  } else
    stop(
      "Invalid number of inital parameters (`init`).",
      " See the section 'MLE' in ?ergmito.",
      call. = FALSE
      )
  

  # Checking optim parameters --------------------------------------------------
  if (!length(optim.args$control))
    optim.args$control <- list()
  optim.args$control$fnscale <- -1
  
  # For BFGS 
  if (!length(optim.args$method)) 
    optim.args$method <- "BFGS"
    
  if (optim.args$method == "BFGS" && !length(optim.args$control$reltol)) {
      optim.args$control$reltol <- .Machine$double.eps*2
  }
  
  # For L-BFGS-B
  if (optim.args$method == "L-BFGS-B") {
    if (!length(optim.args$lower)) optim.args$lower <- -50
    if (!length(optim.args$upper)) optim.args$upper <-  50
    if (!length(optim.args$control$factr))
      optim.args$control$factr <- 1e2
  }

  # Setting arguments for optim
  optim.args$fn <- formulae$loglik
  if (use.grad) 
    optim.args$gr <- formulae$grad
  optim.args$stats <- stats$statmat
  optim.args$hessian <- TRUE
  
  ans <- vector("list", ntries)
  for (i in 1:ntries) {
    optim.args$par <- init[i, , drop=TRUE]
    ans <- do.call(stats::optim, optim.args)
    
    if (ans$convergence == 0)
      break
  }
  
  if ((i == ntries) && !ans$convergence)
    warning("The optim function did not converged.", call. = FALSE)

  # If denegeracy is plausible, then the solution may be close to infinite
  if (attr(degeneracy, "degenerate")) {
    
    optim.args$par <- ans$par
    optim.args$par[attr(degeneracy, "which")] <-
      optim.args$par[attr(degeneracy, "which")]*2.0
    
    # Rerunning the optimization funtion with a parameter near infinite
    # to see if the likelihood function increases value. If it does, then
    # we will update the values
    ans1 <- do.call(stats::optim, optim.args)
    
    if (ans$value < ans1$value) {
      
      # This should come with a warning
      warning("A correction has been made in order to achieve the right MLE.",
              " In this case, the theoretical estimates are near +-Inf, hence",
              " the MLE may not exist.", call. = FALSE, immediate. = TRUE)
      
      ans <- ans1
      ans$par[attr(degeneracy, "which")] <- 
        Inf * sign(ans$par[attr(degeneracy, "which")])
    }
    
  }
  
  # Capturing the names of the parameters
  pnames         <- colnames(formulae$target.stats)
  names(ans$par) <- pnames
  covar.         <- -MASS::ginv(ans$hessian)
  
  # If we reached this level, then we shouldn't be reporting any meaningful
  # for that parameter
  if (any(is.infinite(ans$par))) {
    covar.[attr(degeneracy, "which"), ] <- NA
    covar.[, attr(degeneracy, "which")] <- NA
  }
  
  dimnames(covar.) <- list(pnames, pnames)
  
  # Capturing model
  if (!inherits(model, "formula"))
    model <- eval(model)
  
  # Null loglik
  ll0 <- formulae$loglik(rep(0, length(pnames)))
  
  ans <- structure(
    list(
      call       = match.call(),
      coef       = ans$par,
      iterations = ans$counts["function"],
      mle.lik    = structure(ans$value, class="logLik", df=length(ans$par)),
      null.lik   = structure(ll0, class="logLik", df=0),
      covar      = covar.,
      coef.init  = init,
      formulae   = formulae,
      nobs       = NA,
      network    = eval(model[[2]], envir = ergmitoenv),
      init       = init,
      optim.out  = ans,
      degeneracy = degeneracy
    ),
    class = c("ergmito")
    )
  
  ans$nobs <- nvertex(ans$network)
  ans$nobs <- sum(ans$nobs*(ans$nobs - 1))
  
  ans
  
}

#' @export
#' @rdname ergmito
print.ergmito <- function(x, ...) {
  
  cat("\nERGMito estimates\n")
  if (length(x$degeneracy) && attr(x$degeneracy, "degenerate"))
    cat("Note: Degenerate or near-degenerate model. The MLE may not exists.\n")
  if (x$optim.out$convergence != 0)
    cat("Note: The optimzation did not converged.\n")
    
  print(structure(unclass(x), class="ergm"))
  invisible(x)
  
}

#' @export
#' @rdname ergmito
summary.ergmito <- function(object, ...) {

  # Computing values
  sdval <- sqrt(diag(vcov(object)))
  z     <- coef(object)/sdval
  
  # Generating table
  ans <- structure(
    list(
      coefs = data.frame(
      Estimate     = coef(object),
      `Std. Error` = sdval,
      `z value`    = z,
      `Pr(>|z|)`   = 2*stats::pnorm(-abs(z)),
      row.names    = names(coef(object)),
      check.names  = FALSE
    ),
    aic         = stats::AIC(object),
    bic         = stats::BIC(object),
    model       = deparse(object$formulae$model),
    degeneracy  = object$degeneracy,
    convergence = object$optim.out$convergence
    ),
    class = "ergmito_summary"
  )
  
  ans
}

#' @export
#' @rdname ergmito
print.ergmito_summary <- function(x, ...) {

  cat("\nERGMito estimates\n")
  if (length(x$degeneracy) && attr(x$degeneracy, "degenerate"))
    cat("Note: Degenerate or near-degenerate model. The MLE may not exists\n")
  if (x$convergence != 0)
    cat("Note: The optimzation did not converged.\n")
  cat("\nformula: ", x$model, "\n\n")
  print(x$coefs)
  invisible(x)
    
}


# Methods ----------------------------------------------------------------------
#' @export
#' @rdname ergmito
#' @importFrom stats coef logLik vcov nobs
coef.ergmito <- function(object, ...) {
  
  object$coef
  
}

#' @export
#' @rdname ergmito
logLik.ergmito <- function(object, ...) {
  
  object$mle.lik
  
}

#' @export
#' @rdname ergmito
nobs.ergmito <- function(object, ...) {
  
  object$nobs
  
}

#' @export
#' @rdname ergmito
vcov.ergmito <- function(object, ...) {
  
  object$covar
  
}

