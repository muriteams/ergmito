#' Estimation of ERGMs using Maximum Likelihood Estimation (MLE)
#' 
#' As a difference from [ergm::ergm][ergm], `ergmito` uses the exact log-likelihood
#' function for fitting the model. This implies that all the `2^(n*(n-1))` 
#' graphs are generated for computing the normalizing constant of the ERGM
#' model. As a rule of thumb, directed graphs with more than 5 vertices
#' should not be fitted using MLE, but instead MC-MLE as implemented in the
#' ergm package. The same applies for un-directed graphs with more than 8
#' vertices..
#' 
#' @param x,object An object of class `ergmito`
#' @param model Model to estimate. See [ergm::ergm]. The only difference with
#' `ergm` is that the LHS can be a list of networks.
#' @param gattr_model A formula. Model specification for graph attributes. This
#' is useful when using multiple networks.
#' @param optim.args List. Passed to [stats::optim].
#' @param target.stats A matrix of target statistics (see [ergm::ergm]).
#' @template stats
#' @param init A numeric vector. Sets the starting parameters for the
#' optimization routine. Default is a vector of zeros.
#' @param use.grad Logical. When `TRUE` passes the gradient function to `optim`.
#' This is intended for testing only (internal use).
#' @param ntries Integer scalar. Number of tries to estimate the MLE (see details).
#' @param keep.stats Logical scalar. When `TRUE` (the default), the matrices
#' and vectors associated with the sufficient statistics will be returned.
#' Otherwise the function discards them. This may be useful for saving memory
#' space when estimating multiple models.
#' @param ... Further arguments passed to the method. In the case of `ergmito`,
#' `...` are passed to [ergmito_formulae].
#' 
#' @seealso The function [plot.ergmito] for post-estimation diagnostics.
#' 
#' @return An list of class `ergmito`:
#' 
#' - `call`          The program call.
#' - `coef`          Named vector. Parameter estimates.
#' - `iterations`    Integer. Number of times the log-likelihood was evaluated
#'   (see [stats::optim]).
#' - `mle.lik`       Numeric. Final value of the objective function.
#' - `null.lik`      Numeric. Final value of the objective function for the null model.
#' - `covar`         Square matrix of size `length(coef)`. Variance-covariance matrix
#'   computed using the exact hessian as implemented in [exact_hessian].
#' - `coef.init`     Named vector of length `length(coef)`. Initial set of parameters
#'   used in the optimization.
#' - `formulae`      An object of class [ergmito_loglik][ergmito_formulae].
#' - `nobs`          Integer scalar. Number of networks in the model.
#' - `network`       Networks passed via `model`.
#' - `optim.out`,`optim.args` Results from the optim call and arguments passed to it.
#' - `status`,`note` Convergence code. See [check_convergence]
#' - `best_try`      Integer scalar. Index of the run with the highest log-likelihood value.
#' - `history`       Matrix of size `ntries * (k + 1)`. History of the parameter
#'   estimates and the reached log-likelihood values.
#' - `timer`         Vector of times (for benchmarking). Each unit marks the starting
#'   point of the step.
#'   
#' @section MLE:
#' 
#' Maximum Likelihood Estimates are obtained using the [stats::optim] function.
#' The default method for maximization is `BFGS` using both the log-likelihood
#' function and its corresponding gradient.
#'  
#' Another important factor to consider is the existence of the MLE estimates
#' As shown in Handcock (2003), if the observed statistics are near the border
#' if the support function (e.g. too many edges or almost none), then, even if
#' the MLE estimates exists, the optimization function may not be able to reach
#' the optima. Moreover, if the target (observed) statistics live in the boundary,
#' then the MLE estimates do not exists. In general, this should not be an issue
#' in the context of the pooled model, as the variability of observed statistics
#' should be enough to avoid those situations.
#' 
#' The function `ergmito` will try to identify possible cases of non-existence,
#' of the MLE, and if identified then try to re estimate the model parameters using
#' larger values than the ones obtained, if the log-likelihood is greater, then it is 
#' assumed that the model is degenerate and the corresponding values will be
#' replaced with either `+Inf` or  `-Inf`. By default, this behavior is checked
#' anytime that the absolute value of the estimates is greater than 5, or the
#' sufficient statistics were flagged as potentially outside of the interior of
#' the support (close to zero or to its max).
#' 
#' In the case of `ntries`, the optimization is repeated that number of times,
#' each time perturbing the `init` parameter by adding a Normally distributed
#' vector. The result which reaches the highest log-likelihood will be the one
#' reported as parameter estimates. This feature is intended for testing only.
#' Anecdotally, `optim` reaches the max in the first try.
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
#' # Example 3: Likelihood ratio test using the lmtest R package
#' 
#' if (require(lmtest)) {
#'   data(fivenets)
#'   model1 <- ergmito(fivenets ~ edges + nodematch("female"))
#'   model2 <- ergmito(fivenets ~ edges + nodematch("female") + mutual)
#'   
#'   lrtest(model1, model2)
#'   # Likelihood ratio test
#'   # 
#'   # Model 1: fivenets ~ edges + nodematch("female") 
#'   # Model 2: fivenets ~ edges + nodematch("female") + mutual
#'   #   #Df  LogLik Df  Chisq Pr(>Chisq)
#'   # 1   2 -34.671                     
#'   # 2   3 -34.205 1 0.9312     0.3346
#' }
#' 
#' 
#' @importFrom stats optim terms rnorm
#' @importFrom MASS ginv
#' @name ergmito
NULL

ERGMITO_DEFAULT_OPTIM_CONTROL <- list(
  reltol = .Machine$double.eps ^ 3/4
)

#' @export
#' @rdname ergmito
ergmito <- function(
  model,
  gattr_model   = NULL,
  stats.weights = NULL,
  stats.statmat = NULL,
  optim.args    = list(),
  init          = NULL,
  use.grad      = TRUE,
  target.stats  = NULL,
  ntries        = 1L,
  keep.stats    = TRUE,
  ...
  ) {

  # Keeping track of time
  timer_start <- Sys.time()

  # Generating the objective function
  ergmitoenv <- environment(model)

  formulae   <- ergmito_formulae(
    model,
    gattr_model   = gattr_model, 
    target.stats  = target.stats,
    stats.weights = stats.weights,
    stats.statmat = stats.statmat,
    env           = ergmitoenv,
    ...
    )
  timer <- c(ergmito_formulae = difftime(Sys.time(), timer_start, units = "secs"))
  
  # Verifying existence of MLE
  timer0 <- Sys.time()
  support <- check_support(
    formulae$target.stats,
    formulae$stats.statmat
    )
  timer <- c(timer, check_support = difftime(Sys.time(), timer0, units = "secs"))
  
  npars  <- formulae$npars
  
  # Checking the values of the initial parameters, if an undefined value is passed
  # then replace it with a very large but tend to infinite value
  if (!length(init)) 
    init <- rep(0, npars)

  # Checking optim parameters --------------------------------------------------
  if (!length(optim.args$control))
    optim.args$control <- list()
  optim.args$control$fnscale <- -1
  
  for (n in names(ERGMITO_DEFAULT_OPTIM_CONTROL))
    if (!length(optim.args$control[[n]]))
      optim.args$control[[n]] <- ERGMITO_DEFAULT_OPTIM_CONTROL[[n]]
  
  # For BFGS 
  if (!length(optim.args$method)) 
    optim.args$method <- "BFGS"
    
  # Passed (and default) other than the functions
  optim.args0 <- optim.args
  
  # Setting arguments for optim
  optim.args$fn <- formulae$loglik
  if (use.grad) 
    optim.args$gr <- formulae$grad
  optim.args$hessian       <- FALSE
  optim.args$par           <- init
  
  # Will try to solve the problem more than once... if needed
  ntry <- 1L
  history <- matrix(
    NA, nrow = ntries, ncol = formulae$npars + 1,
    dimnames = list(1L:ntries, c(formulae$term.names, "value"))
    )
  
  timer0 <- Sys.time()
  while (ntry <= ntries) {
    
    # Maximizign the likelihood and storing the value
    cur_ans <- do.call(stats::optim, optim.args)
    
    if ((ntry == 1L) || (ans$value < cur_ans$value)) {
      ans      <- cur_ans
      best_try <- ntry
    }
     
    
    # Storing the current value
    history[ntry, ] <- c(cur_ans$par, cur_ans$value)
    
    # We don't need to do this again.
    if (ntries == 1L)
      break
    
    # Resetting the parameters for the optimization, now this time we start
    # from the init parameters + some random value
    optim.args$par <- stats::rnorm(formulae$npars, -2, 2)
    
    ntry <- ntry + 1
    
  }
  timer <- c(timer, optim = difftime(Sys.time(), timer0, units = "secs"))
  
  # Checking the convergence
  timer0 <- Sys.time()
  estimates <- check_convergence(
    optim_output = ans,
    model        = formulae,
    support      = support
    )
  timer <- c(
    timer,
    chec_covergence = difftime(Sys.time(), timer0, units = "secs")
    )
  
  # Capturing model
  if (!inherits(model, "formula"))
    model <- eval(model)
  
  # Null loglik
  ll0 <-formulae$loglik(params = rep(0, length(estimates$par)))
  
  ans <- structure(
    list(
      call       = match.call(),
      coef       = estimates$par,
      iterations = ans$counts["function"],
      mle.lik    = structure(estimates$ll, class="logLik", df=length(ans$par)),
      null.lik   = structure(ll0, class="logLik", df=0),
      covar      = estimates$vcov,
      coef.init  = init,
      formulae   = formulae,
      nobs       = NA,
      network    = eval(model[[2]], envir = ergmitoenv),
      optim.out  = ans,
      optim.args = optim.args0,
      status     = estimates$status,
      note       = estimates$note,
      best_try   = best_try,
      history    = history
    ),
    class = c("ergmito")
    )
  
  if (!keep.stats) {
    ans$formulae$stats.weights <- NULL
    ans$formulae$stats.statmat <- NULL
  }
  
  ans$nobs <- nvertex(ans$network)
  ans$nobs <- sum(ans$nobs*(ans$nobs - 1))
  
  timer <- c(timer, total = difftime(Sys.time(), timer_start, units = "secs"))
  ans$timer <- timer
  ans
  
}

#' @export
#' @rdname ergmito
print.ergmito <- function(x, ...) {
  
  cat("\nERGMito estimates\n")
  if (length(x$note))
    cat(sprintf("note: %s\n", x$note))
    
  print(structure(unclass(x), class="ergm"))
  invisible(x)
  
}

#' @export
#' @rdname ergmito
summary.ergmito <- function(object, ...) {

  # Computing values
  sdval <- sqrt(diag(vcov(object)))
  z     <- coef(object)/sdval
  
  is_boot <- inherits(object, "ergmito_boot")
  
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
    note        = object$note,
    R           = ifelse(is_boot, object$R, 1L)
    ),
    class = c("ergmito_summary", if (is_boot) "ergmito_summary_boot" else  NULL)
  )
  
  ans
}

#' @export
#' @rdname ergmito
print.ergmito_summary <- function(
  x,
  ...
  ) {

  cat("\nERGMito estimates\n")
  
  if (x$R > 1L)
    cat("\n(bootstrapped model with ", x$R, " replicates.)\n")
  
  if (length(x$note))
    cat(sprintf("note: %s\n", x$note))
  cat("\nformula: ", x$model, "\n\n")
  stats::printCoefmat(
    x$coefs,
    signif.stars  = TRUE,
    signif.legend = TRUE
    )
  
  cat(paste("AIC:", format(x$aic), 
            "  ", "BIC:", format(x$bic), 
            "  ", "(Smaller is better.)", "\n", sep = " "))
  
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
#' @param solver Function. Used to compute the inverse of the hessian matrix. When
#' not null, the variance-covariance matrix is recomputed using that function.
#' By default, `ergmito` uses [MASS::ginv].
#' @rdname ergmito
vcov.ergmito <- function(object, solver = NULL, ...) {
  
  if (is.null(solver))
    return(object$covar)
  
  structure(
    - solver(object$optim.out$hessian),
    dimnames = dimnames(object$covar)
  )
  
}

#' @export
#' @rdname ergmito
formula.ergmito <- function(x, ...) {
  
  x$formulae$model
  
}
