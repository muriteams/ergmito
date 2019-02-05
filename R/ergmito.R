
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
#' @param stats List as returned by [ergm::ergm.allstats]. When this is provided,
#' the function does not call `ergm.allstats`, which can be useful in simulations.
#' @param ... Further arguments passed to the method. In the case of `ergmito`,
#' `...` are passed to [ergmito_formulae].
#' 
#' @seealso The function [plot.ergmito] for post-estimation diagnostics.
#' 
#' @return An list of class `ergmito`:
#' 
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
#' @examples 
#' 
#' # Generating a small graph
#' set.seed(12)
#' n <- 4
#' net <- sna::rgraph(n, tprob = .7)
#' 
#' model <- net ~ edges + mutual + balance
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
#' @importFrom stats optim terms rnorm
#' @importFrom MASS ginv
ergmito <- function(
  model,
  gattr_model = NULL,
  stats       = NULL,
  optim.args  = list(),
  ...
  ) {
  

  # Generating the objective function
  ergmitoenv <- environment(model)
  formulae   <- ergmito_formulae(
    model, gattr_model = gattr_model, stats = stats, env = ergmitoenv, ...)

  npars  <- formulae$npars

  if (!length(optim.args$control))
    optim.args$control <- list()
  
  optim.args$control$fnscale <- -1
  
  if (!length(optim.args$method)) optim.args$method <- "L-BFGS-B"
  if (optim.args$method == "L-BFGS-B") {
    if (!length(optim.args$lower)) optim.args$lower <- -100
    if (!length(optim.args$upper)) optim.args$upper <-  100
    if (!length(optim.args$control$factr))
      optim.args$control$factr <- 1e2
  }

  # Setting arguments for optim
  optim.args$par     <- (init <- stats::rnorm(npars, sd = .1))
  optim.args$fn      <- formulae$loglik
  optim.args$stats   <- stats$statmat
  optim.args$hessian <- TRUE
  
  ans <- do.call(stats::optim, optim.args)

  # Capturing the names of the parameters
  pnames         <- colnames(formulae$stats)
  names(ans$par) <- pnames
  covar.         <- -MASS::ginv(ans$hessian)
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
      network    = eval(model[[2]], envir = ergmitoenv)
    ),
    class="ergmito"
    )
  
  ans$nobs <- nvertex(ans$network)
  ans$nobs <- sum(ans$nobs*(ans$nobs - 1))
  
  ans
  
}

#' @export
#' @rdname ergmito
print.ergmito <- function(x, ...) {
  
  cat("\nERGMito estimates\n")
  print(structure(unclass(x), class="ergm"))
  invisible(x)
  
}

#' @export
#' @rdname ergmito
summary.ergmito <- function(object, ...) {
  cat("\nERGMito estimates\n")
  cat("\nformula: ", deparse(object$formulae$model), "\n\n")
  # Computing values
  sdval <- sqrt(diag(vcov(object)))
  z     <- coef(object)/sdval
  
  # Generating table
  data.frame(
    Estimate     = coef(object),
    `Std. Error` = sdval,
    `z value`    = z,
    `Pr(>|z|)`   = 2*pnorm(-abs(z)),
    row.names    = names(coef(object)),
    check.names  = FALSE
  )
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

