
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
#' @param control List. Passed to [stats::optim].
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
  control = list(maxit = 1e3, reltol=1e-100),
  stats   = NULL,
  ...
  ) {
  
  # Generating the objective function
  ergmitoenv <- environment(model)
  formulae <- ergmito_formulae(model, stats = stats, env = ergmitoenv, ...)

  npars  <- formulae$npars
  
  control$fnscale <- -1
  
  ans <- stats::optim(
    par     = (init <- stats::rnorm(npars)),
    method  = "BFGS",
    fn      = formulae$loglik,
    # gr      = formulae$grad,
    stats   = stats$statmat,
    control = control,
    hessian = TRUE
  )
  
  # Capturing the names of the parameters
  pnames         <- attr(stats::terms(formulae$model), "term.labels")
  names(ans$par) <- pnames
  covar.         <- -MASS::ginv(ans$hessian)
  dimnames(covar.) <- list(pnames, pnames)
  
  # Capturing model
  if (!inherits(model, "formula"))
    model <- eval(model)
  
  structure(
    list(
      call          = match.call(),
      coef          = ans$par,
      iterations    = ans$counts["function"],
      loglikelihood = ans$value,
      covar         = covar.,
      coef.init     = init,
      formulae      = formulae,
      network       = eval(model[[2]], envir = ergmitoenv)
    ),
    class="ergmito"
    )
  
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
  summary(unclass(object))
}

# Methods ----------------------------------------------------------------------
#' @export
#' @rdname ergmito
#' @importFrom stats coef logLik vcov
coef.ergmito <- function(object, ...) {
  
  object$coef
  
}

#' @export
#' @rdname ergmito
logLik.ergmito <- function(object, ...) {
  
  structure(object$loglikelihood, class = "logLik", df = length(coef(object)))
  
}

#' @export
#' @rdname ergmito
vcov.ergmito <- function(object, ...) {
  
  object$covar
  
}
