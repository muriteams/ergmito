
#' Estimation of ERGMs using Exact likelihood functions
#' 
#' As a difference from [ergm::ergm][ergm], `lergm` uses the exact log-likelihood
#' function for fitting the model. This implies that all the `2^(n*(n-1))` 
#' graphs are generated for computing the normalizing constant of the ERGM
#' model. This implies that models with up to 5 nodes are relatively simple
#' to fit, since more than that can become infeasible.
#' 
#' @param x,object An object of class `lergm`
#' @param model Model to estimate. See [ergm::ergm]. The only difference with
#' `ergm` is that the LHS can be a list of networks.
#' @param control List. Passed to [stats::optim].
#' @param stats List as returned by [ergm::ergm.allstats]. When this is provided,
#' the function does not call `ergm.allstats`, which can be useful in simulations.
#' @param ... Further arguments passed to the method. In the case of `lergm`,
#' `...` are passed to [lergm_formulae].
#' 
#' @seealso The function [plot.lergm] for post-estimation diagnostics.
#' 
#' @return An list of class `lergm`:
#' 
#' - `coef`          Named vector. Parameter estimates.
#' - `iterations`    Integer. Number of times the loglikelihood was evaluated
#'   (see [stats::optim]).
#' - `loglikelihood` Numeric. Final value of the objective function.
#' - `covar`         Square matrix of size `length(coef)`. Variance-covariance matrix
#' - `coef.init`     Named vector of length `length(coef)`. Initial set of parameters
#'   used in the optimization.
#' - `formulae`      An object of class [lergm_loglik][lergm_formulae].
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
#' ans_lergm <- lergm(model)
#' ans_ergm  <- ergm(model)
#' 
#' # The lergm should have a larger value
#' ergm.exact(ans_lergm$coef, model)
#' ergm.exact(ans_ergm$coef, model)
#' 
#' summary(ans_lergm)
#' summary(ans_ergm)
#' @importFrom stats optim terms rnorm
#' @importFrom MASS ginv
lergm <- function(
  model,
  control = list(maxit = 1e3, reltol=1e-100),
  stats   = NULL,
  ...
  ) {
  
  # Generating the objective function
  lergmenv <- environment(model)
  formulae <- lergm_formulae(model, stats = stats, env = lergmenv, ...)

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
      network       = eval(model[[2]], envir = lergmenv)
    ),
    class="lergm"
    )
  
}

#' @export
#' @rdname lergm
print.lergm <- function(x, ...) {
  
  cat("\nLittle ERGM estimates\n")
  print(structure(unclass(x), class="ergm"))
  invisible(x)
  
}

#' @export
#' @rdname lergm
summary.lergm <- function(object, ...) {
  cat("\nLittle ERGM estimates\n")
  summary(unclass(object))
}

# Methods ----------------------------------------------------------------------
#' @export
#' @rdname lergm
#' @importFrom stats coef logLik vcov
coef.lergm <- function(object, ...) {
  
  object$coef
  
}

#' @export
#' @rdname lergm
logLik.lergm <- function(object, ...) {
  
  structure(object$loglikelihood, class = "logLik", df = length(coef(object)))
  
}

#' @export
#' @rdname lergm
vcov.lergm <- function(object, ...) {
  
  object$covar
  
}
