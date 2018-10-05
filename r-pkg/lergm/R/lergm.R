#' @importFrom Rcpp sourceCpp
#' @useDynLib lergm, .registration = TRUE
NULL

#' Estimation of ERGMs using Exact likelihood functions
#' 
#' As a difference from [ergm::ergm][ergm], `lergm` uses the exact log-likelihood
#' function for fitting the model. This implies that all the `2^(n*(n-1))` 
#' graphs are generated for computing the normalizing constant of the ERGM
#' model. This implies that models with up to 5 nodes are relatively simple
#' to fit, since more than that can become infeasible.
#' 
#' @param formula,control,offset See [ergm::ergm].
#' @param allstats_force Logical. Passed to [ergm::ergm.allstats].
#' @param astats List as returned by [ergm::ergm.allstats]. When this is provided,
#' the function does not call `ergm.allstats`, which can be useful in simulations.
#' @return An object of class `ergm` and `lergm`.
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
lergm <- function(
  formula,
  control        = list(maxit = 1e3, reltol=1e-30),
  allstats_force = TRUE,
  offset         = NULL,
  stats          = NULL
  ) {
  
  # Calculate centered stats
  if (!length(stats))
    stats <- ergm::ergm.allstats(formula, zeroobs = TRUE, force = allstats_force)
  
  npars  <- ncol(stats$statmat)
  
  if (!length(offset))
    offset <- rep(FALSE, npars)
  
  control$fnscale <- -1
  
  ans <- stats::optim(
    par     = (init <- rnorm(npars)),
    fn      = exact_loglik,
    gr      = exact_loglik_gr,
    method  = "BFGS",
    weights = stats$weights,
    stats   = stats$statmat,
    control = control,
    hessian = TRUE
  )
  
  colnames(stats$statmat) <- as.character(
    statnet.common::list_rhs.formula(formula))
  names(ans$par) <- colnames(stats$statmat)
  
  lergm_class(
    coef          = ans$par,
    sample        = stats$statmat,
    sample.obs    = NULL,
    iterations    = ans$counts["function"],
    loglikelihood = ans$value,
    gradient      = NULL,
    covar         = -MASS::ginv(ans$hessian),
    network       = ergm::ergm.getnetwork(formula),
    coef.init     = init,
    formula       = formula,
    estimate      = "MLE",
    offset        = offset
    )
  
}

#' @export
print.lergm <- function(x, ...) {
  
  cat("\nLittle ERGM estimates\n")
  print(structure(unclass(x), class="ergm"))
  invisible(x)
  
}

#' @export
summary.lergm <- function(object, ...) {
  cat("\nLittle ERGM estimates\n")
  summary(structure(unclass(object), class="ergm"))
}

lergm_class <- function(
  coef,
  sample,
  sample.obs,
  iterations,
  loglikelihood,
  gradient,
  covar,
  network,
  coef.init,
  formula,
  estimate,
  mle.lik = loglikelihood,
  newnetwork = network,
  offset     = offset
) {
  
  n <- network$gal$n
  n <- n*(n-1)
  
  structure(
    list(
      coef = coef,
      sample = sample,
      sample.obs = sample.obs,
      iterations = iterations,
      loglikelihood = loglikelihood,
      gradient = gradient,
      covar = covar,
      network = network,
      coef.init = coef.init,
      formula = formula,
      estimate = estimate,
      control = ergm::control.ergm(),
      mle.lik = structure(loglikelihood, nobs = n, class="logLik", df=n*(n-1)/2-length(coef)),
      null.lik = NULL,
      constrained = NULL,
      constraints = ~.,
      newnetwork = newnetwork,
      est.cov = covar,
      target.stat = summary(formula),
      offset      = offset
    ),
    class=c("lergm", "ergm")
  )
  
}

# Optimization
exact_loglik <- function(params, weights, stats) {
  
  - log(weights %*% exp(stats %*% params))
  
}

exact_loglik_gr <- function(params, weights, stats) {
  
  exp_sum <- exp(stats %*% params)
  
  - 1/log(weights %*% exp_sum)[1]*(t(stats) %*% (exp_sum*weights))
  
  
}
# 
# model <- net ~ edges + mutual + balance
# 
# ans_lergm <- lergm(model)
# ans_ergm  <- ergm(model, control=control.ergm(MCMC.effectiveSize = 3000))
# 
# ergm.exact(ans_lergm$coef, model)
# # ergm.exact(ans_ergm$coef, model)
# 
# summary(ans_lergm)
# summary(ans_ergm)

# # Goodness of fit
# op <- par(mfrow=c(2, 3))
# plot(
#   gof(
#     ans_ergm, 
#     coef = ans_lergm$coef, 
#     control = control.gof.ergm(nsim=1e3)
#   ),
#   main="LERGM")
# plot.new()
# plot(
#   gof(
#     ans_ergm, 
#     coef = ans_ergm$coef, 
#     control = control.gof.ergm(nsim=1e3)
#   ),
#   main="ERGM")
# par(op)