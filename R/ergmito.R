#' Estimation of ERGMs using Maximum Likelihood Estimation (MLE)
#' 
#' `ergmito` uses Maximum Likelihood Estimation (MLE) to fit Exponential Random
#' Graph Models for single or multiple small networks, the later using
#' pooled-data MLE. To do so we use exact likelihoods, which implies fully
#' enumerating the support of the model. Overall, the exact likelihood
#' calculation is only possible when dealing with directed (undirected) networks
#' size 5 (7). In general, directed (undirected) graphs with more than 5 (7)
#' vertices should not be fitted using MLE, but instead other methods such as
#' the MC-MLE algorithm or the Robbins-Monro Stochastic Approximation algorithm,
#' both of which are available in the ergm R package.The workhorse function of
#' `ergmito` is the `ergm` package function [ergm::ergm.allstats()].
#' 
#' @param object An object of class `ergmito`
#' @param model Model to estimate. See [ergm::ergm]. The only difference with
#' `ergm` is that the LHS can be a list of networks.
#' @param model_update A \code{\link[stats:formula]{formula}}. this can be used to
#' apply transformations, create interaction effects, add offset terms, etc. 
#' (see examples below and more details in [ergmito_formulae]).
#' @param optim.args List. Passed to [stats::optim].
#' @param target_stats A matrix of target statistics (see [ergm::ergm]).
#' @template stats
#' @template offset
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
#' @details 
#' The support of the sufficient statistics is calculated using ERGM's
#' [ergm::ergm.allstats()] function.
#' 
#' @seealso The function [plot.ergmito()] and [gof_ergmito()] for post-estimation
#' diagnostics.
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
#' Methods [base::print()], [base::summary()], [stats::coef()], [stats::logLik()],
#' [stats::nobs()], [stats::vcov()], [stats::AIC()], \code{\link[stats:AIC]{stats::BIC()}},
#' [stats::confint()], and  [stats::formula()] are available. 
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
#' 
#' @examples 
#' 
#' # Generating a small graph
#' set.seed(12)
#' n <- 4
#' net <- rbernoulli(n, p = .3)
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
#' # Example 4: Adding an reference term for edge-count ----------------------
#' 
#' # Simulating networks of different sizes
#' set.seed(12344)
#' nets <- rbernoulli(c(rep(4, 10), rep(5, 10)), c(rep(.2, 10), rep(.1, 10)))
#' 
#' # Fitting an ergmito under the Bernoulli model
#' ans0 <- ergmito(nets ~ edges)
#' summary(ans0)
#' # 
#' # ERGMito estimates
#' # 
#' # formula:
#' #   nets ~ edges
#' # 
#' #       Estimate Std. Error z value  Pr(>|z|)    
#' # edges -1.68640    0.15396 -10.954 < 2.2e-16 ***
#' # ---
#' # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#' # AIC: 279.3753    BIC: 283.1436    (Smaller is better.) 
#' 
#' 
#' # Fitting the model including a reference term for networks of size 5.
#' # Notice that the variable -n- and other graph attributes can be used
#' # with -model_update-.
#' ans1 <- ergmito(nets ~ edges, model_update = ~ I(edges * (n == 5)))
#' summary(ans1)
#' # 
#' # ERGMito estimates
#' # 
#' # formula:
#' #   nets ~ edges + I(edges * (n == 5))
#' # 
#' #                     Estimate Std. Error z value  Pr(>|z|)    
#' # edges               -1.18958    0.21583 -5.5116 3.556e-08 ***
#' # I(edges * (n == 5)) -0.90116    0.31250 -2.8837   0.00393 ** 
#' # ---
#' # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#' # AIC: 272.9916    BIC: 280.5282    (Smaller is better.) 
#' 
#' # The resulting parameter for the edge-count is smaller for networks
#' # of size five
#' plogis(coef(ans1)[1])   # 0.23
#' plogis(sum(coef(ans1))) # 0.11
#' 
#' # We can see that in this case the difference in edge-count matters.
#' if (require(lmtest)) {
#' 
#'   lrtest(ans0, ans1)
#'   # Likelihood ratio test
#'   # 
#'   # Model 1: nets ~ edges
#'   # Model 2: nets ~ edges + I(edges * (n == 5))
#'   # #Df  LogLik Df  Chisq Pr(>Chisq)   
#'   # 1   1 -138.69                        
#'   # 2   2 -134.50  1 8.3837   0.003786 **
#'   #   ---
#'   #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#' }
#' 
#' @importFrom stats optim terms rnorm
#' @importFrom MASS ginv
#' @name ergmito
NULL

ERGMITO_DEFAULT_OPTIM_CONTROL <- list(
  reltol = sqrt(.Machine$double.eps)
)

#' @export
#' @rdname ergmito
ergmito <- function(
  model,
  model_update   = NULL,
  stats_weights = NULL,
  stats_statmat = NULL,
  optim.args    = list(),
  init          = NULL,
  use.grad      = TRUE,
  target_stats  = NULL,
  ntries        = 1L,
  keep.stats    = TRUE,
  target_offset = NULL,
  stats_offset  = NULL,
  ...
  ) {

  # Keeping track of time
  timer <- c(start = unname(proc.time()["elapsed"]))

  # Generating the objective function
  ergmitoenv <- environment(model)
  
  formulae   <- ergmito_formulae(
    model,
    model_update  = model_update, 
    target_stats  = target_stats,
    stats_weights = stats_weights,
    stats_statmat = stats_statmat,
    target_offset = target_offset,
    stats_offset  = stats_offset,
    env           = ergmitoenv,
    ...
  )
  
  timer <- c(timer, ergmito_formulae = unname(proc.time()["elapsed"]))
  
  # Verifying existence of MLE
  support <- check_support(
    formulae$target_stats,
    formulae$stats_statmat
  )
  timer <- c(timer, check_support = unname(proc.time()["elapsed"]))
  
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
  
  optim.args$hessian <- FALSE

    
  # Setting arguments for optim
  optim.args$fn   <- formulae$loglik
  if (use.grad) 
    optim.args$gr <- formulae$grad
  
  optim.args$par  <- init
    
  # Will try to solve the problem more than once... if needed
  ntry <- 1L
  history <- matrix(
    NA,
    nrow = ntries,
    ncol = npars + 1,
    dimnames = list(
      1L:ntries,
      c(formulae$term_names, "value")
      )
  )
  
  # Passed (and default) other than the functions
  optim.args0 <- optim.args
  
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
    optim.args$par <- stats::rnorm(npars, -2, 2)
    ntry <- ntry + 1
    
  }

  timer <- c(timer, optim = unname(proc.time()["elapsed"]))
  
  # Checking the convergence
  estimates <- check_convergence(
    optim_output = ans,
    model        = formulae,
    support      = support
  )
  timer <- c(timer, chec_covergence = unname(proc.time()["elapsed"]))
  
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
    
    # We do it in this fashion as 
    to_keep <- setdiff(
        names(formulae),
        c(
          "stats_weights", "stats_statmat",
          "loglik", "grad",
          "target_offset", "stats_offset"
          )
        )

    fenvir <- environment(formulae$loglik)
    formulae <- formulae[to_keep]
    
    # We have to warn the user about this issue if he tries to reuse these functions
    formulae$loglik <- formulae$grad <- function(...) {
      stop(
        "As the option keep.stats = FALSE, this functions are no longer ",
        "available to the user. Re-run the model using keep.stats = TRUE, if you ",
        "want to use the loglikelihood or gradient functions of this model.",
        call. = FALSE
        )
    }
    ans$formulae <- formulae
    
    # Emptying environment, just to be safe
    rm(list = ls(envir = fenvir, all.names = TRUE), envir = fenvir)

  }
  
  # Counting cells, this will depend on whether nets are directed or not
  sizes    <- nvertex(ans$network)
  ans$nobs <- sizes * (sizes - 1)/
    ifelse(is_directed(ans$network), 1, 2)
  
  ans$nobs <- sum(ans$nobs)
  
  ans$timer <- diff(timer)
  ans$timer <- c(ans$timer, total = sum(ans$timer))
  ans
  
}
