
#' Bootstrap of lergm
#' 
#' @param x Either a formula or an object of class [lergm].
#' @param ... Additional arguments passed to the method.
#' @param R Integer. Number of replicates
#' @param ncpus Integer Number of CPUs to use.
#' @param cl An object of class `cluster` (see [parallel::makePSOCKcluster])
#'  
#' @export
#' @importFrom parallel makePSOCKcluster stopCluster clusterEvalQ clusterExport
#' @importFrom stats update.formula var
#' @examples 
#' 
#' # Simulating 20 bernoulli networks of size 4
#' nets <- replicate(20, rbernoulli(4), simplify = FALSE)
lergm_boot <- function(x, ..., R, ncpus = 1L, cl = NULL) UseMethod("lergm_boot")


#' @export
#' @rdname lergm_boot
lergm_boot.formula <- function(x, ..., R, ncpus = 1L, cl = NULL) {
  
  # First run of the model
  x0 <- lergm(x, ...)
  
  # Now we do the bootstrap
  lergm_boot(x = x0, R = R, ncpus = ncpus, cl = cl)
  
}

#' @export
#' @rdname lergm_boot
lergm_boot.lergm <- function(x, ..., R, ncpus = 1L, cl = NULL) {
  
  n <- nnets(x)
  
  # Getting the sample, and baseline model
  IDX    <- replicate(n = R, sample.int(n, n, TRUE), simplify = FALSE)
  model0 <- stats::update.formula(x$formulae$model, nets. ~ .)
  stats0 <- x$formulae$stats
  nets0  <- x$network
  obs_stats0 <- x$formulae$obs_stats
  
  # Creating the cluster and setting the seed
  if (ncpus > 1L && !length(cl)) {
    
    # Setting up the cluster
    cl <- parallel::makePSOCKcluster(ncpus)
    parallel::clusterEvalQ(cl, library(lergm))
    parallel::clusterExport(cl, c("model0", "stats0", "nets0"), envir = environment())
    
    parallel::clusterSetRNGStream(cl)

    # We must stop it afterwards
    on.exit(parallel::stopCluster(cl))
    
  }
  
  # Running the estimation process
  boot_estimates <- if (ncpus > 1L) {
    
    parallel::parLapply(cl, IDX, function(idx) {
      
      nets.  <- nets0[idx]
      lergm(model0, stats = stats0[idx], obs_stats = obs_stats0[idx])
      
    })
    
  } else {
    
    lapply(IDX, function(idx) {
      
      nets.  <- nets0[idx] 
      lergm(model0, stats = stats0[idx], obs_stats = obs_stats0[idx])
      
    })
    
  }
  
  # Computing variance/covariance matrix
  coefs <- do.call(rbind, lapply(boot_estimates, stats::coef))
  
  x$covar  <- stats::var(coefs)
  x$call   <- match.call()
  x$R      <- R
  x$sample <- IDX
  x$dist   <- coefs
  
  class(x) <- c("lergm_boot", class(x))
  
  x
  
}

#' @export
#' @rdname lergm_boot
print.lergm_boot <- function(x, ...) {
  
  cat(sprintf("Bootstrapped %i replicates:\n", x$R))
  print(structure(unclass(x), class="lergm"))
  
}

