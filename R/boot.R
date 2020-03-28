
#' Bootstrap of ergmito
#' 
#' @param x Either a formula or an object of class [ergmito].
#' @param ... Additional arguments passed to the method.
#' @param R Integer. Number of replicates
#' @param ncpus Integer Number of CPUs to use.
#' @param cl An object of class `cluster` (see [parallel::makePSOCKcluster])
#'  
#' @export
#' @importFrom parallel makePSOCKcluster stopCluster clusterEvalQ clusterExport
#' @importFrom stats update.formula var
#' 
#' @details The resulting sample of parameters estimates is then used to compute
#' the variance-covariance matrix of the model. Cases in which `Inf`/`NaN`/`NA`
#' values were returned are excluded from the calculation.
#' 
#' @return An object of class `ergmito_boot` and [ergmito]
#' 
#' @examples 
#' 
#' # Simulating 20 bernoulli networks of size 4
#' nets <- replicate(20, rbernoulli(4), simplify = FALSE)
ergmito_boot <- function(x, ..., R, ncpus = 1L, cl = NULL) UseMethod("ergmito_boot")


#' @export
# @rdname ergmito_boot
ergmito_boot.formula <- function(x, ..., R, ncpus = 1L, cl = NULL) {
  
  # First run of the model
  x0 <- ergmito(x, ...)
  
  # Now we do the bootstrap
  ergmito_boot(x = x0, R = R, ncpus = ncpus, cl = cl)
  
}

#' @export
# @rdname ergmito_boot
ergmito_boot.ergmito <- function(x, ..., R, ncpus = 1L, cl = NULL) {
  
  n <- nnets(x)
  
  if (n < 2)
    stop(
      "I wish I could do bootstrapping with less than 2, but I can't.",
      call. = FALSE)
  else if (n <= 10)
    warning_ergmito(
      "You are doing bootstrapping with less than 10 networks (and even 10 is too few).",
      call.      = FALSE,
      immediate. = TRUE
      )
  
  # Getting the sample, and baseline model
  IDX           <- replicate(n = R, sample.int(n, n, TRUE), simplify = FALSE)
  model0        <- stats::update.formula(x$formulae$model_final, nets0[idx] ~ .)
  target_stats  <- x$formulae$target_stats
  nets0         <- x$network
  stats_weights <- x$formulae$stats_weights
  stats_statmat <- x$formulae$stats_statmat
  target_offset <- x$formulae$target_offset
  stats_offset  <- x$formulae$stats_offset
  
  # Creating the cluster and setting the seed
  if (ncpus > 1L && !length(cl)) {
    
    # Setting up the cluster
    cl <- parallel::makePSOCKcluster(ncpus)
    parallel::clusterEvalQ(cl, library(ergmito))
    parallel::clusterExport(
      cl, c(
        "model0", "target_stats", "nets0", "stats_weights", "stats_statmat",
        "target_offset", "stats_offset"
        ),
      envir = environment())
    
    parallel::clusterSetRNGStream(cl)

    # We must stop it afterwards
    on.exit(parallel::stopCluster(cl))
    
  }
  
  # Running the estimation process
  boot_estimates <- if (ncpus > 1L) {
    
    parallel::parLapply(cl, IDX, function(idx) {
      
      environment(model0) <- environment()
      ans <- tryCatch(ergmito(
        model0,
        stats_weights = stats_weights[idx],
        stats_statmat = stats_statmat[idx],
        target_stats  = target_stats[idx,,drop=FALSE],
        target_offset = target_offset[idx],
        stats_offset  = stats_offset[idx]
        ), error = function(e) e
        )
      
      if (inherits(ans, "error"))
        return(NULL)
      
      stats::coef(ans)
      
    })
    
  } else {
    
    lapply(IDX, function(idx) {
      
      environment(model0) <- environment()
      ans <- tryCatch(suppressWarnings(ergmito(
        model0,
        stats_weights = stats_weights[idx],
        stats_statmat = stats_statmat[idx],
        target_stats  = target_stats[idx,,drop=FALSE],
        target_offset = target_offset[idx],
        stats_offset  = stats_offset[idx]
        )), error = function(e) e)
      
      if (inherits(ans, "error"))
        return(NULL)
      
      stats::coef(ans)
    })
    
  }
  
  # Computing variance/covariance matrix
  coefs <- do.call(rbind, boot_estimates)
  
  # Tagging finite results
  are_finate <- unique(which(!is.finite(coefs), arr.ind = TRUE)[, 1L, drop = FALSE])
  are_finate <- setdiff(seq_along(IDX), are_finate)

  if (length(are_finate) < 2) {
    stop("At most one of the replicates had finite estimates (not Inf/NA/NaN).",
         call. = FALSE)
  } 
  
  x$covar  <- stats::var(coefs[are_finate, , drop = FALSE])
  x$call   <- match.call()
  x$R      <- R
  x$sample <- IDX
  x$dist   <- coefs
  
  class(x) <- c("ergmito_boot", class(x))
  
  x
  
}

#' @export
# @rdname ergmito_boot
print.ergmito_boot <- function(x, ...) {
  
  cat(sprintf("Bootstrapped %i replicates:\n", x$R))
  print(structure(unclass(x), class="ergmito"))
  
  invisible(x)
  
}

