#' Goodness of Fit diagnostics for ERGMito models
#' 
#' @param object An object of class [ergmito].
#' @param probs Numeric vector. Quantiles to plot (passed to [stats::quantile]).
#' @param R Integer scalar. Number of simulations to generate (passed to [sample]).
#' @param ... Further arguments passed to [stats::quantile].
#' 
#' @details 
#' The Goodness of Fit function uses the fitted ERGMito to calculate a given confidence
#' interval for a set of sufficient statistics. By default (and currently the
#' only available option), this is done on the sufficient statistics specified
#' in the model.
#' 
#' In detail, the algorithm is executed as follow:
#' 
#' For every network in the list of networks do:
#' 
#' 1. Calculate the probability of observing each possible graph in its support
#' using the fitted model.
#' 
#' 2. Draw `R` samples from each set of parameters using the probabilities computed
#' in (1).
#' 
#' 3. Using the `quantile` function, calculate the desired quantiles of the 
#' sufficient statistics.
#' 
#' The plot method is particularly convenient since it graphically shows whether
#' the target statistics of the model (observed statistics) fall within the 
#' simulated range.
#' 
#' @return An object of class `ergmito_gof`. This is a list with the following
#' components:
#' - `ci` A list of matrices of length `nnets(object)` with the corresponding
#' confidence intervals for the statistics of the model.
#' - `target.stats` A matrix of the target statistics.
#' - `ergmito.probs` A list of numeric vectors of length `nnets(object)` with the
#' probabilities associated to each possible structure of network.
#' - `probs` The value passed via `probs`.
#' - `model` The fitted model.
#' - `term.names` Character vector. Names of the terms used in the model.
#' - `quantile.args` A list of the values passed via `...`.
#' 
#' @examples 
#' # Fitting the fivenets model
#' data(fivenets, package = "ergmito")
#' fit <- ergmito(fivenets ~ edges + nodematch("female"))
#' 
#' # Calculating the gof
#' ans <- gof(fit)
#' 
#' # Looking at the results
#' ans
#' plot(ans)
#' @name ergmito_gof
#' @importFrom ergm gof
NULL

#' @export
#' @rdname ergmito_gof
gof <- function(...) UseMethod("gof")

#' @export
#' @rdname ergmito_gof
gof.ergmito <- function(
  object,
  probs = c(.05, .95),
  R     = 50000L,
  ...
  ) {
  
  res <- vector(mode = "list", length = nnets(object))
  pr  <- res
  ran <- res
  
  for (i in seq_len(length(res))) {
    
    # Computing the probability of observing each class of networks
    pr[[i]] <- exp(exact_loglik(
      x       = object$formulae$stats[[i]]$statmat,
      statmat = object$formulae$stats[[i]]$statmat,
      params  = stats::coef(object),
      weights = object$formulae$stats[[i]]$weights
    ))*object$formulae$stats[[i]]$weights
    
    # Calculating the quantiles. First, we need to make some room to the data
    # to be stored
    res[[i]] <- matrix(
      ncol = length(stats::coef(object)), nrow = 2,
      dimnames = list(
        sprintf("%4.1f%%", probs*100),
        names(stats::coef(object))
        )
      )
    
    # Makeing space for storing ranges of sufficient statistics
    ran[[i]] <- matrix(
      nrow = length(coef(object)), ncol = 3L,
      dimnames = list(names(coef(object)), c("min", "p50","max"))
      )
    
    # Generating confidence intervals for the statistics in the model
    for (k in seq_len(ncol(res[[i]]))) {
      
      # Sampling from the distribution (in the future we could do this
      # analytically instead)
      res_i <- sample(
        object$formulae$stats[[i]]$statmat[, k], size = R, prob = pr[[i]],
        replace = TRUE)
      
      # Calculating mins, p50, and max (this will be useful for plotting)
      ran[[i]][k, "p50"] <- stats::quantile(res_i, probs = 0.5)
      ran[[i]][k, c("min", "max")] <- range(object$formulae$stats[[i]]$statmat[, k])
      
      res[[i]][, k] <- stats::quantile(res_i, probs = probs, ...)
      
    }
    
  }
  
  structure(
    list(
      ci            = res,
      target.stats  = object$formulae$target.stats,
      ergmito.probs = pr,
      probs         = probs,
      model         = object$formulae$model,
      term.names    = colnames(res[[i]]),
      ranges        = ran,
      quantile.args = list(...)
      ),
    class = "ergmito_gof"
  )
  
}

#' @export
#' @rdname ergmito_gof
#' @details The print method tries to copy (explicitly) the print method of the
#' `gof` function from the `ergm` R package.
print.ergmito_gof <- function(x, ...) {
  
  K <- length(x$term.names)
  for (k in seq_len(K)) {
    cat("\nGoodness-of-fit for ", x$term.names[k], "\n\n")
    
    # Creating the table to be printed
    lower <- sapply(x$ci, "[", i = 1L, j = k, drop = TRUE)
    upper <- sapply(x$ci, "[", i = length(x$probs), j = k, drop = TRUE)
    minx  <- sapply(x$ranges, "[", i = k, j = "min", drop = TRUE)
    p50   <- sapply(x$ranges, "[", i = k, j = "p50", drop = TRUE)
    maxx  <- sapply(x$ranges, "[", i = k, j = "max", drop = TRUE)
    obs   <- x$target.stats[, k]
    
    tab <- data.frame(cbind(obs, minx, p50, maxx, lower, upper))
    dimnames(tab) <- list(
      paste("net", 1:length(x$ci)),
      c("obs", "min", "p50", "max", rownames(x$ci[[1]]))
      )
      
    print(tab)
    cat("\n")
    
  }
  
  invisible(x)
  
}

#' @export
#' @param x An object of class `ergmito_gof`.
#' @param y Ignored.
#' @param main Title of the plot
#' @rdname ergmito_gof
plot.ergmito_gof <- function(
  x,
  y    = NULL,
  main = "Goodness-of-fit Statistics",
  ...
  ) {
  
  K <- length(x$term.names)
  op <- graphics::par(
    mfcol = c(K, 1L), 
    mai   = graphics::par("mai")*c(.05, 1, .05, 1),
    omi   = graphics::par("mai")*c(1, 0, 1, 0)
    )
  on.exit(graphics::par(op))
  for (k in seq_len(K)) {
    
    lower <- sapply(x$ci, "[", i = 1L, j = k, drop = TRUE)
    upper <- sapply(x$ci, "[", i = length(x$probs), j = k, drop = TRUE)
    obs   <- x$target.stats[, k]
    
    # The ranges are obtained from the
    ran_min <- sapply(x$ranges, "[", i = k, j = "min")
    ran_max <- sapply(x$ranges, "[", i = k, j = "max")
    
    graphics::plot(obs, ylab = x$term.names[k], ylim = range(c(ran_min, ran_max)), xlab = "Network N",
         xaxt = ifelse(k != K, "n", "s"))
    graphics::lines(x = lower, lty = "dashed")
    graphics::lines(x = upper, lty = "dashed")
    graphics::lines(x = ran_min, lty = "solid", lwd=2, col="gray")
    graphics::lines(x = ran_max, lty = "solid", lwd=2, col="gray")
    
  }
  
  graphics::par(op)
  graphics::title(
    main =paste0(
      main, "\n",
      sprintf(
        "(Quantiles calculated at %4.1f%% and %4.1f%%)",
        x$probs[1]*100, x$probs[length(x$probs)]*100
        )
      ),
    xlab = "Network #",
    sub  = deparse(x$model)
    )
  
}


