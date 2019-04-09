#' Goodness of Fit diagnostics for ERGMito models
#' 
#' @param object An object of class [ergmito].
#' @param probs Numeric vector. Quantiles to plot (see details).
#' @param sim_ci Logical scalar. If `FALSE`, the default, it will compute the
#' quantiles analytically, otherwise it samples from the ERGM distribution.
#' @param R Integer scalar. Number of simulations to generate (passed to [sample]).
#' This is only used if `sim_ci = TRUE`.
#' @param GOF Formula. Additional set of parameters to perform the GOF.
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
#' 2. If `sim_ci = TRUE`, draw `R` samples from each set of parameters using the
#' probabilities computed. Then using the `quantile` function, calculate the desired
#' quantiles of the sufficient statistics. Otherwise, compute the quantiles using
#' the analytic quantiles using the full distribution.' 
#' 
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
#' ans <- gof_ergmito(fit)
#' 
#' # Looking at the results
#' ans
#' plot(ans)
#' @name ergmito_gof
NULL

#' Given a vector of probabilities, find the feasible lower and upper bound
#' @param probs A numeric vector (should add up to 1).
#' @param lower Numeric scalar. Lower bound within 0 and .5.
#' @param upper Numeric scalar. Upper bound within .5 and 1.
#' @return A named vector with the feasible upper and lower bounds. The names of
#' the elements are mapped to the probabilities.
#' @noRd
get_feasible_ci_bounds <- function(x, probs, lower, upper) {
  
  # Checking range
  if (lower > .5)
    stop("`lower` cannot be more than .5", call. = FALSE)
  if (upper < .5)
    stop("`lower` cannot be less than .5", call. = FALSE)
  
  # Collapsing the ranges, first we need to sort these accordignly.
  # Notice that we need to do this since the probs passed to this
  # function are computed based on the graph, not in the statistic
  # itself.
  idx   <- order(x)
  x     <- x[idx]
  probs <- probs[idx]
  probs <- stats::aggregate(probs ~ x, FUN = sum)
  x     <- unname(probs[, 1L, drop = TRUE])
  probs <- unname(probs[, 2L, drop = TRUE])
  
  cumprobs <- cumsum(probs)
  n <- length(cumprobs)
  
  if (abs(cumprobs[n] - 1.0) > 1e-10)
    stop("`probs` does not add up to 1.", call. = FALSE)

  # Finding the feasible lower bound
  for (i in 1L:n)
    if (lower < cumprobs[i]) {
      
      if (i > 1L)
        i <- i - 1L
      
      break
    }
  
  # Finding the feasible upper bound
  for (j in n:1L)
    if (upper > cumprobs[j]) {
      
      if (j < n)
        j <- j + 1L
      
      break
    }
      
  c(x[c(i, j)], cumprobs[c(i, j)])
  
}



#' @export
#' @rdname ergmito_gof
gof_ergmito <- function(
  object,
  GOF    = NULL,
  probs  = c(.05, .95),
  sim_ci = FALSE,
  R      = 50000L,
  ...
  ) {
  
  res <- vector(mode = "list", length = nnets(object))
  pr  <- res
  ran <- res
  
  if (!is.null(GOF)) {
    # Deparsing formula
    if (!inherits(GOF, "formula"))
      stop("`GOF` must be an object of class formula.", call. = FALSE)
    
    # Analyzing terms
    test <- any(attr(stats::terms(GOF), "term.labels") %in% 
                  attr(stats::terms(object$formulae$model), "term.labels"))
    if (test)
      stop(
        "The `GOF` argument must specify terms others than those included in the model.",
        call. = FALSE
      )
    
    # This formula includes everything
    GOF <- stats::as.formula(gsub(".*[~]", "~ . + ", deparse(GOF)))
    GOF <- stats::update.formula(object$formulae$model, GOF)
  }
  
  target.stats <- NULL
  
  for (i in seq_len(length(res))) {
    
    # Has the user provided a formula? In this case we need to use an alternative
    # matrix of `statmat` and `weights`. We are restricted to whatever allstats
    # allows generating (e.g. distance is not included)
    if (!is.null(GOF)) {
      
      GOF <- stats::update.formula(GOF, object$network[[i]] ~ .)
      environment(GOF) <- environment()
      
      statmat. <- ergm::ergm.allstats(GOF, zeroobs = FALSE) 
      weights. <- statmat.$weights
      statmat. <- statmat.$statmat
      
      target.stats <- rbind(target.stats, summary(GOF))
      
    } else {
      
      statmat. <- object$formulae$stats[[i]]$statmat
      weights. <- object$formulae$stats[[i]]$weights
      
      target.stats <- rbind(target.stats, object$formulae$target.stats[i, ])
    }
    
    # Computing the probability of observing each class of networks
    pr[[i]] <- exp(exact_loglik(
      x       = statmat.[, names(stats::coef(object)), drop = FALSE],
      statmat = statmat.[, names(stats::coef(object)), drop = FALSE],
      params  = stats::coef(object),
      weights = weights.
    ))*weights.
    
    # Calculating the quantiles. First, we need to make some room to the data
    # to be stored
    res[[i]] <- matrix(
      nrow = ncol(statmat.), ncol = 4,
      dimnames = list(
        colnames(statmat.),
        c("lower-q", "upper-q", "lower-p", "upper-p")
        )
      )
    
    # Makeing space for storing ranges of sufficient statistics
    ran[[i]] <- matrix(
      nrow = nrow(res[[i]]), ncol = 3L,
      dimnames = list(colnames(statmat.), c("min", "mean","max"))
      )
    
    # Generating confidence intervals for the statistics in the model. We do this
    # by calculating the exact CI based on the data.
    for (k in seq_len(nrow(res[[i]]))) {
      
      # Sampling from the distribution (in the future we could do this
      # analytically instead)
      if (sim_ci) {
        
        res_i <- sample(statmat.[, k], size = R, prob = pr[[i]], replace = TRUE)
        
        res_i <- unname(c(
          do.call(
            stats::quantile,
            c(list(x = res_i, probs = probs), list(...))
            ),
          probs))
        
        
      } else {
        
        # Finding the feasible bounds, based on the ordering of this data
        res_i <- get_feasible_ci_bounds(
          x     = statmat.[,k], 
          probs = pr[[i]],
          lower = probs[1],
          upper = probs[2]
          )
        
      }
      
      # Calculating mins, mean, and max (this will be useful for plotting)
      ran[[i]][k, "mean"] <- sum(statmat.[, k] * pr[[i]])
      ran[[i]][k, c("min", "max")] <- range(statmat.[, k])
      
      # This res_i vector contains the lower/upper bounds and the associated
      # probabiities.
      res[[i]][k, ] <- res_i
      
    }
    
  }
  
  structure(
    list(
      ci            = res,
      target.stats  = target.stats,
      ergmito.probs = pr,
      probs         = probs,
      model         = if (is.null(GOF))
        object$formulae$model
      else {
        as.formula(gsub(".*[~]", " ~ ", deparse(GOF)))
        },
      term.names    = rownames(res[[i]]),
      ranges        = ran,
      sim_ci        = sim_ci,
      quantile.args = list(...)
      ),
    class = "ergmito_gof"
  )
  
}

#' @export
#' @param digits Number of digits to used when printing
#' @rdname ergmito_gof
#' @details The print method tries to copy (explicitly) the print method of the
#' `gof` function from the `ergm` R package.
print.ergmito_gof <- function(x, digits = 2L, ...) {
  
  K <- length(x$term.names)
  for (k in seq_len(K)) {
    cat("\nGoodness-of-fit for", x$term.names[k], "\n\n", sep=" ")
    
    # Creating the table to be printed
    lower  <- sapply(x$ci, "[", i = k, j = "lower-q", drop = TRUE)
    upper  <- sapply(x$ci, "[", i = k, j = "upper-q", drop = TRUE)
    lowerp <- sapply(x$ci, "[", i = k, j = "lower-p", drop = TRUE)
    upperp <- sapply(x$ci, "[", i = k, j = "upper-p", drop = TRUE)
    minx   <- sapply(x$ranges, "[", i = k, j = "min", drop = TRUE)
    mean   <- sapply(x$ranges, "[", i = k, j = "mean", drop = TRUE)
    maxx   <- sapply(x$ranges, "[", i = k, j = "max", drop = TRUE)
    obs    <- x$target.stats[, k]
    
    tab <- data.frame(cbind(obs, minx, mean, maxx, lower, upper, lowerp, upperp))
    dimnames(tab) <- list(
      paste("net", 1:length(x$ci)),
      c("obs", "min", "mean", "max", "lower", "upper", "lower prob.", "upper prob.")
      )
      
    print(tab, digits = digits)
    cat("\n")
    
  }
  
  if (!x$sim_ci)
    cat(
      "Note: Exact confidence intervals where used.",
      "This implies that the requestes CI may differ from the one used",
      "(see ?gof_ergmito).\n\n"
      )
  else
    cat(
      "Note: Approximated confidence intervals where used (see ?gof_ergmito).\n\n"
      )
  
  
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
    
    lower <- sapply(x$ci, "[", i = k, j = "lower-q", drop = TRUE)
    upper <- sapply(x$ci, "[", i = k, j = "upper-q", drop = TRUE)
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
        "(Quantiles to cover at least a %4.1f%% CI)",
        x$probs[length(x$probs)]*100 - x$probs[1]*100 
        )
      ),
    xlab = "Network #",
    sub  = deparse(x$model)
    )
  
}


