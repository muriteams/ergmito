
#' Utility functions to query network dimensions
#' @param x Either an object of class [ergmito], [network], or [matrix].
#' @export
nvertex <- function(x) UseMethod("nvertex")

#' @export
#' @rdname nvertex
nvertex.network <- function(x) {
  
  x$gal$n
  
}

#' @export
#' @rdname nvertex
nvertex.matrix <- function(x) {
  
  ncol(x)
  
}

#' @export
#' @rdname nvertex
nvertex.list <- function(x) {
  
  sapply(x, nvertex)
  
}

#' @export
#' @rdname nvertex
nvertex.ergmito <- function(x) {
  
  if (nnets(x) == 1)
    nrow(x$network)
  else
    sapply(x$network, nrow)
  
}

#' @export
#' @rdname nvertex
nnets <- function(x) UseMethod("nnets")

#' @export
#' @rdname nvertex
nnets.list <- function(x) length(x)

#' @export
#' @rdname nvertex
nnets.matrix <- function(x) 1L

#' @export
#' @rdname nvertex
nnets.network <- function(x) 1L

#' @export
#' @rdname nvertex
nnets.ergmito <- function(x) {
  
  x$formulae$nnets
  
}

#' Extract function to be used with the `texreg` package.
#'
#' If available, this function can be used to generate nice looking tables
#' of little ERGMS.
#' @param model An object of class `ergmito`.
#' @param include.aic,include.bic,include.loglik See [texreg::extract].
#' @param include.nnets Logical. When true, it adds the Number of networks used
#' to the list of gof statistics. This can be useful when running pooled models.
#' @param ... Further arguments passed to [summary.ergmito].
#' @export
#'
extract.ergmito <- function(
  model, include.aic = TRUE, include.bic = TRUE, include.loglik = TRUE,
  include.nnets = TRUE,
  ...
  ) {

  if (!requireNamespace("texreg", quietly = TRUE))
    stop("Need to install the `texreg` package.", call. = FALSE)


  # Copied from texreg::extract.ergm

  s <- summary(model, ...)
  coefficient.names <- rownames(s$coefs)
  coefficients <- s$coefs[, 1]
  standard.errors <- s$coefs[, "Std. Error"]
  significance <- s$coefs[, "Pr(>|z|)"]
  
  # Adding gof statistics ------------------------------------------------------
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE && !is.null(s$aic)) {
    aic <- s$aic
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE && !is.null(s$bic)) {
    bic <- s$bic
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE && !is.null(model$mle.lik)) {
    lik <- model$mle.lik[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  # Number of networks 
  if (include.nnets) {
    gof         <- c(gof, nnets(model$network))
    gof.names   <- c(gof.names, "# Networks")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- texreg::createTexreg(coef.names = coefficient.names, coef = coefficients,
                     se = standard.errors, pvalues = significance, gof.names = gof.names,
                     gof = gof, gof.decimal = gof.decimal)
  return(tr)

}

setMethod(
  "extract", signature = className("ergmito", "ergmito"),
  definition = extract.ergmito
  )
