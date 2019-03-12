
#' Utility functions to query network dimensions
#' @param x Either an object of class [ergmito], [network], [formula], or [matrix].
#' @param ... Further arguments passed to the method. Currently only `nedges.network`
#' receives arguments (see [network::network.edgecount]).
#' @export
nvertex <- function(x) UseMethod("nvertex")

#' @export
#' @rdname nvertex
nedges <- function(x, ...) UseMethod("nedges")

#' @export
#' @rdname nvertex
nedges.network <- function(x, ...) {
  network::network.edgecount(x, ...)
}

#' @export
#' @rdname nvertex
nedges.list <- function(x, ...) {
  sapply(x, nedges, ...)
}

#' @export
#' @rdname nvertex
nedges.matrix <- function(x, ...) {
  sum(x != 0)
}

#' @export
#' @rdname nvertex
nedges.ergmito <- function(x, ...) {
  nedges(x$network, ...)
}

#' @export
#' @rdname nvertex
nedges.formula <- function(x, ...) {
  nedges(eval(x[[2]]), envir = environment(x))
}

#' @export
#' @rdname nvertex
nvertex.network <- function(x) {
  
  network::network.size(x)
  
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
  
  nvertex(x$network)
  
}

#' @export
#' @rdname nvertex
nvertex.formula <- function(x) {
  nvertex(eval(x[[2]], envir = environment(x)))
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

#' @export
#' @rdname nvertex
nnets.formula <- function(x) {
  
  nnets(eval(x[[2]], envir = environment(x)))
  
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


#' An alternative to `as.matrix` to retrieve adjacency matrix fast
#' 
#' This function does not perform significant checks. Furthermore, this function
#' won't keep the row/col names.
#' 
#' @param x An object to be coerced as an adjacency matrix.
#' @export
#' 
as.adjmat <- function(x) UseMethod("as.adjmat")

#' @export
#' @rdname as.adjmat
as.adjmat.network <- function(x) {
  
  n   <- nvertex(x)
  ans <- matrix(0, nrow = n, ncol = n)
  
  ties <- x$mel[sapply(x$mel, length) > 0L]
  
  ans[cbind(
    sapply(ties, "[[", "outl"),
    sapply(ties, "[[", "inl")
  )] <- 1L
  
  ans

}

#' @export
#' @rdname as.adjmat
as.adjmat.list <- function(x) {
  lapply(x, as.adjmat)
}

#' @export
#' @rdname as.adjmat
as.adjmat.formula <- function(x) {
  
  as.adjmat(eval(x[[2]], envir = environment(x)))
  
}