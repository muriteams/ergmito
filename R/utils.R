
#' Utility functions to query network dimensions
#' @param x Either an object of class [lergm], [network], or [matrix].
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
nvertex.lergm <- function(x) {
  
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
nnets.lergm <- function(x) {
  
  x$formulae$nnets
  
}