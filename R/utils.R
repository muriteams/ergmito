
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

#' An alternative to `as.matrix` to retrieve adjacency matrix fast
#' 
#' This function does not perform significant checks. Furthermore, this function
#' won't keep the row/col names.
#' 
#' @param x An object to be coerced as an adjacency matrix.
#' @export
#' 
as_adjmat <- function(x) UseMethod("as_adjmat")

#' @export
#' @rdname as_adjmat
as_adjmat.network <- function(x) {
  
  n   <- nvertex(x)
  ans <- matrix(0L, nrow = n, ncol = n)
  
  if (x$gal$mnext == 1)
    return(ans)
  
  ties <- x$mel[sapply(x$mel, length) > 0L]
  
  ans[cbind(
    sapply(ties, "[[", "outl"),
    sapply(ties, "[[", "inl")
  )] <- 1L
  
  ans

}

#' @export
#' @rdname as_adjmat
as_adjmat.matrix <- function(x) x

#' @export
#' @rdname as_adjmat
as_adjmat.list <- function(x) {
  lapply(x, as_adjmat)
}

#' @export
#' @rdname as_adjmat
as_adjmat.formula <- function(x) {
  
  as_adjmat(eval(x[[2]], envir = environment(x)))
  
}

make_chunks <- function(N, chunk_size) {
  
  if (chunk_size > N)
    return(list(from=1, to=N))
  
  chunks <- seq(0, N, by = chunk_size)
  chunks <- list(from = chunks[-length(chunks)] + 1, to = chunks[-1])
  chunks$to[length(chunks$to)] <- N
  
  chunks
  
}