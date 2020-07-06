
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
# @rdname nvertex
nedges.network <- function(x, ...) {
  network::network.edgecount(x, ...)
}

#' @export
# @rdname nvertex
nedges.list <- function(x, ...) {
  sapply(x, nedges, ...)
}

#' @export
# @rdname nvertex
nedges.matrix <- function(x, ...) {
  sum(x != 0)
}

#' @export
# @rdname nvertex
nedges.ergmito <- function(x, ...) {
  nedges(x$network, ...)
}

#' @export
# @rdname nvertex
nedges.formula <- function(x, ...) {
  nedges(eval(x[[2]]), envir = environment(x))
}

#' @export
# @rdname nvertex
nvertex.network <- function(x) {
  
  network::network.size(x)
  
}

#' @export
# @rdname nvertex
nvertex.matrix <- function(x) {
  
  ncol(x)
  
}

#' @export
# @rdname nvertex
nvertex.list <- function(x) {
  
  sapply(x, nvertex)
  
}

#' @export
# @rdname nvertex
nvertex.ergmito <- function(x) {
  
  nvertex(x$network)
  
}

#' @export
# @rdname nvertex
nvertex.formula <- function(x) {
  nvertex(eval(x[[2]], envir = environment(x)))
}

#' @export
#' @rdname nvertex
nnets <- function(x) UseMethod("nnets")

#' @export
# @rdname nvertex
nnets.list <- function(x) length(x)

#' @export
# @rdname nvertex
nnets.matrix <- function(x) 1L

#' @export
# @rdname nvertex
nnets.network <- function(x) 1L

#' @export
# @rdname nvertex
nnets.ergmito <- function(x) {
  
  x$formulae$nnets
  
}

#' @export
# @rdname nvertex
nnets.formula <- function(x) {
  
  nnets(eval(x[[2]], envir = environment(x)))
  
}

#' @export
#' @rdname nvertex
#' @param check_type Logical scalar. When checking for whether the network is
#' directed or not, we can ask the function to return with an error if what we
#' are checking is not an object of class network, otherwise it simply returns
#' false.
#' @return `is_directed` checks whether the passed networks are directed using
#' the function \code{\link[network:network.indicators]{is.directed}}. In the case of multiple networks,
#' the function returns a logical vector. Only objects of class `network` can be
#' checked, otherwise, if `check_type = FALSE`, the function returns `TRUE` by default.
#' @examples 
#' set.seed(771)
#' net <- lapply(rbernoulli(c(4, 4)), network::network, directed = FALSE)
#' is_directed(net)
#' is_directed(net[[1]])
#' is_directed(net ~ edges)
#' \dontrun{
#'   is_directed(net[[1]][1:4, 1:4], check_type = TRUE) # Error
#' }
#' is_directed(net[[1]][1:4, 1:4])
is_directed <- function(x, check_type = FALSE) UseMethod("is_directed")

#' @export
# @rdname nvertex
is_directed.network <- function(x, check_type = FALSE) network::is.directed(x)

#' @export
# @rdname nvertex
is_directed.list <- function(x, check_type = FALSE) {
  sapply(x, is_directed, check_type = check_type)
}

#' @export
# @rdname nvertex
is_directed.default <- function(x, check_type = FALSE) {
  
  if (check_type) 
    stop(
      "Only objects of class `network` or `ergmito` can be checked for directed.",
      call. = FALSE
      )
  else if (inherits(x, "list"))
    return(rep(TRUE, length(x)))
  else
    return(TRUE)
}

#' @export
# @rdname nvertex
is_directed.ergmito <- function(x, check_type = FALSE) {
  is_directed(x$network, check_type = check_type)
}

#' @export
# @rdname nvertex
is_directed.formula <- function(x, check_type = FALSE) {
  is_directed(eval(x[[2]], envir = environment(x)), check_type = check_type)
}

#' An alternative to `as.matrix` to retrieve adjacency matrix fast
#' 
#' This function does not perform significant checks. Furthermore, this function
#' won't keep the row/col names.
#' 
#' @param x An object to be coerced as an adjacency matrix.
#' @export
#' @examples 
#' 
#' set.seed(1231)
#' x <- matrix_to_network(rbernoulli(rep(5, 100)))
#' benchmarkito(
#'   as_adjmat = as_adjmat(x),
#'   as.matrix = lapply(x, as.matrix)
#' )
as_adjmat <- function(x) UseMethod("as_adjmat")

#' @export
# @rdname as_adjmat
as_adjmat.network <- function(x) {
  
  n   <- nvertex(x)
  ans <- matrix(0L, nrow = n, ncol = n)
  
  if (x$gal$mnext == 1)
    return(ans)
  
  ties <- x$mel[sapply(x$mel, length) > 0L]
  
  if (length(ties))
    ans[cbind(
      sapply(ties, "[[", "outl"),
      sapply(ties, "[[", "inl")
    )] <- 1L
  
  # In the case of undirected networks, we need to pass the
  # edges specifically
  if (!network::is.directed(x)) {
    ans <- ans + t(ans)
    ans[ans != 0] <- 1
  }
  
  ans

}

#' @export
# @rdname as_adjmat
as_adjmat.matrix <- function(x) x

#' @export
# @rdname as_adjmat
as_adjmat.list <- function(x) {
  lapply(x, as_adjmat)
}

#' @export
# @rdname as_adjmat
as_adjmat.formula <- function(x) {
  
  as_adjmat(eval(x[[2]], envir = environment(x)))
  
}

make_chunks <- function(N, chunk_size) {
  
  if (N < 1)
    stop(
      "N should be an integer greater than 0.",
      call. = FALSE
      )
  
  if (chunk_size > N)
    return(list(from=1, to=N))
  
  chunks <- seq(0, N, by = chunk_size)
  chunks <- list(from = chunks[-length(chunks)] + 1, to = chunks[-1])
  chunks$to[length(chunks$to)] <- N
  
  chunks
  
}

#' We don't need to warn during the 
#' @noRd
warning_ergmito <- function(...) {
  
  if (getOption("ergmito_warning", TRUE) == FALSE)
    return()
  
  warning(...)
}
