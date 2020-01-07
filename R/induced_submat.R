#' Extract a submatrix from a network
#' 
#' This is similar to [network::get.inducedSubgraph]. The main difference is that
#' the resulting object will always be a list of matrices, and it is vectorized.
#' 
#' @param x Either a list or single matrices or network objects.
#' @param v Either a list or a single integer vector of vertices to subset.
#' @param ... Currently ignored.
#' 
#' @details Depending on the lengths of `x` and `v`, the function can take the
#' following strategies:
#' 
#' - If both are of the same size, then it will match the networks and the vector
#' of indices.
#' 
#' - If `length(x) == 1`, then it will use that single network as a baseline
#' for generating the subgraphs.
#' 
#' - If `length(v) == 1`, then it will generate the subgraph using the same set
#' of vertices for each network.
#' 
#' - If both have more than one element, but different sizes, then the function
#' returns with an error.
#' 
#' @examples 
#' x <- rbernoulli(100)
#' induced_submat(x, c(1, 10, 30:50))
#' 
#' x <- rbernoulli(c(20, 20))
#' induced_submat(x, c(1:10))
#' 
#' @return 
#' A list of matrices as a result of the subsetting.
#' @export
induced_submat <- function(x, v, ...) UseMethod("induced_submat")

#' @rdname induced_submat
#' @export
induced_submat.list <- function(x, v, ...) {
  
  if (!is.list(v))
    v <- list(v)
  
  x <- as_adjmat(x)
  
  induced_submat.(x, v)
  
}

#' @rdname induced_submat
#' @export
induced_submat.network <- function(x, v, ...) {
  
  x <- list(as_adjmat(x))
  if (!is.list(v))
    v <- list(v)
  induced_submat.(x, v, ...)
  
}

#' @rdname induced_submat
#' @export
induced_submat.matrix <- function(x, v, ...) {
  
  if (!is.list(v))
    v <- list(v)
  
  induced_submat.(list(x), v, ...)
  
}