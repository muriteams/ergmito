#' Manipulation of network objects
#' 
#' This function implements a vectorized version of [network::network]`.adjmat`.
#' It allows us to turn regular matrices into network objects quickly.
#' 
#' @param x Either a single square matrix (adjacency matrix), or a list of these.
#' @param directed Logical scalar, if `FALSE` then the function only checks the
#' upper diagonal of the matrix assuming it is undirected.
#' @param loops Logical scalar. When `FALSE` (default) it will skip the diagonal
#' of the adjacency matrix.
#' @param hyper,multiple,bipartite Currently Ignored. Right now all the network objects
#' created by this function set these parameters as `FALSE`.
#' @param ... Further arguments passed to the method.
#' 
#' @details This version does not support adding the name parameter yet. The 
#' function in the network package includes the name of the vertices as an
#' attribute.
#' 
#' Just like in the network function, `NA` are checked and added accordingly, i.e.
#' if there is an `NA` in the matrix, then the value is recorded as a missing edge.
#' 
#' @return An object of class `network`. This is a list with the following elements:
#' - `mel` *Master Edge List*: A named list with length equal to the number of edges in
#' the network. The list itself has 3 elements: `inl` (tail), `outl` (head), and 
#' `atl` (attribute). By default `atl`, a list itself, has a single element: `na`.
#' 
#' - `gal` *Graph Attributes List*: a named list with the following elements:
#'   - `n` Number of nodes
#'   - `mnext` Number of edges + 1
#'   - `directed`,`hyper`,`loops`,`multiple`,`bipartite` The arguments passed to
#'   the function.
#'
#' - `val` *Vertex Attributes List*
#' 
#' - `iel` *In Edgest List*
#' 
#' - `oel` *Out Edgest List*
#' 
#' @examples 
#' set.seed(155)
#' adjmats  <- rbernoulli(rep(5, 20))
#' networks <- matrix_to_network(adjmats)
#' 
#' @export
matrix_to_network <- function(
  x, ...
  ) UseMethod("matrix_to_network")

#' @export
#' @rdname matrix_to_network
matrix_to_network.matrix <- function(
  x,
  directed  = rep(TRUE, length(x)),
  hyper     = rep(FALSE, length(x)),
  loops     = rep(FALSE, length(x)),
  multiple  = rep(FALSE, length(x)),
  bipartite = rep(FALSE, length(x)),
  ...
  ) {
  
  matrix_to_network_cpp(
    x = list(x),
    directed = directed,
    hyper = hyper,
    loops = loops,
    multiple = multiple,
    bipartite = bipartite
    )[[1L]]
  
}

#' @export
# @rdname matrix_to_network
matrix_to_network.list <- function(
  x,
  directed  = rep(TRUE, length(x)),
  hyper     = rep(FALSE, length(x)),
  loops     = rep(FALSE, length(x)),
  multiple  = rep(FALSE, length(x)),
  bipartite = rep(FALSE, length(x)),
  ...
  ) {
  
  # Checking the length of the attributes
  if (is.logical(directed) && (length(directed) == 1L))
    directed <- rep(directed, nnets(x))
  if (is.logical(hyper) && (length(hyper) == 1L))
    hyper <- rep(hyper, nnets(x))
  if (is.logical(loops) && (length(loops) == 1L))
    loops <- rep(loops, nnets(x))
  if (is.logical(multiple) && (length(multiple) == 1L))
    multiple <- rep(multiple, nnets(x))
  if (is.logical(bipartite) && (length(bipartite) == 1L))
    bipartite <- rep(bipartite, nnets(x))
  
  # Checking matching lengths
  if (length(directed) != nnets(x))
    stop("The length of -directed- doesn't matches the number of networks in -x-.")
  if (length(hyper) != nnets(x))
    stop("The length of -hyper- doesn't matches the number of networks in -x-.")
  if (length(loops) != nnets(x))
    stop("The length of -loops- doesn't matches the number of networks in -x-.")
  if (length(multiple) != nnets(x))
    stop("The length of -multiple- doesn't matches the number of networks in -x-.")
  if (length(bipartite) != nnets(x))
    stop("The length of -bipartite- doesn't matches the number of networks in -x-.")
  
  # Checking all are any of the types
  if (inherits(x, "network"))
    return(x)
  
  # Is this a list of networks?
  if (all(sapply(x, inherits, what = "network")))
    return(x)
  
  if (!all(sapply(x, inherits, what = "matrix")))
    stop(
      "When passing lists, all objects have to be either a list of network ",
      "objects, or a list of matrices.", call. = FALSE
      )
  
  # To save memory, we do this by chunks
  chunks <- make_chunks(nnets(x), 50000)
  
  res <- vector("list", nnets(x))
  for (i in seq_along(chunks$from)) {
    
    res[chunks$from[i]:chunks$to[i]] <- matrix_to_network_cpp(
      x         = x[chunks$from[i]:chunks$to[i]],
      directed  = directed[chunks$from[i]:chunks$to[i]],
      hyper     = hyper[chunks$from[i]:chunks$to[i]],
      loops     = loops[chunks$from[i]:chunks$to[i]],
      multiple  = multiple[chunks$from[i]:chunks$to[i]],
      bipartite = bipartite[chunks$from[i]:chunks$to[i]]
    )
    
  }
  
  res
  
}

#' 
#' #' Alter network objects
#' #' @param x,attrvalie,attrname see [network::set.vertex.attribute]
#' #' @export
#' add_vertex_attr <- function(x, attrvalue, attrname) {
#'   
#'   if (is.network(x))
#'     x <- list(x)
#'   
#'   if (!all(sapply(x, inherits, "network")))
#'     stop("`x` must be an object of class network.", call. = FALSE)
#'   
#'   add_vertex_attr_cpp(x, attrvalue, attrname, names(x[[1]]$val[[1]]))
#'   
#'   # x
#'   
#' }
