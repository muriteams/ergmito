#' Block-diagonal models using `ergm`
#' 
#' These two functions are use to go back and forth from a pooled ergm vs a
#' blockdiagonal model, the latter to be fitted using [ergm::ergm]. 
#' 
#' @param x In the case of `blockdiagonalize`, a list of networks or matrices.
#' For `splitnetwork` a single network object with a vertex attribute that can
#' be used to split the data.
#' @param attrname Name of the attribute that holds the block ids.
#' @examples 
#' \dontrun{
#' library(ergm)
#' data(fivenets)
#' 
#' fivenets2 <- blockdiagonalize(fivenets, attrname = "block") # A network with
#' ans0 <- ergm(
#'   fivenets2 ~ edges + nodematch("female"),
#'   constrains = ~blockdiag("block")
#'   )
#' ans1 <- ergmito(fivenets ~ edges + nodematch("female"))
#' 
#' }
#' @export
blockdiagonalize <- function(x, attrname = "block") {
  
  sizes <- nvertex(x)
  N     <- sum(sizes)
  
  # Extracting the attributes (if needed)
  if (is.list(x) && network::is.network(x[[1]])) {
    
    # Vertex attributes
    VATTRS <- lapply(x, network::list.vertex.attributes)
    VATTRS <- lapply(VATTRS, setdiff, y = "na") # Known problem in network
    VATTRS <- mapply(
      function(net, a)
        as.data.frame(structure(
          lapply(a, network::get.vertex.attribute, x = net),
          names = a
          )),
      net = x, a = VATTRS, SIMPLIFY = FALSE
      )
    
    VATTRS <- do.call(rbind, VATTRS)
    if (ncol(VATTRS))
      VATTRS <- as.list(VATTRS)
    else
      VATTRS <- NULL
    
    # Edge attributes
    EATTRS <- lapply(x, network::list.edge.attributes)
    EATTRS <- lapply(EATTRS, setdiff, y = "na") # Known problem in network
    EATTRS <- mapply(
      function(net, a) 
        structure(
          lapply(a, network::get.edge.attribute, x = net),
          names = a
          ),
      net = x, a = EATTRS
    )
    
    EATTRS <- do.call(rbind, EATTRS)
    if (ncol(EATTRS))
      EATTRS <- as.list(EATTRS)
    else
      EATTRS <- NULL
    
  } else {
    VATTRS <- NULL
    EATTRS <- NULL
  }
  
  ADJMAT <- matrix(0, nrow = N, ncol = N)
  IDX    <- cumsum(sizes)
  
  ADJMAT_LIST  <- as_adjmat(x)
  VATTRS[[attrname]] <- vector("integer", N)
  for (i in seq_along(IDX)) {
    block <- (IDX[i] + 1L - sizes[i]):IDX[i]
    VATTRS[[attrname]][block] <- i
    ADJMAT[block,][,block] <- ADJMAT_LIST[[i]]
  }
  
  # Creating the network object
  NET <- network::as.network(ADJMAT)
  
  # Adding attributes
  for (n in names(VATTRS)) {
    
    if (n == "name")
      network::set.vertex.attribute(NET, "name_original", as.vector(VATTRS[[n]]))
    else
      network::set.vertex.attribute(NET, n, as.vector(VATTRS[[n]]))
    
  }
  
  for (n in names(EATTRS)) {
    
    if (n == "name")
      network::set.edge.attribute(NET, "name_original", EATTRS[n])
    else
      network::set.edge.attribute(NET, n, EATTRS[n])
    
  }
  
  return(NET)
  
}

#' @export
#' @rdname blockdiagonalize
splitnetwork <- function(x, attrname) {
  
  vertex_block <- network::get.vertex.attribute(x, attrname)
  blocks       <- sort(unique(vertex_block))

  L <- vector("list", length(blocks))
  for (i in seq_along(blocks)) {
    ids <- which(vertex_block == blocks[i])
    L[[i]] <- network::get.inducedSubgraph(x, ids)
  }
  
  L
  
}

# library(ergmito)
# data("fivenets")
# blockdiagonalize(fivenets)
# (ans <- ergm_blockdiag(fivenets ~ mutual))
# ans$stats.weights
