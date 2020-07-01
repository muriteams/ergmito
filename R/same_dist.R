same_dist_wrap <- function(x, what, map01 = NULL, map10 = NULL) {
  structure(
    x,
    what  = what,
    class = "ergmito_same_dist",
    map01 = map01,
    map10 = map10
    )
}

`%notin%` <- function(x, table) any(!(x %in% table))

#' Compare pairs of networks to see if those came from the same distribution
#' 
#' If two networks are of the same size, and their vertex attributes are
#' equal in terms of set comparison, then we say those came from the same
#' distribution
#' 
#' @param net0,net1 Networks to be compared.
#' @param attrnames Character vector. (optional) Names of the vertex attributes
#' to be be compared on. This is ignored in the matrix case.
#' @param ... Ignored.
#' 
#' @details This function is used during the call of [ergmito_formulae] to
#' check whether the function can recycle previously computed statistics for
#' the likelihood function. In the case of models that only contain structural
#' terms, i.e. attribute less, this can save significant amount of computing
#' time and memory.
#' 
#' @return A logical with an attribute `what`. `TRUE` meaning that the two
#' networks come from the same distribution, and `FALSE` meaning that they do
#' not. If `FALSE` the `what` attribute will be equal to either `"size"` or
#' the name of the attribute that failed the comparison.
#' 
#' @examples 
#' data(fivenets)
#' same_dist(fivenets[[1]], fivenets[[2]]) # Yes, same size
#' same_dist(fivenets[[1]], fivenets[[2]], "female") # No, different attr dist
#' @export
same_dist <- function(net0, net1, attrnames = NULL, ...) UseMethod("same_dist")

#' @export
# @rdname same_dist
same_dist.matrix <- function(net0, net1, attrnames = NULL, ...) {
  
  if (!inherits(net1, "matrix"))
    stop(
      "Cannot compare two objects with differnt classes. ",
      "`net0` is of class matrix while `net1` of class '",
      class(net1), "'.",
      call. = FALSE
    )
  
  if (nvertex(net0) != nvertex(net1))
    return(same_dist_wrap(FALSE, "size"))
  
  return(same_dist_wrap(TRUE, NULL))
  
}

#' @export
# @rdname same_dist
same_dist.network <- function(net0, net1, attrnames = NULL, ...) {
  
  if (!network::is.network(net1))
    stop(
      "Cannot compare two objects with differnt classes. ",
      "`net0` is of class network while `net1` of class '",
      class(net1), "'.",
      call. = FALSE
      )

  # 1. Size
  if (nvertex(net0) != nvertex(net1))
    return(same_dist_wrap(FALSE, "size"))
  
  # 2. Attributes
  if (length(attrnames) && !is.null(attrnames)) {
    
    # First question: Are the requested attributes present?
    anames0 <- network::list.vertex.attributes(net0)
    if (attrnames %notin% anames0)
      stop(
        "One or more attributes listed on `attrnames` ",
        "are not preset in `net0`.", call. = FALSE
        )
    
    anames1 <- network::list.vertex.attributes(net0)
    if (attrnames %notin% anames1)
      stop(
        "One or more attributes listed on `attrnames` ",
        "are not preset in `net1`.", call. = FALSE
      )
    
    # Now we compare, joint attribute distribution. While may not be necessary
    # all the time, it is safer to do so.
    a0 <- matrix(nrow = nvertex(net0), ncol = length(attrnames))
    a1 <- matrix(nrow = nvertex(net0), ncol = length(attrnames))
    for (i in seq_along(attrnames)) {
      
      a0[,i] <- network::get.vertex.attribute(net0, attrnames[i])
      a1[,i] <- network::get.vertex.attribute(net1, attrnames[i])
      
    }
    
    # Sorting the arguments accordingly
    ord0 <- do.call("order", lapply(1:length(attrnames), function(i) a0[,i]))
    ord1 <- do.call("order", lapply(1:length(attrnames), function(i) a1[,i]))

    if (any(a0[ord0,] != a1[ord1,]))
      return(same_dist_wrap(FALSE, attrnames))
    
    # Returning permutations, this can be useful if one wants to do something
    # for one and map it ot the other... this was hard to figure out!
    n <- nvertex(net0)
    ord0_ <- order(ord0)
    ord1_ <- order(ord1)
    return(
      same_dist_wrap(
        TRUE, attrnames,
        map01 = ord0[ord1_],
        map10 = ord1[ord0_]
        )
    )
    
    
  }
  
  same_dist_wrap(TRUE, attrnames)
    
}
