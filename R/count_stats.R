#' Count Network Statistics
#' 
#' This function is similar to what [ergm::summary_formula] does, but it provides
#' a fast wrapper suited for matrix class objects (see benchmark in the examples).
#' @param X List of square matrices. (networks)
#' @param terms Character vector with the names of the statistics to calculate.
#' Currently, the only available statistics are: '\Sexpr{paste(ergmito::AVAILABLE_STATS(), collapse="', '")}'.
#' @param ... Passed to the method.
#' @param attrs A list of vectors. This is used when `term` has a nodal attribute
#' such as `nodeicov(attrname="")`.
#' @export
#' @return A matrix of size `length(X) * length(terms)` with the corresponding
#' counts of statistics.
#' @examples 
#' # DGP 
#' set.seed(123199)
#' x <- rbernoulli(rep(5, 10))
#' ans0 <- count_stats(x, c("mutual", "edges"))
#' 
#' # Calculating using summary_formula
#' ans1 <- lapply(x, function(i) {
#'   ergm::summary_formula(i ~ mutual + edges)
#' })
#' 
#' ans1 <- do.call(rbind, ans1)
#' 
#' # Comparing
#' all.equal(unname(ans0), unname(ans1))
#' 
#' # count_stats is vectorized (and so faster)
#' bm <- benchmarkito(
#'   count_stats = count_stats(x, c("mutual", "edges")),
#'   lapply      = lapply(x, function(i) {
#'   ergm::summary_formula(i ~ mutual + edges)
#' }), times = 50
#' )
#' 
#' plot(bm)
#' 
count_stats <- function(X, ...) UseMethod("count_stats")

#' @export
#' @rdname count_stats
AVAILABLE_STATS <- function() count_available()


#' @export
#' @rdname count_stats
count_stats.formula <- function(X, ...) {
  
  # Retrieving networks
  LHS <- eval(X[[2]], envir = environment(X))
  
  if (inherits(LHS, "matrix") | inherits(LHS, "network"))
    LHS <- list(LHS)

  # Checking which ones are undirected
  are_undirected <- !is_directed(LHS)
  are_undirected <- which(are_undirected)
  
  if (length(are_undirected))
    stop(
      "Counting statistics with count_stats in undirected networks is not ",
      "supported yet. The following networks in the formula are undirected: ",
      paste(are_undirected, collapse = ", "), ".", call. = FALSE
      )
    
  # Analyzing the formula (with network as a reference)
  ergm_model <- analyze_formula(X, LHS)
  
  # Can we do it?
  available <- which(!(ergm_model$term_names %in% count_available()))
  if (length(available))
    stop(
      "The following term(s)s are not available in count_stats: ",
      paste(ergm_model$names[available], collapse = ", "),
      ".", call. = FALSE
      )
  
  # Capturing attributes
  for (a in seq_along(ergm_model$term_attrs)) {
    ergm_model$attrs[[a]] <- if (!length(ergm_model$term_attrs[[a]]))
      double(0)
    else {
      
      # This check is important, for now. Future versions may include more
      # complex terms that hold more than one attribute.
      if (length(ergm_model$term_attrs[[a]]) > 1L)
        stop(
          "For now, terms with more than one attribute are not supported on. ",
          "The current model you are trying to fit uses the term: ",
          ergm_model$term_passed[a],
          " which includes the following attributes: ",
          paste(ergm_model$term_attrs[[a]], collapse=", "),
          call. = FALSE
          )
        
      lapply(LHS, function(net) {
        network::get.vertex.attribute(net, attrname = ergm_model$term_attrs[[a]])
      })
      
    }
  }
  
  # Coercion is later since we need to check for arguments
  LHS <- as_adjmat(LHS)
  
  # Coercing into the appropiate type
  if (network::is.network(LHS))
    LHS <- list(as_adjmat(LHS))
  else if (is.list(LHS)) {
    
    is_net <- sapply(LHS, network::is.network)
    
    # Coercing into a net
    for (i in which(is_net))
      LHS[[i]] <- as_adjmat(LHS[[i]])
    
  }
  
  out <- matrix(nrow = nnets(LHS), ncol = length(ergm_model$term_names),
                dimnames = list(NULL, ergm_model$term_passed))
  
  for (j in 1:ncol(out)) {
    
    out[, j] <- count_stats(
      X     = LHS,
      terms = ergm_model$term_names[j],
      attrs = ergm_model$attrs[[j]]
      )
    
  }
  
  out

  
}

#' @export
#' @rdname count_stats
count_stats.list <- function(X, terms, attrs = NULL, ...) {
  
  chunks <- make_chunks(length(X), 2e5)
  
  # if (!length(attrs))
  #   attrs <- replicate(length(X), numeric(0), simplify = FALSE)
  
  # Checking the types of objects
  test <- which(!sapply(X, inherits, what = "matrix"))
  if (length(test))
    stop("When `X` is a list, it must be a list of matrices. There are ",
         "some objects that are not: ", paste(test, collapse = ", "), ".",
         call. = FALSE)
  
  ans <- matrix(NA, nrow = length(X), ncol=length(terms))
  all_same_attr <- length(attrs) == 1L
  for (s in seq_along(chunks$from)) {
    
    i <- chunks$from[s]
    j <- chunks$to[s]
    
    for (k in seq_along(terms)) {
      if (!length(attrs))
        ans[i:j, k] <- count_stats.(X[i:j], terms[k], list(double(0L)))
      else if (all_same_attr)
        ans[i:j, k] <- count_stats.(X[i:j], terms[k], attrs)
      else
        ans[i:j, k] <- count_stats.(X[i:j], terms[k], attrs[i:j])
    }
    
  }
  
  ans
  
}

#' Geodesic distance matrix (all pairs)
#' 
#' Calculates the shortest path between all pairs of vertices in a network.
#' This uses the power matrices to do so, which makes it efficient only for
#' small networks.
#' 
#' @param x Either a list of networks (or square integer matrices), an integer
#' matrix, a network, or an ergmito.
#' @param force Logical scalar. If `force = FALSE` (the default) and `nvertex(x) > 100`
#' it returns with an error. To force computation use `force = TRUE`.
#' @param ... Further arguments passed to the method.
#' @param simplify Logical scalar. When `TRUE` it returns a matrix, otherwise, 
#' a list of length `nnets(x)`.
#' 
#' @export
#' @examples 
#' data(fivenets)
#' geodesic(fivenets)
#' 
#' # Comparing with sna
#' if (require("sna")) {
#'   net0 <- fivenets[[1]]
#'   net  <- network::network(fivenets[[1]])
#'   benchmarkito(
#'     ergmito = ergmito::geodesic(net0),
#'     sna     = sna::geodist(net), times = 1000
#'   )
#' }
geodesic <- function(x, force = FALSE, ...) UseMethod("geodesic")

#' @export
#' @rdname geodesic
geodesita <- geodesic

#' @export
# @rdname geodesic
geodesic.list <- function(x, force = FALSE, ...) {
  
  geodesic.(as_adjmat(x), force = force)
  
}

#' @export
#' @rdname geodesic
geodesic.matrix <- function(x, force = FALSE, simplify = FALSE, ...) {
  
  ans <- geodesic.(list(x), force = force)
  if (simplify)
    return(ans[[1]])
  ans
}

#' @export
#' @rdname geodesic
geodesic.network <- function(x, force = FALSE, simplify = FALSE, ...) {
  
  ans <- geodesic.(list(as_adjmat(x)), force = force)
  if (simplify)
    return(ans[[1]])
  
  ans
}
