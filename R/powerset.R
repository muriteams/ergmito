
#' Power set of Graphs of size `n`
#' 
#' Generates the set of all possible binary networks of size `n`.
#' 
#' @param n Integer. Number of edges.
#' @param force Logical scalar. When `TRUE` it generates the power set for `n>5`,
#' otherwise it returns with error.
#' @param directed Logical scalar. Whether to generate the power set of directed
#' or undirected graphs,
#' @param chunk_size Number of matrices to process at a time. If n = 5, then
#' stack memory on the computer may overflow if `chunk_size` is relatively large.
#' @examples 
#' powerset(2)
#' powerset(4, directed = FALSE)
#' @export
powerset <- function(
  n,
  directed   = TRUE,
  force      = FALSE,
  chunk_size = 2e5
  ) {
  
  # Calculating power sets
  sets <- .powerset(n, force, directed)
  
  N <- ifelse(directed, 2^(n*(n-1)), 2^(n*(n-1)/2))
  chunks <- make_chunks(N, chunk_size = chunk_size)
  
  ans <- vector("list", N)
  for (s in seq_along(chunks$from)) {

    i <- chunks$from[s]
    j <- chunks$to[s]
    
    ans[i:j] <- wrap_powerset(sets, from=i-1, to=j, n=n)
    
  }
  
  # Reflecting the matrix
  if (!directed)
    ans <- lapply(ans, function(x) x + t(x))
  
  ans
}
