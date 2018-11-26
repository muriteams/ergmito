make_chunks <- function(N, chunk_size) {
  
  if (chunk_size > N)
    return(list(from=1, to=N))
  
  chunks <- seq(0, N, by = chunk_size)
  chunks <- list(from = chunks[-length(chunks)] + 1, to = chunks[-1])
  chunks$to[length(chunks$to)] <- N
  
  chunks
  
}


#' Power set of Directed Graphs of size `n`
#' @param n Integer. Number of edges.
#' @param force Logical. When `TRUE` it generates the powerset for `n>5`, otherwise
#' it returns with error.
#' @param chunk_size Number of matrices to process at a time. If n = 5, then
#' stack memory on the computer may overflow if `chunk_size` is relatively large.
#' @examples 
#' powerset(2)
#' @export
powerset <- function(
  n,
  force      = FALSE,
  chunk_size = 2e5
  ) {
  
  # Calculating power sets
  sets <- .powerset(n, force)
  
  N <- 2^(n*(n-1))
  chunks <- make_chunks(N, chunk_size = chunk_size)
  
  ans <- vector("list", N)
  for (s in seq_along(chunks$from)) {

    i <- chunks$from[s]
    j <- chunks$to[s]
    
    ans[i:j] <- wrap_powerset(sets, from=i-1, to=j, n=n)
    
  }
  
  ans
}
