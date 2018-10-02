#' Power set of Directed Graphs of size `n`
#' @param n Integer. Number of edges.
#' @param force Logical. When `TRUE` it generates the powerset for `n>5`, otherwise
#' it returns with error.
#' @param as_matrix Logical. When `TRUE` it returns a list of matrices of size `n`,
#' otherwise it returns a list of vectors with the not-zero locations of the
#' matrix in column-major form.
#' @param mc_cores Integer. Passed to [parallel::mclapply]
#' @examples 
#' powerset(2)
#' @export
powerset <- function(n, force=FALSE, as_matrix=TRUE, mc_cores=2L) {
  
  s <- .powerset(n, force)
  m0 <- matrix(0, ncol=n, nrow=n)
  if (as_matrix) {
    return(parallel::mclapply(s, function(m) {
      m0[m + 1L] <- 1
      m0
      }, mc.cores = mc_cores))
  }
  
  parallel::mclapply(s, "+", 1, mc.cores=mc_cores)
  
}
