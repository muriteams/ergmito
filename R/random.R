#' Random bernoulli graph
#' @param n Integer vector. Size of the graph. If `length(n) > 1`, then it will
#' a list of random graphs.
#' @param p Probability of a tie
#' @return A square matrix of size `n` with zeros in the diagonal.
#' @examples 
#' # A graph of size 4
#' rbernoulli(4)
#' 
#' # 3 graphs of various sizes
#' rbernoulli(c(3, 4, 2))
#' @export
#' @importFrom stats runif
rbernoulli <- function(n, p=.5) {
  
  if (length(n) > 1)
    return(Map(rbernoulli, n = n, p = p))
  
  net <- matrix(as.integer(stats::runif(n^2) > p), ncol=n)
  diag(net) <- 0
  net
  
}