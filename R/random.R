#' Random bernoulli graph
#' @param n Integer. Size of the graph
#' @param p Probability of a tie
#' @return A square matrix of size `n` with zeros in the diagonal.
#' @examples 
#' # A graph of size 4
#' rbernoulli(4)
#' @export
#' @importFrom stats runif
rbernoulli <- function(n, p=.5) {
  
  net <- matrix(as.integer(stats::runif(n^2) > p), ncol=n)
  diag(net) <- 0
  net
  
}