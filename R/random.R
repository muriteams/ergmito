#' Random Bernoulli graph
#' 
#' @param n Integer vector. Size of the graph. If `length(n) > 1`, then it will
#' a list of random graphs.
#' @param p Probability of a tie. This may be either a scalar, or a vector of the
#' same length of `n`.
#' @return If `n` is a single number, a square matrix of size `n` with zeros in
#' the diagonal. Otherwise it returns a list of `length(n)` square matrices of
#' sizes equal to those specified in `n`.
#'  
#' @examples 
#' # A graph of size 4
#' rbernoulli(4)
#' 
#' # 3 graphs of various sizes
#' rbernoulli(c(3, 4, 2))
#' 
#' # 3 graphs of various sizes and different probabilities
#' rbernoulli(c(3, 4, 6), c(.1, .2, .3))
#' @export
#' @importFrom stats runif
rbernoulli <- function(n, p=.5) {
  
  if (length(n) > 1)
    return(Map(rbernoulli, n = n, p = p))
  
  net <- matrix(as.integer(stats::runif(n^2) < p), ncol=n)
  diag(net) <- 0
  net
  
}