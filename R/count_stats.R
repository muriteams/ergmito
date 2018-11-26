#' Count Network Statistics
#' 
#' This function is similar to what [ergm::summary_formula] does, but it provides
#' a fast wrapper suited for matrix class objects.
#' @param X List of square matrices. (networks)
#' @param terms Character vector with the names of the statistics to calculate.
#' Currently, the only available statistics are: `edges` and `mutual`.
#' @export
#' @return A matrix of size `length(X) * length(terms)` with the corresponding
#' counts of statistics.
#' @examples 
#' # DGP 
#' x <- powerset(5)
#' ans0 <- count_stats(x[1:20], c("mutual", "edges"))
#' 
#' # Calculating using summary_formula
#' fm <- x[[i]] ~ mutual + edges
#' ans1 <- lapply(1:20, function(i) {
#'   environment(fm) <- environment()
#'   ergm::summary_formula(fm)
#' })
#' 
#' ans1 <- do.call(rbind, ans1)
#' 
#' # Comparing
#' all.equal(unname(ans0), unname(ans1))
#' 
count_stats <- function(X, terms) {
  
  chunks <- make_chunks(length(X), 2e5)
  
  ans <- matrix(NA, nrow = length(X), ncol=length(terms))
  
  for (s in seq_along(chunks$from)) {
    
    i <- chunks$from[s]
    j <- chunks$to[s]
    
    ans[i:j,] <- count_stats.(X[i:j], terms)
    
  }
  
  ans
  
}
