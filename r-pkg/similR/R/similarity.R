#' Similarity and Distance between pairs of binary matrices
#' 
#' Multiple measurements of similarity and distance between pairs of binary
#' matrices as listed in Choi et al (2010).
#' 
#' 
#' @param statistic Character. Name of the similarity index to be using.
#' @param normalized Logical. When `TRUE` it returns the normalized hamming distance,
#' which ranges between 0 and 1 (currently only used in `statistic="hamming"`).
#' @param firstonly Logical. When `TRUE`, the comparison is done as the first
#' matrix to all only.
#' @param exclude_j Logical. When `TRUE`, the comparison between matrices `i` and
#' `j` is done after the jth column and rows are removed from each.
#' @template matrix
#' @details 
#' All of the available statistics are based on a 2x2 contingency matrix counting
#' matches and missmatches between each pair of matrices (R1, R2).
#'
#' \tabular{llcc}{
#'      \tab   \tab R2 \tab   \cr
#'      \tab   \tab 1  \tab 0 \cr
#'   R1 \tab 1 \tab a  \tab b \cr
#'      \tab 0 \tab c  \tab d  
#' }
#'
#' A complete list of the statistics available is available in the **Similarity**
#' and **Distance** sections.
#' 
#' `distance` is just an alias for `similarity`.
#' 
#' @export
#' @return A matrix of size `n*(n - 1)/2` by `length(statistic)`, where columns
#' 1 and 2 indicate the id if `i` and `j`, and the reminder columns are the
#' corresponding distances/similarities. 
#' 
#' @references 
#' 
#' Choi, S. S., Cha, S. H., & Tappert, C. C. (2010). A survey of binary similarity
#' and distance measures. Journal of Systemics, Cybernetics and Informatics, 8(1), 43-48.
#' 
#' Krackhardt, D. (1990). Assessing the political landscape: Structure, cognition,
#' and power in organizations. Administrative science quarterly, 342-369.
#' 
#' Gower, J. C., & Legendre, P. (1986). Metric and Euclidean properties of
#' dissimilarity coefficients. Journal of classification, 3(1), 5-48.
#' 
#' @examples 
#' # Getting all networks of size 3
#' data(powerset03)
#'  
#' # We can compute it over the entire set
#' head(similarity(powerset03, statistic="s14"))
#' 
#' # Or over two pairs
#' head(similarity(powerset03[[1]], powerset03[[2]], powerset03[[3]], statistic="s14"))
#' 
#' # We can compute multiple distances at the same time
#' ans <- similarity(powerset03, statistic=c("hamming", "dennis", "jaccard"))
#' head(ans)
#' 
#' @seealso The [statistics] object contains a list with the available statistics
#' for convenience.
#' 
similarity <- function(M, ..., statistic, normalized = TRUE, firstonly=FALSE, exclude_j = FALSE)
  UseMethod("similarity")


#' @export
#' @rdname similarity
similarity.list <- function(M, ..., statistic, normalized = TRUE, firstonly=FALSE, exclude_j = FALSE) {
  
  # Checking the dots
  dots <- list(...)
  if (length(dots))
    stop(
      "If `M` is a list, all arguments (matrices) should be provided via `M` and not `...`.",
      call. = FALSE
      )
  
  if (length(statistic) > 1) {
    ans <- lapply(statistic, .similarity, M = M, normalized = normalized,
                  firstonly = firstonly, exclude_j = exclude_j)
    ans <- do.call(
      cbind,
      c(ans[1], lapply(ans[-1], "[", i=, j=-c(1,2)))
      )
  } else {
    ans <- .similarity(M = M, statistic = statistic, normalized = normalized,
                       firstonly = firstonly, exclude_j = exclude_j)
  }
  
  colnames(ans) <- c("i", "j", statistic)
  
  ans
  
}

#' @export
#' @rdname similarity
similarity.matrix <- function(M, ..., statistic, normalized=TRUE, firstonly=FALSE, exclude_j = FALSE) {
  
  dots <- list(...)
  if (!length(dots))
    stop(
      "Only one matrix was provided. `similarity` can only be computed between 2 matrices.",
      call. = FALSE
      )
  
  similarity(M = c(list(M), dots), statistic=statistic, normalized=normalized,
             firstonly = firstonly, exclude_j = exclude_j)
  
}
