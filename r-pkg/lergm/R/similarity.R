#' Similarity indexes
#' 
#' Includes the S14 statistic and hamming distance between two or more
#' pairs of matrices.
#' 
#' @noMd
#' @param M Either a list of matrices of size `n` (need not to be square), or
#' a single matrix of size `n` (see details).
#' @param statistic Character. Name of the similarity index to be using.
#' @param bool Logical. When `TRUE` it returns the normalized hamming distance,
#' which ranges between 0 and 1.
#' @param ... More matrices to be passed to the function.
#' @details 
#' If `M` is a matrix, then the function requires the user to pass at least a
#' second matrix via de `...` notation.
#' 
#' In the case of the `S14` function, following Krackhardt's 1989:
#' \tabular{lcc}{
#'      \tab R2 \tab   \cr
#'   R1 \tab a  \tab b \cr
#'      \tab c  \tab d  
#' }
#' 
#' \deqn{%
#' \sqrt{\left(\frac{a}{(a + c)} - \frac{b}{(b + d)}\right)\times\left(\frac{a}{(a + b)} - \frac{c}{(c + d)}\right)}
#' }{%
#' S14 = [(a/(a + c) - b/(b + d))*(a/(a + b) - c/(c + d))]^(1/2)
#' }
#' 
#' Which is an statistic lying between 0 and 1.
#' @export
#' @return A matrix of size `n*(n - 1)/2` by 3, where columns 1 and 2 indicate
#' ego and alter, and the third column the `S14` statistic.
#' 
#' @references 
#' Krackhardt, "Assessing the Political Landscape: Structure, Cognition, and Power in Organizations",
#' Administrative Science Quarterly, Vol. 35, No. 2 (Jun., 1990), pp. 342-369
#' 
#' @examples 
#' # Getting all networks of size 3
#' nets <- powerset(3)
#' 
#' # We can compute it over the entire set
#' head(similarity(nets, statistic="S14"))
#' 
#' # Or over two pairs
#' similarity(nets[[1]], nets[[2]], nets[[3]], statistic="S14")
#' 
#' # We can do the same with the hamming distance
#' similarity(nets, statistic="hamming")
#' 
similarity <- function(M, ..., statistic, normalized = TRUE) UseMethod("similarity")

#' @export
#' @rdname similarity
similarity.list <- function(M, ..., statistic, normalized = TRUE) {
  
  # Checking the dots
  dots <- list(...)
  if (length(dots))
    stop(
      "If `M` is a list, all arguments (matrices) should be provided via `M` and not `...`.",
      call. = FALSE
      )
  
  if (length(statistic) > 1) {
    ans <- lapply(statistic, .similarity, M = M, normalized = normalized)
    ans <- do.call(
      cbind,
      c(ans[1], lapply(ans[-1], "[", i=, j=-c(1,2)))
      )
    colnames(ans) <- c("i", "j", statistic)
    ans
  }
  else
    .similarity(M = M, statistic = statistic, normalized = normalized)
  
}

#' @export
#' @rdname similarity
similarity.matrix <- function(M, ..., statistic, normalized=TRUE) {
  
  dots <- list(...)
  if (!length(dots))
    stop(
      "Only one matrix was provided. `similarity` can only be computed between 2 matrices.",
      call. = FALSE
      )
  
  similarity(M = c(list(x), dots), statistic=statistic, normalized=normalized)
  
}
