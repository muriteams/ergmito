#' S14 similarity index
#' @noMd
#' @param x Either a list of matrices of size `n` (need not to be square), or
#' a single matrix of size `n` (see details).
#' @details 
#' If `x` is a matrix, then the function requires the user to pass at least a
#' second matrix via de `...` notation.
#' Following Krackhardt's 1989:
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
#' head(S14(nets))
#' 
#' # Or over a couple of pairs
#' S14(nets[[1]], nets[[2]])
#' 
S14 <- function(x, ...) UseMethod("S14")

#' @export
#' @rdname S14
S14.list <- function(x, ...) {
  
  # Checking the dots
  dots <- list(...)
  if (length(dots))
    stop(
      "If `x` is a list, all arguments (matrices) should be provided via `x` and not `...`.",
      call. = FALSE
      )
  
  .S14(x)
  
}

#' @export
#' @rdname S14
S14.matrix <- function(x, ...) {
  
  dots <- list(...)
  if (!length(dots))
    stop(
      "Only one matrix was provided. S14 can only be computed between 2 matrices.",
      call. = FALSE
      )
  
  S14(c(list(x), dots))
  
}