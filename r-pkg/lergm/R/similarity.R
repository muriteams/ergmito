#' Similarity indexes
#' 
#' Includes the S14 statistic and hamming distance between two or more
#' pairs of matrices.
#' 
#' @noMd
#' @param x Either a list of matrices of size `n` (need not to be square), or
#' a single matrix of size `n` (see details).
#' @param bool Logical. When `TRUE` it returns the normalized hamming distance,
#' which ranges between 0 and 1.
#' @param ... More matrices to be passed to the function.
#' @details 
#' If `x` is a matrix, then the function requires the user to pass at least a
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
#' head(S14(nets))
#' 
#' # Or over two pairs
#' S14(nets[[1]], nets[[2]], nets[[3]])
#' 
#' # We can do the same with the hamming distance
#' hamming(nets)
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

#' @rdname S14
#' @export
#' @details The hamming distance is computed as the total number of entries that
#' are different between a single pair of individuals. When `normalized = TRUE`
#' such count is divided by `n*(n-1)` so it ranges between 0 and 1.
hamming <- function(x, ..., normalized=FALSE) UseMethod("hamming")

#' @rdname S14
#' @export
hamming.matrix <- function(x, ..., normalized=FALSE) {
  
  dots <- list(...)
  if (!length(dots))
    stop(
      "Only one matrix was provided. S14 can only be computed between 2 matrices.",
      call. = FALSE
    )
  
  hamming(c(list(x), dots), normalized = normalized)
  
}

#' @rdname S14
#' @export
hamming.list <- function(x, ..., normalized = FALSE) {
  
  # Checking the dots
  dots <- list(...)
  if (length(dots))
    stop(
      "If `x` is a list, all arguments (matrices) should be provided via `x` and not `...`.",
      call. = FALSE
    )
  
  .hamming(x, normalized = normalized)
  
}

#' @rdname S14
#' @export
starwid <- function(x, ..., normalized=FALSE) UseMethod("starwid")

#' @rdname S14
#' @export
starwid.matrix <- function(x, ..., normalized=FALSE) {
  
  dots <- list(...)
  if (!length(dots))
    stop(
      "Only one matrix was provided. S14 can only be computed between 2 matrices.",
      call. = FALSE
    )
  
  starwid(c(list(x), dots), normalized = normalized)
  
}

#' @rdname S14
#' @export
starwid.list <- function(x, ..., normalized = FALSE) {
  
  # Checking the dots
  dots <- list(...)
  if (length(dots))
    stop(
      "If `x` is a list, all arguments (matrices) should be provided via `x` and not `...`.",
      call. = FALSE
    )
  
  .starwid(x, normalized = normalized)
  
}

#' @rdname S14
#' @export
dennis <- function(x, ..., normalized=FALSE) UseMethod("dennis")

#' @rdname S14
#' @export
dennis <- function(x, ..., normalized=FALSE) {
  
  dots <- list(...)
  if (!length(dots))
    stop(
      "Only one matrix was provided. S14 can only be computed between 2 matrices.",
      call. = FALSE
    )
  
  dennis(c(list(x), dots), normalized = normalized)
  
}

#' @rdname S14
#' @export
dennis <- function(x, ..., normalized = FALSE) {
  
  # Checking the dots
  dots <- list(...)
  if (length(dots))
    stop(
      "If `x` is a list, all arguments (matrices) should be provided via `x` and not `...`.",
      call. = FALSE
    )
  
  .dennis(x, normalized = normalized)
  
}
