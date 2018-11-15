#' Locally Aggregated Structures (LAS)
#' 
#' Generates LAS matrices using different criteria as described in Krackhardt (1990).
#' 
#' @template matrix
#' @param rule Character. Rule used to create the LAS matrix. Can be
#' any of `"intersection"`, `"union"`, or `"threshold"` (see the *Rule* section).
#' @param threshold Number of minimum agreements to set the tie to 1. This is used
#' for the case when `rule = "threshold"`.
#' @references
#' Krackhardt, D. (1987). Cognitive social structures. Social networks, 9(2), 109-134.
#' @aliases LAS CS Concensus-Structure
#' @examples 
#' data(powerset03)
#' x <- powerset03[1:3]
#' 
#' # It only exists if i and j declare it
#' las(x, rule="i")
#' 
#' # Same as above
#' las(x[[1]], x[[2]], x[[3]], rule="i")
#' 
#' # It exists if either i or j declare it
#' las(x, rule="u")
#' 
#' # Not enough agreement
#' las(x, rule="t", threshold=.9)
#' @export
las <- function(M, ..., rule, threshold=.5) UseMethod("las")

#' @export
#' @rdname las
LAS <- las

#' @export
#' @rdname las
las.list <- function(M, ..., rule, threshold=.5) {
  
  # Checking the dots
  dots <- list(...)
  if (length(dots))
    stop(
      "If `M` is a list, all arguments (matrices) should be provided via `M` and not `...`.",
      call. = FALSE
    )
  
  .las(M = M, rule = rule, threshold = threshold)
  
}

#' @export
#' @rdname las
las.matrix <- function(M, ..., rule, threshold=.5) {
  
  dots <- list(...)
  if (!length(dots))
    stop(
      "Only one matrix was provided. `las` can only be computed between 2 matrices.",
      call. = FALSE
    )
  
  las(M = c(list(M), dots), rule=rule, threshold=threshold)
  
}