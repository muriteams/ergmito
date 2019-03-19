
#' @importFrom Rcpp sourceCpp
#' @useDynLib ergmito, .registration = TRUE
#' @importFrom ergm ergm.allstats
#' @importFrom network network.size
#' @importFrom stats pnorm model.matrix update.formula optim AIC BIC rnorm terms
#'  coef coefficients var
#' @importFrom methods setClass setGeneric setMethod setRefClass
#' @importFrom utils capture.output
NULL