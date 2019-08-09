
#' @importFrom Rcpp sourceCpp
#' @useDynLib ergmito, .registration = TRUE
#' @importFrom ergm ergm.allstats ergm
#' @importFrom network network.size is.network get.vertex.attribute
#' list.vertex.attributes
#' @importFrom stats pnorm model.matrix update.formula optim AIC BIC rnorm terms
#'  coef coefficients var quantile aggregate as.formula printCoefmat formula
#'  simulate
#' @importFrom methods setClass setGeneric setMethod setRefClass
#' @importFrom utils capture.output
NULL
