#' Extract function to be used with the `texreg` package.
#'
#' To be used with the \CRANpkg{texreg} package. This function can be used to
#' generate nice looking tables of ERGMitos estimates.
#' 
#' @param model An object of class `ergmito`.
#' @param include.aic,include.bic,include.loglik See [texreg::extract].
#' @param include.nnets Logical. When true, it adds the Number of networks used
#' to the list of gof statistics. This can be useful when running pooled models.
#' @param include.convergence Logical. When true it, adds the convergence value
#' of the [stats::optim] function (0 means convergence).
#' @param ... Further arguments passed to [summary.ergmito].
#' @export
#' @examples 
#' 
#' library(texreg)
#' data(fivenets)
#' ans <- ergmito(fivenets ~ edges + nodematch("female"))
#' screenreg(ans)
#' 
#' @importFrom texreg extract
extract.ergmito <- function(
  model, include.aic = TRUE, include.bic = TRUE, include.loglik = TRUE,
  include.nnets = TRUE,
  include.convergence = TRUE,
  ...
) {
  
  if (!requireNamespace("texreg", quietly = TRUE))
    stop("Need to install the `texreg` package.", call. = FALSE)
  
  
  # Copied from texreg::extract.ergm
  
  s <- summary(model, ...)
  coefficient.names <- rownames(s$coefs)
  coefficients <- s$coefs[, 1]
  standard.errors <- s$coefs[, "Std. Error"]
  significance <- s$coefs[, "Pr(>|z|)"]
  
  # Adding gof statistics ------------------------------------------------------
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE && !is.null(s$aic)) {
    aic <- s$aic
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE && !is.null(s$bic)) {
    bic <- s$bic
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE && !is.null(model$mle.lik)) {
    lik <- model$mle.lik[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  
  # Number of networks 
  if (include.nnets) {
    gof         <- c(gof, nnets(model$network))
    gof.names   <- c(gof.names, "Num. networks")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  # Convergence
  if (include.convergence) {
    gof         <- c(gof, model$optim.out$convergence)
    gof.names   <- c(gof.names, "Convergence")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- texreg::createTexreg(coef.names = coefficient.names, coef = coefficients,
                             se = standard.errors, pvalues = significance, gof.names = gof.names,
                             gof = gof, gof.decimal = gof.decimal)
  return(tr)
  
}


setMethod(
  "extract", signature = className("ergmito", "ergmito"),
  definition = extract.ergmito
)
