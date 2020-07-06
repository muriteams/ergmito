# Methods ----------------------------------------------------------------------

# @rdname ergmito
#' @export
print.ergmito <- function(x, ...) {
  
  cat("\nERGMito estimates\n")
  if (length(x$note))
    cat(sprintf("note: %s\n", x$note))
  
  print(structure(unclass(x), class="ergm"))
  
  invisible(x)
  
}

# @rdname ergmito
#' @export
summary.ergmito <- function(object, ...) {
  
  # Computing values
  sdval <- sqrt(diag(vcov(object)))
  z     <- coef(object)/sdval
  
  is_boot <- inherits(object, "ergmito_boot")
  
  # Generating table
  ans <- structure(
    list(
      coefs = data.frame(
        Estimate     = coef(object),
        `Std. Error` = sdval,
        `z value`    = z,
        `Pr(>|z|)`   = 2*stats::pnorm(-abs(z)),
        row.names    = names(coef(object)),
        check.names  = FALSE
      ),
      aic         = stats::AIC(object),
      bic         = stats::BIC(object),
      model       = object$formulae$model_final,
      note        = object$note,
      R           = ifelse(is_boot, object$R, 1L),
      n           = nnets(object)
    ),
    class = c("ergmito_summary", if (is_boot) "ergmito_summary_boot" else  NULL)
  )
  
  ans
}

# @rdname ergmito
#' @export
print.ergmito_summary <- function(
  x,
  ...
) {
  
  cat("\nERGMito estimates (MLE)\n")
  cat("This model includes", x$n, "networks.\n")
  if (x$R > 1L)
    cat("\n(bootstrapped model with ", x$R, " replicates.)\n")
  
  if (length(x$note))
    cat(sprintf("note: %s\n", x$note))
  
  cat("\nformula:\n  ")
  print(x$model)
  cat("\n")
  
  stats::printCoefmat(
    x$coefs,
    signif.stars  = TRUE,
    signif.legend = TRUE
  )
  
  cat(paste("AIC:", format(x$aic), 
            "  ", "BIC:", format(x$bic), 
            "  ", "(Smaller is better.)", "\n", sep = " "))
  
  invisible(x)
  
}

# @rdname ergmito
#' @export
#' @importFrom stats coef logLik vcov nobs formula
coef.ergmito <- function(object, ...) {
  
  object$coef
  
}

# @rdname ergmito
#' @export
logLik.ergmito <- function(object, ...) {
  
  object$mle.lik
  
}

# @rdname ergmito
#' @export
nobs.ergmito <- function(object, ...) {
  
  object$nobs
  
}

#' @export
#' @param solver Function. Used to compute the inverse of the hessian matrix. When
#' not null, the variance-covariance matrix is recomputed using that function.
#' By default, `ergmito` uses [MASS::ginv].
#' @rdname ergmito
vcov.ergmito <- function(object, solver = NULL, ...) {
  
  if (is.null(solver))
    return(object$covar)
  
  structure(
    - solver(object$optim.out$hessian),
    dimnames = dimnames(object$covar)
  )
  
}

# @rdname ergmito
#' @export
formula.ergmito <- function(x, ...) {
  
  x$formulae$model_final
  
}
