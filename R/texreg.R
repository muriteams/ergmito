#' Extract function to be used with the `texreg` package.
#'
#' To be used with the \CRANpkg{texreg} package. This function can be used to
#' generate nice looking tables of ERGMitos estimates.
#' 
#' @param model An object of class `ergmito`.
#' @param include.aic,include.bic,include.loglik See [texreg::extract].
#' @param include.nnets Logical. When true, it adds the Number of networks used
#' to the list of gof statistics. This can be useful when running pooled models.
#' @param include.offset Logical. When equal to `TRUE`, it adds one line per
#' offset term to the table, omiting sd and significance.
#' @param include.convergence Logical. When true it, adds the convergence value
#' of the [stats::optim] function (0 means convergence).
#' @param include.timing Logical, When true it will report the elapsed time
#' in seconds.
#' @param ... Further arguments passed to the [base::summary()] of [ergmito].
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
  model,
  include.aic = TRUE,
  include.bic = TRUE,
  include.loglik = TRUE,
  include.nnets = TRUE,
  include.offset = TRUE,
  include.convergence = TRUE,
  include.timing      = TRUE,
  ...
) {
  
  if (!requireNamespace("texreg", quietly = TRUE))
    stop("Need to install the `texreg` package.", call. = FALSE)
  
  
  # Copied from texreg::extract.ergm
  
  s                 <- summary(model, ...)
  coefficient.names <- rownames(s$coefs)
  coefficients      <- s$coefs[, 1]
  standard.errors   <- s$coefs[, "Std. Error"]
  significance      <- s$coefs[, "Pr(>|z|)"]
  
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
  
  # Checking offset
  model_terms  <- terms(formula(model))
  if (include.offset && length(noffsets <- attr(model_terms, "offset"))) {
    
    offsets <- rownames(attr(model_terms, "factor"))[attr(model_terms, "offset")]
    offsets <- paste("(offset)", gsub("^offset[(]|[)]$", "", offsets))

    coefficient.names <- c(coefficient.names, offsets)
    coefficients      <- c(coefficients, rep(1, length(offsets)))
    standard.errors   <- c(standard.errors, rep(NA, length(offsets)))
    significance      <- c(significance, rep(NA, length(offsets)))
    
  }
  
  # Checking boot
  if (inherits(model, "ergmito_boot")) {
    
    # How many replicates
    gof         <- c(gof, model$R)
    gof.names   <- c(gof.names, "N replicates")
    gof.decimal <- c(gof.decimal, FALSE)
    
    # How many of R we used
    gof         <- c(gof, model$nvalid)
    gof.names   <- c(gof.names, "N Used replicates")
    gof.decimal <- c(gof.decimal, FALSE)
    
  }
  
  # Adding elapsed time
  if (include.timing) {
    
    gof <- c(
      gof,
      ifelse(
        inherits(model, "ergmito_boot"),
        model$timer_boot["total"],
        model$timer["total"]
        )
      )
    
    gof.names   <- c(gof.names, "Time (seconds)")
    gof.decimal <- c(gof.decimal, TRUE)
      
  }
  
  return(
    texreg::createTexreg(
      coef.names  = coefficient.names,
      coef        = coefficients,
      se          = standard.errors,
      pvalues     = significance,
      gof.names   = gof.names,
      gof         = gof,
      gof.decimal = gof.decimal
      )
  )
  
}


setMethod(
  "extract", signature = className("ergmito", "ergmito"),
  definition = extract.ergmito
)

setMethod(
  "extract", signature = className("ergmito_boot", "ergmito"),
  definition = extract.ergmito
)
