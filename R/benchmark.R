#' Utility to benchmark expression in R
#' 
#' This is just an internal utility included in the package which is not designed
#' to be accurate. If you need accurate benchmarks, you should take a look at the
#' \CRANpkg{microbenchmark} and \CRANpkg{bench} R packages.
#' 
#' @param ... List of expressions to benchmark.
#' @param times Integer scalar. Number of replicates.
#' @param rand_ord Logical. When `TRUE`, the expressions are executed in a random
#' order.
#' @details 
#' The print method prints a summary including quantiles, relative elapsed times,
#' and the order (fastest to slowest).
#' 
#' @return A data frame of class `ergmito_benchmark` with `times` rows and the
#' following columns:
#' 
#' - `id` Integer. Id of the expression.
#' - `expr` Factor. Expression executed.
#' - `user.self`, `sys.self`, `elapsed`, `sys.child` time in seconds (see [proc.time()]).
#' 
#' Includes the following attributes: `ncalls`, `call`, `label`, and `expr`.
#' 
#' @examples 
#' bm <- benchmarkito(
#'   exp(1:100000),
#'   sqrt(1:100000),
#'   times = 20
#' )
#' 
#' plot(bm)
#' print(bm)
#' 
#' @export
benchmarkito <- function(..., times = 100, rand_ord = TRUE) {
 args  <- as.list(match.call())[-1L]
 n_args <- names(args)
 if (is.null(n_args))
   names(args) <- sapply(args, deparse)
 else {
   are_unnamed <- (n_args == "") | is.na(n_args)
   names(args)[are_unnamed] <- sapply(args[are_unnamed], deparse)
 }

 args <- args[setdiff(names(args), c("times", "rand_ord"))]
 ncalls <- length(args)
 ord <- rep(1:ncalls, times = times)
 if (rand_ord)
   ord <- sample(ord)

 ans  <- lapply(ord, function(i) {
   t0 <- proc.time()
   eval(args[[i]])
   proc.time() - t0
   })
 ans  <- cbind(
   data.frame(id = ord, expr = names(args)[ord]),
   as.data.frame(do.call(rbind, ans))
 )
 
 structure(
   ans,
   ncalls = ncalls,
   call   = match.call(),
   label  = names(args),
   expr   = sapply(args, deparse),
   class  = c("ergmito_benchmark", "data.frame")
 )
 
}

#' @export
print.ergmito_benchmark <- function(x, ...) {
  
  averages <- lapply(1:attr(x, "ncalls"), function(i) {
    x. <- x[x$id == i, setdiff(colnames(x), c("id", "expr"))]
    
    # Centiles
    q. <- rbind(stats::quantile(x.$elapsed, probs = c(.05, .5, .95)))
    cbind(
      rbind(colMeans(as.matrix(x.))),
      q.
      )
  })
  averages <- as.data.frame(do.call(rbind, averages))
  
  # Relative with respect to the fastest
  fastest           <- which.min(averages[, "elapsed"])[1]
  averages$relative <- averages$elapsed / averages$elapsed[fastest]
  averages$ord      <- rank(averages$elapsed, ties.method = "random")
  
  cat("call:\n  ")
  print(attr(x, "call"))
  cat("\nTimes (in seconds, see proc.time()):\n\n")
  print(data.frame(label = attr(x, "label"), averages, check.names = FALSE))
  invisible(x)
  
}

#' @export
plot.ergmito_benchmark <- function(x, y = NULL, ...) {
  
  ans <- do.call(cbind, lapply(1:attr(x, "ncalls"), function(i) {
    x[x$id == i, "elapsed"]
  }))
  
  colnames(ans) <- attr(x, "label")
  graphics::boxplot(ans, ...)
  
}
