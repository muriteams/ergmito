compute_mfrow <- function(k) {
  
  if (k == 1) c(1,1)
  else if (k == 2) c(1, 2)
  else if (k < 5) c(2, 2)
  else if (k < 7) c(3, 2)
  else if (k < 10) c(3, 3)
  else if (k < 13) c(4, 3)
  else c(4, 4)
  
}


#' Function to visualize the optimization surface
#' 
#' General diagnostics function. This function allows to visualize the surface
#' to be maximize at around a particular point.
#' 
#' It calculates the surface coordinates for each pair of parameters included
#' in the ERGMito.
#' 
#' @param x An object of class [ergmito].
#' @param y,... Ignored.
#' @param domain A list.
#' @param plot. Logical. When `TRUE` (default), the function will call [graphics::image]
#' and plot all possible combination of parameters.
#' @param par_args Further arguments to be passed to [graphics::par]
#' @param image_args Further arguments to be passed to [graphics::image]
#' @param extension Numeric. Range value of the function.
#' @return A list of length `choose(length(object$coef), 2)` (all possible
#' combinations of pairs of parameters), each with the following elements:
#' - `z` A matrix
#' - `z` A vector
#' - `y` A vector
#' - `xlab` A string. Name of the ERGM parameter in the x-axis.
#' - `ylab` A string. Name of the ERGM parameter in the y-axis.
#' 
#' The list is returned invisible. 
#' @export
#' @examples 
#' 
#' set.seed(12)
#' x1 <- sna::rgraph(4)
#' x2 <- sna::rgraph(5)
#' 
#' ans <- ergmito(list(x1, x2) ~ edges + mutual + balance)
#' 
#' plot(ans)
#' 
#' @seealso The [ergmito] function.
#' @importFrom graphics image par
#' @importFrom viridisLite viridis
#' @importFrom utils combn
plot.ergmito <- function(
  x,
  y = NULL,
  domain     = NULL,
  plot.      = TRUE,
  par_args   = list(),
  image_args = list(),
  extension  = 10,
  ...
  ) {
  
  if (!inherits(x, "ergmito"))
    stop("This function only accepts objects of class `ergmito`.", call. = FALSE)
  
  f <- x$formulae$loglik
  k <- x$formulae$npars
  
  # Over what domain should we calculate this?
  if (!length(domain)) {
   
    if (extension<0)
      stop("`extension` should be a positive value.", call. = FALSE)
    
    # Local region
    domain <- lapply(x$coef, function(z) {
      seq(z - extension, z + extension, length.out = 100)
    })
    
  }
  
  # Range of values to go through
  averages <- sapply(domain, mean)
  npoints  <- sapply(domain, length)
  
  # Pairs of
  parcombs <- utils::combn(k, 2, simplify = FALSE)
  
  # Setting up plotting device
  if (plot.) {
    
    if (length(par_args$mfrow) == 0)
      par_args$mfrow <- compute_mfrow(length(parcombs))
    
    op <- do.call(graphics::par, par_args)
    on.exit(graphics::par(op))
    
  }
  
  Z <- vector("list", length(parcombs))
  for (p in seq_along(parcombs)) {
    
    # Obtaining the indices
    i0 <- parcombs[[p]][1]
    j0 <- parcombs[[p]][2]
    
    # Creating the matrix space for this
    Z[[p]] <- list(
      z = matrix(NA, nrow = npoints[i0], ncol = npoints[j0])
    )
    par0 <- coef(x)
    
    for (i in 1:npoints[i0]) 
      for (j in 1:npoints[j0]) {
        
        par0[i0] <- domain[[i0]][i]
        par0[j0] <- domain[[j0]][j]
        
        # Calculating the loglikelihood for that set of observations
        Z[[p]]$z[i, j] <- f(par0)
        
      }
    
    Z[[p]]$x <- domain[[i0]]
    Z[[p]]$y <- domain[[j0]]
    Z[[p]]$xlab <- names(par0)[i0]
    Z[[p]]$ylab <- names(par0)[j0]
    
    # Should I plot?
    if (plot.) {
      do.call(graphics::image, c(Z[[p]], list(col = viridisLite::viridis(30))))
     
      graphics::points(x = coef(x)[i0], y = coef(x)[j0], col="red", cex=1.5,
             pch = 20)
       
    }
    
  }
  
    
  invisible(Z)
  
}


