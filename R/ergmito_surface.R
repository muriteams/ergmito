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
#' @param breaks Integer scalar. Number of splits per dimension.
#' @param params_labs Named vector. Alternative labels for the parameters. It 
#' should be in the form of `c("orignial name" = "new name")`.
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
#' x <- rbernoulli(c(4, 4, 5))
#' 
#' ans <- ergmito(x ~ edges + balance)
#' 
#' plot(ans)
#' 
#' @seealso The [ergmito] function.
#' @importFrom graphics image par plot title lines
#' @importFrom utils combn
#' @importFrom stats setNames
plot.ergmito <- function(
  x,
  y          = NULL,
  domain     = NULL,
  plot.      = TRUE,
  par_args   = list(),
  image_args = list(),
  breaks     = 50L,
  extension  = 4L,
  params_labs = stats::setNames(names(coef(x)), names(coef(x))),
  ...
  ) {
  
  # Checking that the alternative labels are defined for all parameters
  # in the model.
  if (!all(names(params_labs) %in% names(coef(x)))) {
    stop("When specifying `param_labs`, all parameters should be specified.",
         call. = FALSE)
  }
  params_labs <- params_labs[names(coef(x))]
  
  if (!inherits(x, "ergmito"))
    stop("This function only accepts objects of class `ergmito`.", call. = FALSE)
  
  f <- function(p) {
    with(
      x$formulae,
      loglik(
        params        = p,
        target.stats  = target.stats,
        stats.weights = stats.weights,
        stats.statmat = stats.statmat
        )
    )}
  
  k <- x$formulae$npars
  
  # Over what domain should we calculate this?
  if (!length(domain)) {
   
    if (extension<0)
      stop("`extension` should be a positive value.", call. = FALSE)
    
    # Local region
    domain <- lapply(x$coef, function(z) {
      seq(z - extension, z + extension, length.out = breaks)
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
    # if (length(par_args$mar) == 0)
    #   par_args$mar <- graphics::par("mar") * .5
    
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
        
        # Calculating the log-likelihood for that set of observations
        Z[[p]]$z[i, j] <- f(par0)
        
      }
    
    Z[[p]]$x <- domain[[i0]]
    Z[[p]]$y <- domain[[j0]]
    Z[[p]]$xlab <- params_labs[i0]
    Z[[p]]$ylab <- params_labs[j0]
    
    # Should I plot?
    if (plot.) {
      
      # dput(viridisLite::viridis(100))
      viridis30 <- c("#440154FF", "#450558FF", "#46085CFF", "#470D60FF", "#471063FF", 
                     "#481467FF", "#481769FF", "#481B6DFF", "#481E70FF", "#482173FF", 
                     "#482576FF", "#482878FF", "#472C7AFF", "#472F7CFF", "#46327EFF", 
                     "#453581FF", "#453882FF", "#443B84FF", "#433E85FF", "#424186FF", 
                     "#404587FF", "#3F4788FF", "#3E4A89FF", "#3D4D8AFF", "#3C508BFF", 
                     "#3B528BFF", "#39558CFF", "#38598CFF", "#375B8DFF", "#355E8DFF", 
                     "#34608DFF", "#33638DFF", "#32658EFF", "#31688EFF", "#2F6B8EFF", 
                     "#2E6D8EFF", "#2D708EFF", "#2C718EFF", "#2B748EFF", "#2A768EFF", 
                     "#29798EFF", "#287C8EFF", "#277E8EFF", "#26818EFF", "#26828EFF", 
                     "#25858EFF", "#24878EFF", "#238A8DFF", "#228D8DFF", "#218F8DFF", 
                     "#20928CFF", "#20938CFF", "#1F968BFF", "#1F998AFF", "#1E9B8AFF", 
                     "#1F9E89FF", "#1FA088FF", "#1FA287FF", "#20A486FF", "#22A785FF", 
                     "#24AA83FF", "#25AC82FF", "#28AE80FF", "#2BB07FFF", "#2EB37CFF", 
                     "#31B67BFF", "#35B779FF", "#39BA76FF", "#3DBC74FF", "#41BE71FF", 
                     "#47C06FFF", "#4CC26CFF", "#51C56AFF", "#56C667FF", "#5BC863FF", 
                     "#61CA60FF", "#67CC5CFF", "#6DCD59FF", "#73D056FF", "#78D152FF", 
                     "#7FD34EFF", "#85D54AFF", "#8CD646FF", "#92D741FF", "#99D83DFF", 
                     "#A0DA39FF", "#A7DB35FF", "#ADDC30FF", "#B4DE2CFF", "#BBDE28FF", 
                     "#C2DF23FF", "#C9E020FF", "#D0E11CFF", "#D7E219FF", "#DDE318FF", 
                     "#E4E419FF", "#EBE51AFF", "#F1E51DFF", "#F7E620FF", "#FDE725FF"
      )
      
      do.call(graphics::image, c(Z[[p]], list(col = viridis30)))
     
      graphics::points(x = coef(x)[i0], y = coef(x)[j0], col="red", cex=1.5,
             pch = 20)
       
    }
    
  }
  
    
  invisible(Z)
  
}


