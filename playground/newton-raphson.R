nr <- function(x, fn, d, H, ..., tol = 1e-20, maxiter = 100) {
  
  d0 <- d(x, ...)
  d1 <- d0
  i  <- 0L
  while (i < maxiter) {
    
    i  <- i + 1L
    
    # Updating the value
    x <- x - solve(H(x, ...)) %*% d0
      
    # Computing the values
    d1 <- d(x, ...)

    if (norm(d1) < tol)
      break
    
    message("Current value of d(x) = ", d1)
    
    d0 <- d1
  }
  
  list(
    par         = as.vector(x),
    value       = fn(x, ...),
    counts      = c(`function` = 0, gradient = i),
    convergence = ifelse(i == maxiter, 1, 0),
    message     = NULL,
    hessian     = H(x, ...)
  )
  
}


library(ergmito)

data("fivenets")
set.seed(13331)
nets <- rbernoulli(rep(5, 20), .2)
f <- ergmito_formulae(nets ~ edges + mutual)

ans0 <- nr(
  c(1, 1),
  d  = f$grad,
  H  = f$hess, 
  fn = f$loglik,
  stats.weights = f$stats.weights,
  stats.statmat = f$stats.statmat,
  target.stats  = f$target.stats
  )
ans0
ans1 <- ergmito(nets ~ edges + ttriad)
ans1$optim.out
ans0$value - ans1$optim.out$value
