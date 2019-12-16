nr <- function(x, fn, d, H, ..., tol = 1e-5, maxiter = 100) {
  
  time0 <- Sys.time()
  
  i  <- 0L
  while (i < maxiter) {
    
    i  <- i + 1L
    
    # Updating the value
    d0 <- d(x, ...)
    x1 <- x - solve(H(x, ...)) %*% d0 # solve(H(x, ...)) %*% d0
    d1 <- d(x1, ...)
      
    # Computing the values
    if (norm(d1 - d0) < tol)
      break
    
    message("Current value of x = ", norm(d1 - d0))
    x <- x1
    
  }
  
  time1 <- Sys.time()
  
  list(
    par         = as.vector(x1),
    value       = fn(x1, ...),
    counts      = c(`function` = 0, gradient = i),
    convergence = ifelse(i == maxiter, 1, 0),
    message     = NULL,
    hessian     = H(x1, ...),
    time        = difftime(time1, time0, units = "s")
  )
  
}


library(ergmito)

data("fivenets")
set.seed(1331)
nets <- rbernoulli(rep(5, 40), .2)
model <- fivenets ~ edges + nodematch("female")
f <- ergmito_formulae(model)

ans0 <- nr(
  rep(0, 2),
  d  = f$grad,
  H  = f$hess, 
  fn = f$loglik,
  stats.weights = f$stats.weights,
  stats.statmat = f$stats.statmat,
  target.stats  = f$target.stats
  )
ans0
ans1 <- ergmito(model)
ans1 <- c(ans1$optim.out, list(time = ans1$time["optim"], value = logLik(ans1)))
ans0$value - ans1$value

ans0$time;ans1$time

ans0$counts;ans1$counts

ans0$par;ans1$par