nr <- function(x, fn, d, H, ..., tol = 1e-10, maxiter = 100) {
  
  time0 <- Sys.time()
  
  d0 <- d(x, ...)
  d1 <- d0
  i  <- 0L
  while (i < maxiter) {
    
    i  <- i + 1L
    
    # Updating the value
    x <- x - solve(H(x, ...)) %*% d0 / 2 # solve(H(x, ...)) %*% d0
      
    # Computing the values
    d1 <- d(x, ...)

    if (norm(d1) < tol)
      break
    
    # message("Current value of d(x) = ", norm(d1))
    
    d0 <- d1
  }
  
  time1 <- Sys.time()
  
  list(
    par         = as.vector(x),
    value       = fn(x, ...),
    counts      = c(`function` = 0, gradient = i),
    convergence = ifelse(i == maxiter, 1, 0),
    message     = NULL,
    hessian     = H(x, ...),
    time        = difftime(time1, time0, units = "s")
  )
  
}


library(ergmito)

data("fivenets")
set.seed(1331)
nets <- rbernoulli(rep(5, 40), .2)
model <- nets ~ edges # + mutual
f <- ergmito_formulae(model)

ans0 <- nr(
  rep(0, 1),
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
ans0$value < ans1$value

ans0$time
ans1$time["optim"]
