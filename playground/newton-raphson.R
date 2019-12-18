nr <- function(x, fn, d, H, ..., delta = 1.0, tol = 1e-5, maxiter = 20) {
  
  time0 <- Sys.time()
  d0 <- d(x, ...)
  i  <- 0L
  x0 <- x
  H0 <- solve(H(x0, ...))
  while (i < maxiter) {
    
    i  <- i + 1L

    # Updating the value
    H1 <- tryCatch(solve(H(x0, ...)), error = function(e) e)
    if (inherits(H1, "error")) {
      H1 <- H0
    }
    
    step <- tryCatch(H1 %*% d0 * delta, error =function(e) e)
    if (any(step > 100)) {
      step <- step / max(abs(step)) * 2
    }
    
    x1   <- x0 - step
    d1   <- d(x1, ...)
    
    # Computing the values
    if (any(!is.finite(d1))) {
      
      delta <- delta * .75
      next

    } else if (norm(d1) < tol)
      break
    
    x0 <- x1
    d0 <- d1
    H0 <- H1
  }
  
  time1 <- Sys.time()
  
  list(
    par         = as.vector(x1),
    value       = fn(x1, ...),
    counts      = c(`function` = 0, gradient = i),
    convergence = ifelse(i == maxiter, 1, 0),
    message     = NULL,
    hessian     = H(x1, ...),
    time        = difftime(time1, time0, units = "s"),
    delta       = delta,
    d1          = d1
  )
  
}


library(ergmito)

data("fivenets")
set.seed(1331)
nets <- rbernoulli(rep(5, 100), .8)
model <- nets ~ edges + istar(2) + ostar(2)
f <- ergmito_formulae(model)


ans0 <- nr(
  rep(1, ncol(f$target.stats)),
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

exact_hessian(
  params = ans1$par,
  x      = f$target.stats,
  stats.weights = f$stats.weights,
  stats.statmat = f$stats.statmat,
)
ans1$hessian

exact_gradient(
  params = ans1$par,
  x      = f$target.stats,
  stats.weights = f$stats.weights,
  stats.statmat = f$stats.statmat,
)
