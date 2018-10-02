# The MLE routine depends on the stats4 R package
library(stats4, quietly = TRUE)

# Function to compute a beta mle
beta_mle <- function(Y, X) {
  
  call <- match.call()
  ll <- function(p) {
    ans <- sum(log(dbeta(Y, exp(X %*% p[-1]), exp(p[1]))))
    
    if (is.infinite(ans))
      return(.Machine$double.xmax*sign(ans)*1e-3)
    
    ans
  }
  
  sol <- optim(
    rnorm(ncol(X) + 1),
    ll,
    control=list(fnscale = -1, reltol=1e-20),
    hessian = TRUE
  )
  
  cnames <- colnames(X)
  if (!length(cnames))
    cnames <- paste0("X", 1:ncol(X))
  
  cnames <- c("beta", cnames)
  
  
  new(
    "mle",
    call      = call,
    coef      = structure(sol$par, names = cnames),
    fullcoef  = structure(sol$par, names = cnames),
    vcov      = structure(-MASS::ginv(sol$hessian), dimnames = list(cnames, cnames)),
    min       = sol$value,
    details   = sol,
    minuslogl = function(x) -ll(x),
    nobs      = length(Y),
    method    = "BFGS"
  )
  
}

calc_pval <- function(b, information) {
  
  p <- pnorm(b, sd = sqrt(diag(information)))
  
  if (any(is.na(p)))
    return(rep(NA, length(p)))
  
  ifelse(p > .5, 1-p, p)*2
  
  
}

# Simulates and undirected graph with density d
sim_graph <- function(n, d) {
  
  G    <- matrix(0, ncol=n, nrow=n)
  G[which(runif(n*n) < d)] <- 1L
  diag(G) <- 0L
  
  G
  
}

# With probability p reverts cell content
permute_graph <- function(G, p) {
  
  n <- nrow(G)
  swap <- which(runif(n*n) < (1-p))
  
  G[swap] <- 1 - G[swap]
  
  diag(G) <- 0L
  
  G
  
}

# Computes the hamming distance between the observed and predicted network.
hamming <- function(G0, G, i) {
  
  n <- ncol(G0)-1
  1 - sum(abs(G0[-i,-i] - G[-i, -i]))/(n*(n-1))
  
}

#' Function to simulate a team of size `n`
sim_team <- function(n, dens, prec) {
  
  # Step 1: True network
  G0 <- sim_graph(n, dens)
  
  # Step 2: Permuted versions
  if (length(prec) == 1)
    prec <- rep(prec, n)
  
  G <- lapply(prec, permute_graph, G = G0)
  
  # Step 3: Computing hamming distance
  list(
    G0       = G0,
    G        = G,
    prec     = mean(prec),
    prec_hat = mean(
      unlist(Map(hamming, G0 = replicate(n,G0, simplify = FALSE), G = G, i=1:n))
      ),
    n        = n
  )
}

#' Simulates an experiment with `length(n)` teams
sim_experiment <- function(n, dens, prec, X, theta) {
  
  dat <- Map(sim_team, n = n, dens=dens, prec=prec)
  
  # Extracting the mean precision per team
  mean_prec <- sapply(dat, "[[", "prec")
  
  response <- cbind(X, mean_prec) %*% theta
  response <- rbeta(length(n), exp(response), 1.5)
  
  dat <- Map(function(a, b) c(a, list(response=b)), a = dat, b=response)
  
  dat
  
}


