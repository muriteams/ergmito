powerset <- function(n) {
  
  set   <- 1:(n*(n-1))
  
  sets <- NULL
  for (i in set) {
    for (s in sets)
      sets <- c(sets, list(c(s, i)))
    sets <- c(sets, list(i))
    # print(sets)
  }
  
  c(sets, list(NULL))
  
}

Rcpp::sourceCpp("models/powersets.cpp")

n <- 4
sets <- powerset(n)

library(parallel)
library(sna)
library(ergm)

stats <- parallel::mclapply(sets, function(s) {
  
  ans <- matrix(0, ncol=n, nrow=n)
  ans[s+1] <- 1
  summary(ans ~ edges + meandeg + mutual)
})

S <- do.call(rbind, stats)

set.seed(1123)
net    <- sna::rgraph(n, tprob = .4)
stats0 <- summary(net ~ edges + meandeg + mutual)
ll <- function(params) {
  
  P <- matrix(params, ncol=3, nrow=nrow(S), byrow = TRUE)
  sum(stats0*params) - log(sum(exp(P * S)))
  
}

sol  <- optim(rnorm(3), ll, control=list(fnscale=-1, reltol=1e-25, maxit=2e3));sol
sol2 <- ergm(net ~ edges + meandeg + mutual);summary(sol2)
logLik(sol2)
sol
ll(sol2$coef)
