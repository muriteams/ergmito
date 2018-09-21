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
sets <- powerset(4)
nets <- lapply(sets, function(s) {
  
  ans <- matrix(0, ncol=4, nrow=4)
  ans[s+1] <- 1
  ans
})

library(parallel)
library(sna)
library(ergm)

stats <- mclapply(nets, function(n) {
  summary(n ~ edges + balance + mutual)
}, mc.cores = 4)

S <- do.call(rbind, stats)

set.seed(1123)
net    <- sna::rgraph(4, tprob = .7)
stats0 <- summary(net ~ edges + balance + mutual)
ll <- function(params) {
  
  P <- matrix(params, ncol=3, nrow=nrow(S), byrow = TRUE)
  sum(stats0*params) - sum(P * S)
  
}

sol <- ABCoptim::abc_optim(c(.1, .1, .1), ll, fnscale=-1)
