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

n     <- 3
npars <- 2

sets <- powerset(n)

library(parallel)
library(sna)
library(ergm)



stats <- parallel::mclapply(sets, function(s) {
  
  ans <- matrix(0, ncol=n, nrow=n)
  ans[s+1] <- 1
  summary(ans ~ mutual + edges)
})

S <- do.call(rbind, stats)
Smeans <- colMeans(S)["edges"]
S <- S - matrix(Smeans, nrow = nrow(S), ncol=ncol(S), byrow = TRUE)

set.seed(1155)
net    <- as.matrix(netdiffuseR::rgraph_er(n, p = .75))
stats0 <- summary(net ~ mutual + edges) - Smeans
ll <- function(params) {
  
  P <- matrix(params, ncol=length(params), nrow=nrow(S), byrow = TRUE)
  sum(stats0*params) - log(sum(exp(P * S)))
  
}

gr <- function(params) {
  
  P <- matrix(params, ncol=length(params), nrow=nrow(S), byrow = TRUE)
  sum(stats0) - 1/sum(exp(P * S))*sum(exp(P * S)*S)
  
}

sol  <- optim(par = rnorm(npars), fn = ll, gr = gr, control=list(fnscale=-1, reltol=1e-25, maxit=2e3));sol
sol2 <- ergm(net ~ mutual + edges, control =control.ergm(MCMC.samplesize = 2e3,MCMC.interval = 2e3));summary(sol2) 
ll(sol2$coef)
sol$par
ll(sol$par)
sol$par;sqrt(diag(-solve(optimHess(sol$par, ll))))

ergm.exact(sol2$coef, sol2$formula)
ergm.exact(sol$par, sol2$formula)
ergm.allstats(sol2$formula)
unique(S)
