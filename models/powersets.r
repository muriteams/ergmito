rm(list=ls())

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


set.seed(1155)
net    <- as.matrix(netdiffuseR::rgraph_er(n, p = .75))
stats0 <- summary(net ~ mutual + edges)

S <- do.call(rbind, stats)
S <- S - matrix(stats0, nrow=nrow(S), ncol=ncol(S), byrow = TRUE)
stats0 <- stats0 - stats0

# Optimization
ll <- function(params) {
  
  - log(sum(exp(S %*% params)))
  
}

gr <- function(params) {
  
  # Sum(exp())
  exp_sum <- exp(S %*% params)

  - 1/log(sum(exp_sum))*(t(S) %*% exp_sum)
  

}

sol  <- optim(
  par    = rnorm(npars),
  fn     = ll,
  gr     = gr,
  method = "BFGS",
  control= list(
    fnscale = -1,
    reltol  = 1e-25,
    maxit   = 2e3
    )
  )
sol

# ERGM
sol2 <- ergm(
  net ~ mutual + edges,
  control = control.ergm(MCMC.samplesize = 2e3,MCMC.interval = 2e3))

# Comparing side by side

ll(sol2$coef)
ll(sol$par)
sol$par;sqrt(diag(-solve(optimHess(sol$par, ll))))
summary(sol2) 

ergm.exact(sol2$coef, sol2$formula)
ergm.exact(sol$par, sol2$formula)

# Comparing statistics ---------------------------------------------------------
X0 <- ergm.allstats(sol2$formula)
X0 <- do.call(cbind, X0)[,c(2:3,1)]
X0 <- X0[order(X0[,1], X0[,2]),]

# My version
X1 <- S
X1 <- table(X1[,1], X1[,2])
X1 <- cbind(
  expand.grid(
    as.integer(rownames(X1)),
    as.integer(colnames(X1))
    ),
  as.vector(X1)
)

X1 <- X1[X1[,3] > 0,]
X1 <- X1[order(X1[,1], X1[,2]),]
X0-X1
