# Simulation of a probit

set.seed(1231)
n <- 2000
k <- 4

b <- rnorm(k)
X <- matrix(runif(n*k), ncol=k)
Y <- as.integer(X%*%b + rnorm(n) > 0)

ans0 <- glm(Y~0+X, family = binomial("probit"))

mle <- function(p) {
  X1 <- X[Y==1,]
  X0 <- X[Y==0,]
  sum(log(pnorm(- X0 %*% p[-1]/p[1]))) +
    sum(log(1 - pnorm(-X1 %*% p[-1]/p[1])))
}

ans1 <- optim(c(1, rnorm(k)), mle, control=list(fnscale=-1))
ans0$coefficients
ans1$par[-1]/ans1$par[1]

# Models -----------------------------------------------------------------------

set.seed(1)
b <- c(-2, 4, 1)/10
X <- cbind(
  rbinom(n, 2, .5),
  rbeta(n, 8, 2),
  runif(n, 0, 10)
  )

# Simulating the data
# pr <- plogis(X %*% b)
# pr <- pbeta(exp(X %*% b))

# Generating the data
# Y <- rbinom(n, 10, pr)
# X <- scale(X)
xp <- exp(X %*% b)
Y <- rbeta(n, xp, 1.5)
# Y <- rlogis(n, exp(X %*% b))
hist(Y, breaks=50)
library(stats4)
beta_mle <- function(Y, X) {
  
  call <- match.call()
  ll <- function(p)
    sum(log(dbeta(Y, exp(X %*% p[-1]), exp(p[1]))))
  
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
    vcov      = structure(-solve(sol$hessian), dimnames = list(cnames, cnames)),
    min       = sol$value,
    details   = sol,
    minuslogl = function(x) -ll(x),
    nobs      = length(Y),
    method    = "BFGS"
    )
  
}


ans <- beta_mle(Y, X)

summary(ans)
