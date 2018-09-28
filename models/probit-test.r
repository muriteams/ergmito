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
b <- c(-2, 4, 1)
X <- cbind(
  rbinom(n, 2, .5),
  rbeta(n, 8, 2),
  runif(n, 0, 10)
  )

# Simulating the data
pr <- plogis(X %*% b)

# Generating the data
Y <- rbinom(n, 10, pr)

# Log-likelihood
mle <- function(p) 
  sum(log(dbinom(Y, 10, plogis(X %*% p))))
  

ans1 <- optim(rnorm(3), mle, control=list(fnscale = -1, reltol=1e-20), hessian = TRUE)

# Results
cat(sprintf(
  "%8.4f (%0.2f)",
  ans1$par,
  sqrt(-diag(solve(ans1$hessian)))
), sep="\n")

