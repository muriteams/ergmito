library(magrittr)
library(stats4)

set.seed(65454)

source("models/beta_mle.R")

n_sims  <- 500
n_teams <- 50
n       <- replicate(n_sims, sample(c(3, 4, 5), n_teams, TRUE), simplify = FALSE)
dens    <- replicate(n_sims, runif(n_teams), simplify = FALSE)
prec    <- replicate(n_sims, runif(n_teams), simplify = FALSE)

theta   <- replicate(n_sims, rnorm(3), simplify=FALSE)
X       <- lapply(n, function(n0) {
  cbind(
    X1 = rbinom(n_teams, n0, .5)/n0,
    X2          = rnorm(n_teams)
    )
  })


# Simple example
z <- sim_experiment(
  n     = n[[2]],
  dens  = dens[[2]],
  prec  = prec[[2]],
  X     = X[[2]],
  theta = c(-2, 1, .5)
  )

# Extracting the data
d <- cbind(
  h = sapply(z, "[[", "prec_hat"),
  y = sapply(z, "[[", "response"),
  X = X[[1]]
)

ans <- beta_mle(d[,"y"], d[,-2])
summary(ans)

# Correlation
plot(
  sapply(z, "[[", "prec"),
  sapply(z, "[[", "prec_hat")
)


ans <- parallel::mcmapply(
  sim_experiment,
  n = n, dens=dens, prec=prec, X=X, theta=theta, mc.cores = 8,
  SIMPLIFY = FALSE
  )

# Estimating models ------------------------------------------------------------

mles0 <- parallel::mcmapply(function(dat, x) {
  
  # Extracting the data
  d <- cbind(
    y = sapply(dat, "[[", "response"),
    hamming_distance = sapply(dat, "[[", "prec_hat"),
    X = x
  )
  
  beta_mle(d[,"y"], d[,-1])
  
}, dat = ans, x=X, mc.cores=8, SIMPLIFY=FALSE)

mles1 <- parallel::mcmapply(function(dat, x) {
  
  # Extracting the data
  d <- cbind(
    y = sapply(dat, "[[", "response"),
    group_size = sapply(dat, "[[", "n"),
    X = x
  )
  
  beta_mle(d[,"y"], d[,-1])
  
}, dat = ans, x=X, mc.cores=8, SIMPLIFY=FALSE)


mles2 <- parallel::mcmapply(function(dat, x) {
  
  # Extracting the data
  d <- cbind(
    y = sapply(dat, "[[", "response"),
    hamming_distance = sapply(dat, "[[", "prec_hat"),
    group_size = sapply(dat, "[[", "n"),
    X = x
  )
  
  beta_mle(d[,"y"], d[,-1])
  
}, dat = ans, x=X, mc.cores=8, SIMPLIFY=FALSE)



# Computing pvalues ------------------------------------------------------------

pvals0 <- lapply(mles0, function(model) calc_pval(coef(model), vcov(model)))
pvals0 <- do.call(rbind, pvals0)
boxplot(1 - pvals0, main="Power")

pvals1 <- lapply(mles1, function(model) calc_pval(coef(model), vcov(model)))
pvals1 <- do.call(rbind, pvals1)
boxplot(1 - pvals1, main="Power")

pvals2 <- lapply(mles2, function(model) calc_pval(coef(model), vcov(model)))
pvals2 <- do.call(rbind, pvals2)
boxplot(1 - pvals2, main="Power")


# Correlation plots

# nvar  <- unlist(n)
# score <- lapply(ans, function(a) sapply(a, "[[", "response")) %>%
#   do.call(c, .)
# prec_hat  <- lapply(ans, function(a) sapply(a, "[[", "prec_hat")) %>%
#   do.call(c, .)
# prec  <- lapply(ans, function(a) sapply(a, "[[", "prec")) %>%
#   do.call(c, .)

cor(prec, nvar)
smoothScatter(prec, prec_hat)
