library(microbenchmark)
library(ergmito)

Rcpp::sourceCpp("playground/fast-loglik.cpp")

params <- c(.5, .8)

nets  <- powerset(4)

# Computing statistics per network
model <- n ~ edges + mutual
stats <- parallel::mclapply(nets, function(n) {
  
  environment(model) <- environment()
  summary(model)
  
}, mc.cores = 4L)
stats <- do.call(rbind, stats)

# Summary statistics
allstats <- ergm::ergm.allstats(nets[[1]] ~ edges + mutual, zeroobs = FALSE)

# Function to compute the loglike of all networks
parloglik <- function(x) {
  
  apply(x, 1, function(z) {
    
    ergmito:::exact_loglik(params, z, allstats)
    
  })
  
}

# Are these equal
all(parloglik(stats) ==
  exact_ll(stats, params, allstats$weights, allstats$statmat))

microbenchmark(
  apply = parloglik(stats),
  rcpp1 = exact_ll(stats, params, allstats$weights, allstats$statmat, 1),
  rcpp4 = exact_ll(stats, params, allstats$weights, allstats$statmat, 4)
)


microbenchmark(
  apply = parloglik(stats[1,,drop=FALSE]),
  rcpp1 = exact_ll(stats[1,,drop=FALSE], params, allstats$weights, allstats$statmat, 1),
  rcpp4 = exact_ll(stats[1,,drop=FALSE], params, allstats$weights, allstats$statmat, 4)
)


