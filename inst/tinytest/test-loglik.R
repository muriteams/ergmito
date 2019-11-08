# Generating the data
params <- c(.5, .8)
nets   <- powerset(3)

# Computing statistics per network
model <- n ~ edges + mutual
stats <- lapply(nets, function(n) {
  
  environment(model) <- environment()
  summary(model)
  
})
stats <- do.call(rbind, stats)

# Summary statistics
allstats <- ergm::ergm.allstats(nets[[1]] ~ edges + mutual, zeroobs = FALSE)

# Function to compute the loglike of all networks
ergmloglik <- function(x) {
  
  apply(x, 1, function(z) {
    
    ergmito:::exact_loglik2(params, z, allstats)
    
  })
  
}

ans0 <- ergmloglik(stats)
ans1 <- exact_loglik(stats, params, allstats$weights, allstats$statmat)

expect_equal(sum(exp(ans0)), 1, tol = 1e-5)
expect_equal(sum(exp(ans1)), 1, tol = 1e-5)

expect_equivalent(ans0[], ans1[])
  

# # Are these equal
# all(parloglik(stats) ==
#       exact_ll(stats, params, allstats$weights, allstats$statmat))