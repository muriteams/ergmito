
# Pre-fitted ergm
load(system.file("test-data-for-tests.rda", package="ergmito"))

ans1 <-ergmito(net ~ mutual + edges)

ans0 <- ergm::ergm.exact(coef(ans0), net ~ mutual + edges)
ans1 <- ergm::ergm.exact(coef(ans1), net ~ mutual + edges)

expect_lt(ans0, ans1)
  

data(fivenets)
set.seed(121)
A <- rbinom(5, 1, .3)
for (i in seq_along(fivenets))
  network::set.network.attribute(fivenets[[i]], "y" ,A[i])

set.seed(8871)
net1 <- rbernoulli(4)
net2 <- rbernoulli(5)

set.seed(1717171);ans0 <- ergmito(list(net1, net2) ~ edges + mutual)
set.seed(1717171);ans1 <- ergmito(list(net2, net1) ~ edges + mutual)

expect_equal(coef(ans0), coef(ans1), tol = 1e-4, scale = 1)
expect_equal(vcov(ans0), vcov(ans1), tol = 1e-4, scale = 1)


set.seed(121)
net1 <- rbernoulli(4)
set.seed(1000); ans0 <- ergmito(net1 ~ edges + mutual)
set.seed(1000); ans1 <- ergmito(list(net1, net1) ~ edges + mutual)

expect_equal(coef(ans0), coef(ans1), tol = 1e-4, scale = 1)
expect_equal(vcov(ans0), vcov(ans1)*2, tol = 1e-4, scale = 1) # Vcov is half of it!

expect_output(print(ans0), "ERGMito")
expect_output(print(summary(ans0)), "z value")

net1 <- network::network(net1, directed = FALSE)
options(ergmito_warning = TRUE)
expect_warning(ergmito(net1 ~ edges), "undirected graphs")
options(ergmito_warning = FALSE)

# Errors
expect_error(ergmito(list(net1, rbernoulli(4)) ~ edges), "same type")
expect_error(ergmito(net1 ~ gwdsp(1)), "not supported")

# Offset parameters ------------------------------------------------------------
data(fivenets)

ans <- ergmito(
  fivenets ~ edges + nodematch("female") + ttriad,
  model_update = ~ . - offset(I(ttriple/2)) - ttriple
  )

# expect_output(print(ans), "offset")
expect_output(print(summary(ans)), "offset")
expect_length(coef(ans), 2)

ans <- suppressWarnings(ergmito_boot(ans, R = 10))
expect_length(vcov(ans)[1,], 2)

# expect_output(print(ans), "offset")
expect_output(print(summary(ans)), "offset")
plot(ans)

# Infinite offset
ans0 <- ergmito(fivenets ~ edges, model_update = ~ .+offset(I(ifelse(edges < 3, -Inf, 0))))

f <- ergmito_formulae(fivenets ~ edges)

# Updating everything
ids <- c(2, 3, 4)
f$target_stats  <- f$target_stats[ids, , drop = FALSE]
f$target_offset <- f$target_offset[ids]
f$stats_statmat <- f$stats_statmat[ids]
f$stats_weights <- f$stats_weights[ids]
f$stats_offset  <- f$stats_offset[ids]

for (i in 1:3) {
  
  off <- which(f$stats_statmat[[i]] < 3)
  f$stats_statmat[[i]] <- f$stats_statmat[[i]][-off, , drop = FALSE]
  f$stats_weights[[i]] <- f$stats_weights[[i]][-off]
  f$stats_offset[[i]] <- f$stats_offset[[i]][-off]
  
}

ans1 <- with(f, ergmito(
  fivenets[ids] ~ edges, 
  target_stats  = target_stats,
  target_offset = target_offset,
  stats_weights = stats_weights,
  stats_statmat = stats_statmat,
  stats_offset  = stats_offset
  ))
                                                            
expect_equal(coef(ans1), coef(ans0))
expect_equal(logLik(ans1), logLik(ans0))
expect_equal(vcov(ans1), vcov(ans0))
