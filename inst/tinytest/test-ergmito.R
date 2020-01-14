
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

options(ergmito_warning = TRUE)
expect_warning(ans <- ergmito(fivenets ~ edges + mutual, gattr=~ y), "\"y\"")
options(ergmito_warning = FALSE)
expect_output(print(ans$formulae), "elements by using")

set.seed(8871)
net1 <- rbernoulli(4)
net2 <- rbernoulli(5)

set.seed(1717171);ans0 <- ergmito(list(net1, net2) ~ edges + mutual)
set.seed(1717171);ans1 <- ergmito(list(net2, net1) ~ edges + mutual)

expect_equal(coef(ans0), coef(ans1), tolerance = 1e-4)
expect_equal(vcov(ans0), vcov(ans1), tolerance = 1e-4)


set.seed(121)
net1 <- rbernoulli(4)
set.seed(1000); ans0 <- ergmito(net1 ~ edges + mutual, zeroobs = TRUE)
set.seed(1000); ans1 <- ergmito(list(net1, net1) ~ edges + mutual)

expect_equal(coef(ans0), coef(ans1), tolerance = 1e-4)
expect_equal(vcov(ans0), vcov(ans1)*2, tolerance = 1e-4) # Vcov is half of it!

expect_output(print(ans0), "ERGMito")
expect_output(print(summary(ans0)), "z value")

net1 <- network::network(net1, directed = FALSE)
options(ergmito_warning = TRUE)
expect_warning(ergmito(net1 ~ edges), "undirected graphs")
options(ergmito_warning = FALSE)

# Errors
expect_error(ergmito(list(net1, rbernoulli(4)) ~ edges), "same type")
expect_error(ergmito(net1 ~ gwdsp(1)), "not supported")

