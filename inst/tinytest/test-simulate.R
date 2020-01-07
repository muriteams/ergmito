data(fivenets)
fit <- ergmito(fivenets ~ edges + mutual)

expect_error(simulate(fit, 200, which_networks = 100))
expect_error(simulate(fit, 200, which_networks = -1))

expect_silent(ans <- simulate(fit, 200))
expect_silent(ans <- simulate(fit, 200, cl = NULL, ncores = 2L))

expect_error(ans <- simulate(fit, 200, cl = NULL, ncores = 2L, sizes = c(3,4)))


# Refit with simulated data should yield the same
ans0 <- ergmito(fivenets ~ edges + nodematch("female"))
set.seed(8819)
ans1 <- simulate(ans0, 20, which_networks = 1:5)
ans1 <- ergmito(ans1 ~ edges + nodematch("female"))
expect_equal(coef(ans0), coef(ans1), tol = .25)

  
# Simulations with attribute based network -------------------------------------
ans0  <- new_rergmito(
  fivenets[[3]] ~ edges + nodematch("female"),
  theta = coef(ergmito(fivenets ~ edges + nodematch("female")))
  )
set.seed(12)
s     <- ans0$sample(500, s = 4)
attrs <- lapply(s, network::get.vertex.attribute, attrname = "female")
attrs <- do.call(rbind, attrs)
attrs <- t(attrs) - network::get.vertex.attribute(fivenets[[3]], attrname = "female")
expect_true(all(attrs == 0))

# Recovering parameters
expect_equal(
  coef(ergmito(s ~ edges + nodematch("female"))),
  ans0$theta, tol = .2
)





