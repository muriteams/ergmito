data(fivenets)
fit <- ergmito(fivenets ~ edges + mutual)

expect_error(simulate(fit, 200, which_networks = 100))
expect_error(simulate(fit, 200, which_networks = -1))

expect_silent(ans <- simulate(fit, 200))
  
