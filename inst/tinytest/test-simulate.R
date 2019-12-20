data(fivenets)
fit <- ergmito(fivenets ~ edges + mutual)

expect_error(simulate(fit, 200, which_networks = 100))
expect_error(simulate(fit, 200, which_networks = -1))

expect_silent(ans <- simulate(fit, 200))
expect_silent(ans <- simulate(fit, 200, cl = NULL, ncores = 2L))

expect_error(ans <- simulate(fit, 200, cl = NULL, ncores = 2L, sizes = c(3,4)))
