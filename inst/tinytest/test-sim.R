
# Single network test
set.seed(4532)
nets <- rbernoulli(8)

expect_error(new_rergmito(nets ~ edges, theta = -1, force = FALSE), "force")

set.seed(4532)
nets <- rbernoulli(3)

expect_error(new_rergmito(
  list(nets, network::network(nets, directed = FALSE)) ~ edges, "mix"
))

Sys.unsetenv("ERGMITO_TEST")
set.seed(1);ans0 <- new_rergmito(nets ~ edges + mutual, theta=c(-.5, .5))
Sys.setenv(ERGMITO_TEST = 1)
set.seed(1);ans1 <- new_rergmito(nets ~ edges + mutual, theta=c(-.5, .5), force = TRUE)

expect_equal(ans0$prob, ans1$prob)
expect_equal(sum(ans0$probabilities), 1)

# Multiple network
set.seed(45532)
nets <- rbernoulli(c(3,3,3))

Sys.unsetenv("ERGMITO_TEST")
set.seed(1);ans0 <- new_rergmito(nets ~ edges + mutual, theta = c(-.5, .5))
Sys.setenv(ERGMITO_TEST = 1)
set.seed(1);ans1 <- new_rergmito(nets ~ edges + mutual, theta = c(-.5, .5), force = TRUE)

expect_equal(ans0[[1]]$probabilities, ans1[[1]]$probabilities)
expect_equal(sum(ans0[[1]]$probabilities), 1)

# Multiple sizes
Sys.unsetenv("ERGMITO_TEST")
set.seed(45532)
nets <- rbernoulli(c(3,4))
set.seed(1);ans0 <- new_rergmito(nets ~ edges + mutual, theta=c(-.5,.5))
Sys.setenv(ERGMITO_TEST = 1)

expect_equal(sum(ans0[[1]]$probabilities), 1)
expect_equal(sum(ans0[[2]]$probabilities), 1)

# Sampler with the data
data("fivenets")
mod <- ergmito(fivenets ~ edges + nodematch("female"))
network::delete.vertices(fivenets[[1]], 1)

Sys.unsetenv("ERGMITO_TEST")
ans0 <- new_rergmito(fivenets[[1]] ~ edges + nodematch("female"), theta = coef(mod))
Sys.setenv(ERGMITO_TEST = 1)
ans1 <- new_rergmito(fivenets[[1]] ~ edges + nodematch("female"), theta = coef(mod), force = TRUE)
Sys.unsetenv("ERGMITO_TEST")

expect_equal(ans0$probabilities, ans1$probabilities)
expect_equal(sum(ans0$probabilities), 1)
  

set.seed(4532)
nets <- rbernoulli(3)
ans  <- new_rergmito(nets ~ edges + mutual, theta = c(-.5, .5))

expect_output(print(ans), "ERGMito simulator")
expect_output(print(str(ans)), "language new_rergmito")
expect_equal(ans$get_networks(), powerset(3))
expect_equal(ans$get_networks(), ans[])

Sys.unsetenv("ERGMITO_TEST")


# Undirected networks ----------------------------------------------------------
set.seed(554)
nets <- rbernoulli(4, .1)
nets <- network::network(nets, directed = FALSE)
ans <- ergmito(nets ~ edges + triangles)
coef_ans <- coef(ans)
ans  <- suppressWarnings(
  new_rergmito(
    nets ~ edges + triangles, force = TRUE,
    theta = coef_ans
    ))

expect_length(ans$networks, 2^(4*3/2))
expect_equal(
  as_adjmat(ans$get_networks()),
  powerset(4, directed = FALSE)
  )
expect_equal(sum(ans$probabilities), 1)


