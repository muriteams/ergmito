
# Single network test
set.seed(4532)
nets <- rbernoulli(5)

expect_error(new_rergmito(nets ~ edges, sizes = 8, force = FALSE), "force")

expect_error(new_rergmito(
  list(nets, network::network(nets, directed = FALSE)) ~ edges, "mix"
))

Sys.unsetenv("ERGMITO_TEST")
set.seed(1);ans0 <- new_rergmito(nets ~ edges + mutual, sizes = 3, mc.cores = 1L)
Sys.setenv(ERGMITO_TEST = 1)
set.seed(1);ans1 <- new_rergmito(nets ~ edges + mutual, sizes = 3, mc.cores = 1L, force = TRUE)

expect_equal(ans0$prob, ans1$prob)
expect_equal(sum(ans0$prob$`3`), 1)

# Multiple network
set.seed(45532)
nets <- rbernoulli(c(4,3,4))

Sys.unsetenv("ERGMITO_TEST")
set.seed(1);ans0 <- new_rergmito(nets ~ edges + mutual, sizes = 3, mc.cores = 1L)
Sys.setenv(ERGMITO_TEST = 1)
set.seed(1);ans1 <- new_rergmito(nets ~ edges + mutual, sizes = 3, mc.cores = 1L, force = TRUE)

expect_equal(ans0$prob, ans1$prob)
expect_equal(sum(ans0$prob$`3`), 1)

# Sampler with the data
data("fivenets")
mod <- ergmito(fivenets ~ edges + nodematch("female"))
network::delete.vertices(fivenets[[1]], 1)

Sys.unsetenv("ERGMITO_TEST")
ans0 <- new_rergmito(fivenets[[1]] ~ edges + nodematch("female"), theta = coef(mod), mc.cores = 1L)
Sys.setenv(ERGMITO_TEST = 1)
ans1 <- new_rergmito(fivenets[[1]] ~ edges + nodematch("female"), theta = coef(mod), mc.cores = 1L, force = TRUE)
Sys.unsetenv("ERGMITO_TEST")

expect_equal(ans0$prob, ans1$prob)
expect_equal(sum(ans0$prob$`3`), 1)
  

set.seed(4532)
nets <- rbernoulli(5)
ans  <- new_rergmito(nets ~ edges + mutual, sizes = 3, mc.cores = 1L)

expect_output(print(ans), "ERGMito simulator")
expect_output(print(str(ans)), "language new_rergmito")
expect_error(ans$get_networks(s = 5), "sampling")
expect_error(ans[,5], "sampling")
expect_equal(ans$get_networks(s = 3), powerset(3))
expect_equal(ans$get_networks(s = 3), ans[,3]$`3`)

Sys.unsetenv("ERGMITO_TEST")

# Parallel
ans  <- suppressWarnings({
  new_rergmito(nets ~ edges + mutual, sizes = 3, force = TRUE, ncores = 2)
  })
expect_length(ans$networks$`3`, 2^(2*3))
expect_equal(
  as_adjmat(ans$get_networks(s = 3)),
  powerset(3, directed = TRUE)
)
expect_equal(sum(ans$prob$`3`), 1)


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

expect_length(ans$networks$`4`, 2^(4*3/2))
expect_equal(
  as_adjmat(ans$get_networks(s = 4)),
  powerset(4, directed = FALSE)
  )
expect_equal(sum(ans$prob$`4`), 1)

# Parallel
ans  <- suppressWarnings(
  new_rergmito(
    nets ~ edges + esp(2), force = TRUE, ncores = 2,
    theta = coef_ans
    )
  )

expect_length(ans$networks$`4`, 2^(4*3/2))
expect_equal(
  as_adjmat(ans$get_networks(s = 4)),
  powerset(4, directed = FALSE)
)
expect_equal(sum(ans$prob$`4`), 1)

# Checking errors
