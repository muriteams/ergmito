context("Sampler")

test_that("Using the count_stats vs ergm::summary_formula is equal", {
  
  # Single network test
  set.seed(4532)
  nets <- rbernoulli(5)
  Sys.setenv(ERGMITO_TEST = "")
  set.seed(1);ans0 <- new_rergmito(nets ~ edges + mutual, sizes = 3, mc.cores = 1L)
  Sys.setenv(ERGMITO_TEST = 1)
  set.seed(1);ans1 <- new_rergmito(nets ~ edges + mutual, sizes = 3, mc.cores = 1L, force = TRUE)
  
  expect_equal(ans0$prob, ans1$prob)
  expect_equal(sum(ans0$prob$`3`), 1)
  
  # Multiple network
  set.seed(45532)
  nets <- rbernoulli(c(4,3,4))
  
  Sys.setenv(ERGMITO_TEST = "")
  set.seed(1);ans0 <- new_rergmito(nets ~ edges + mutual, sizes = 3, mc.cores = 1L)
  Sys.setenv(ERGMITO_TEST = 1)
  set.seed(1);ans1 <- new_rergmito(nets ~ edges + mutual, sizes = 3, mc.cores = 1L, force = TRUE)
  
  expect_equal(ans0$prob, ans1$prob)
  expect_equal(sum(ans0$prob$`3`), 1)
  
  # Sampler with the data
  data("fivenets")
  mod <- ergmito(fivenets ~ edges + nodematch("female"))
  network::delete.vertices(fivenets[[1]], 1)

  Sys.setenv(ERGMITO_TEST = "")  
  ans0 <- new_rergmito(fivenets[[1]] ~ edges + nodematch("female"), theta = coef(mod), mc.cores = 1L)
  Sys.setenv(ERGMITO_TEST = 1)
  ans1 <- new_rergmito(fivenets[[1]] ~ edges + nodematch("female"), theta = coef(mod), mc.cores = 1L, force = TRUE)
  Sys.setenv(ERGMITO_TEST = "")
  
  expect_equal(ans0$prob, ans1$prob)
  expect_equal(sum(ans0$prob$`3`), 1)
  
})

test_that("Methods work as expected", {
  
  set.seed(4532)
  nets <- rbernoulli(5)
  ans  <- new_rergmito(nets ~ edges + mutual, sizes = 3, mc.cores = 1L)
  
  expect_output(print(ans), "ERGMito simulator")
  expect_output(print(str(ans)), "language new_rergmito")
  expect_error(ans$get_networks(s = 5), "sampling")
  expect_error(ans[,5], "sampling")
  expect_equal(ans$get_networks(s = 3), powerset(3))
  expect_equal(ans$get_networks(s = 3), ans[,3]$`3`)
  
  
})