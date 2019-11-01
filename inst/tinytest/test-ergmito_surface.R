context("Surface function works just fine")

test_that("Simple call with 2 parameters", {
  
  set.seed(1)
  net <- rbernoulli(rep(5, 5))
  ans <- ergmito(net ~ edges + mutual)
  
  expect_silent(plot(ans))
  
  
})