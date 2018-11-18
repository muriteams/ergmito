context("Likelihood function")

test_that("Higher than ergm", {
  
  # Pre-fitted ergm
  load("test-data-for-tests.rda")
  
  ans1 <-lergm(net ~ mutual + edges)
  
  ans0 <- ergm::ergm.exact(coef(ans0), net ~ mutual + edges)
  ans1 <- ergm::ergm.exact(coef(ans1), net ~ mutual + edges)
  
  expect_lt(ans0, ans1)
  
})

test_that("Order doesn't matter", {
  
  set.seed(8871)
  net1 <- rbernoulli(4)
  net2 <- rbernoulli(5)
  
  set.seed(1717171);ans0 <- lergm(list(net1, net2) ~ edges + mutual)
  set.seed(1717171);ans1 <- lergm(list(net2, net1) ~ edges + mutual)
  
  expect_equal(coef(ans0), coef(ans1))
  expect_equal(vcov(ans0), vcov(ans1))
  
})

test_that("Multiple nets", {
  
  set.seed(121)
  net1 <- rbernoulli(4)
  set.seed(1000); ans0 <- lergm(net1 ~ edges + mutual, zeroobs = FALSE)
  set.seed(1000); ans1 <- lergm(list(net1, net1) ~ edges + mutual)
  
  expect_equal(coef(ans0), coef(ans1))
  expect_equal(vcov(ans0), vcov(ans1)*2) # Vcov is half of it!
  
})

