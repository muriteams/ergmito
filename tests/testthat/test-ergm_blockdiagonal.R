context("Blockdiagonal model using the ergm framework")

test_that("Blockdiagonalizing works as expected", {
  
  data(fivenets)
  suppressMessages(require(network))
  
  set.seed(771)
  fivenets <- lapply(fivenets, function(net) {
    
    net %v% "rand1" <- runif(network.size(net))
    net
    
  })
  
  net1 <- blockdiagonalize(fivenets)
  net1 <- splitnetwork(net1, "block")
  net2 <- blockdiagonalize(as_adjmat(fivenets)) # Regardless of whether it has attrs or not
  net2 <- splitnetwork(net2, "block")
  
  # Structure is preserved
  expect_equal(as_adjmat(net1), as_adjmat(fivenets))
  expect_equal(as_adjmat(fivenets), as_adjmat(net2))
  
  # Attributes are preserved
  expect_equal(
    lapply(net1, get.vertex.attribute, "rand1"),
    lapply(fivenets, get.vertex.attribute, "rand1")
  )
})

test_that("Blockdiagonal model give somewhat the same answer", {
  data(fivenets)
  net0 <- blockdiagonalize(fivenets)
  
  ans0 <- ergmito(fivenets ~ edges + nodematch("female"))
  require(ergm)
  suppressMessages(
    ans1 <- ergm(net0 ~ edges + nodematch("female"), constraints = ~ blockdiag("block"))
  )
  
  expect_equal(coef(ans0), coef(ans1), tol=.5)
  
})