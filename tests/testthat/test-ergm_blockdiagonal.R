context("Blockdiagonal model using the ergm framework")

test_that("Blockdiagonalizing works as expected", {
  
  data(fivenets)
  net1 <- blockdiagonalize(fivenets)
  net1 <- splitnetwork(net1, "block")
  net2 <- blockdiagonalize(as_adjmat(fivenets)) # Regardless of whether it has attrs or not
  net2 <- splitnetwork(net2, "block")
  
  expect_equal(as_adjmat(net1), as_adjmat(fivenets))
  expect_equal(as_adjmat(fivenets), as_adjmat(net2))
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