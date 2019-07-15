context("Network wrappers")

test_that("Matrix to network works", {
  
  pset <- powerset(3)
  
  # List method
  ans0 <- matrix_to_network(pset)
  ans1 <- lapply(pset, as.network)
  
  expect_equal(as_adjmat(ans0), as_adjmat(ans1))
  
})