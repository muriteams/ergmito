context("Utils")

test_that("as.adjmat is equivalent to as.matrix", {
  
  data(sampson, package="ergm")
  data(faux.desert.high, package="ergm")
  
  expect_equivalent(as.adjmat(samplike), as.matrix(samplike))
  expect_equivalent(as.adjmat(faux.desert.high), as.matrix(faux.desert.high))
  
})