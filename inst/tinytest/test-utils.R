
data(sampson, package="ergm")
data(faux.desert.high, package="ergm")

expect_equivalent(as_adjmat(samplike), as.matrix(samplike))
expect_equivalent(as_adjmat(faux.desert.high), as.matrix(faux.desert.high))
  