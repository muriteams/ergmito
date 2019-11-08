pset <- powerset(4)

# List method
ans0 <- matrix_to_network(pset)
ans1 <- lapply(pset, as.network)

expect_equal(as_adjmat(ans0), as_adjmat(ans1))

expect_equal(
  as_adjmat(matrix_to_network(pset[[106]])),
  as_adjmat(as.network(pset[[106]]))
  )
