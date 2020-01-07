
data(sampson, package="ergm")
data(faux.desert.high, package="ergm")

expect_equivalent(as_adjmat(samplike), as.matrix(samplike))
expect_equivalent(as_adjmat(faux.desert.high), as.matrix(faux.desert.high))

# Checking which are directed or undirected
set.seed(713)
nets <- rbernoulli(c(4,4,4))
nets0 <- lapply(nets, network::network, directed = FALSE)
nets1 <- lapply(nets, network::network, directed = TRUE)
expect_equal(is_directed(nets0), rep(FALSE, 3))
expect_equal(is_directed(nets1), rep(TRUE, 3))
expect_true(!is_directed(nets0[[1]]))
expect_true(is_directed(nets0[[1]][1:4,1:4]))
expect_error(is_directed(nets0[[1]][1:4,1:4], check_type = TRUE))

