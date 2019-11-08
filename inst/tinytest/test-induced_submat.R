# One on one
set.seed(5554)
x   <- rbernoulli(100, .1)
net <- network::as.network(x)

ids <- c(1:10, 91:100)

ans0 <- x[ids, ids]
ans1 <- induced_submat(x, ids)
ans2 <- induced_submat(net, ids)

expect_equal(ans0, ans1[[1]])
expect_equivalent(as.matrix(net)[ids, ids], ans2[[1]])
  
# Multiple indices
set.seed(5554)
x    <- rbernoulli(100, .1)
net  <- network::as.network(x)
ids1 <- c(1:10, 91:100)
ids2 <- c(50:60, 1:10)

ans0 <- list(x[ids1, ids1], x[ids2, ids2])
ans1 <- induced_submat(x, list(ids1, ids2))
ans2 <- induced_submat(net, list(ids1, ids2))

# Multiple indices
expect_equal(ans0, ans1)
expect_equivalent(
  list(
    as.matrix(net)[ids1, ids1],
    as.matrix(net)[ids2, ids2]
    ), ans2)
  
set.seed(5554)
x    <- rbernoulli(100, .1)
net  <- network::as.network(x)
ids1 <- c(1:10, 91:100)
ids2 <- c(50:60, 1:10)

ans0 <- list(x[ids1, ids1], x[ids2, ids2])
ans1 <- induced_submat(list(x, x), list(ids1, ids2))
ans2 <- induced_submat(list(net, net), list(ids1, ids2))

expect_equal(ans0, ans1)
expect_equivalent(
  list(
    as.matrix(net)[ids1, ids1],
    as.matrix(net)[ids2, ids2]
  ), ans2)

# Multiple networks, single vector
ans0 <- list(x[ids1, ids1], x[ids1, ids1])
ans1 <- induced_submat(list(x,x), ids1)
ans2 <- induced_submat(list(net,net), ids1)
expect_equal(ans0, ans1)
expect_equal(ans0, ans2)

# Errors
expect_error(induced_submat(list(), list()), "length 0")
expect_error(induced_submat(list(x,x,x), list(ids1, ids2)), "same as")
expect_error(induced_submat(list(matrix(nrow=1, ncol=3)), ids1), "square")
expect_error(induced_submat(x, c(1,1)), "repeated")
expect_error(induced_submat(x, c(-1,1)), "boundary")
expect_error(induced_submat(x, c(1,-1)), "boundary")

