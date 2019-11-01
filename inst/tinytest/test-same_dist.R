
data("fivenets")
expect_true(same_dist(fivenets[[1]], fivenets[[2]]))
expect_false(same_dist(fivenets[[1]], fivenets[[2]], "female"))

expect_error(same_dist(fivenets[[1]], as_adjmat(fivenets[[1]])))
expect_error(same_dist(fivenets[[1]], fivenets[[2]], "females"))
  
data("fivenets")

nets <- rbernoulli(c(4, 5))
expect_true(same_dist(nets[[1]], nets[[1]]))
expect_false(same_dist(nets[[1]], nets[[2]]))

expect_error(same_dist(nets[[1]], fivenets[[1]]))
# expect_error(same_dist(nets[[1]], nets[[2]], "females"))

data("fivenets")

# Single variable
net1 <- fivenets[[1]]
network::set.vertex.attribute(net1, "test", c(.1, .8, .2, .4))

net2 <- network::permute.vertexIDs(
  network::network.copy(net1),
  v = c(2,4,1,3)
  )

ans <- same_dist(net1, net2, "test")

expect_equal(
  network::get.vertex.attribute(net1, "test")[attr(ans, "map01")],
  network::get.vertex.attribute(net2, "test")
)

expect_equal(
  network::get.vertex.attribute(net2, "test")[attr(ans, "map10")],
  network::get.vertex.attribute(net1, "test")
)

# More than 1
net1 <- fivenets[[1]]
network::set.vertex.attribute(net1, "test1", c(.1, .8, .2, .4))
network::set.vertex.attribute(net1, "test2", c(.5, 0, .3, .1))

net2 <- network::permute.vertexIDs(
  network::network.copy(net1),
  v = c(2,4,1,3)
)

ans <- same_dist(net1, net2, c("test1", "test2"))

expect_equal(
  network::get.vertex.attribute(net1, "test1")[attr(ans, "map01")],
  network::get.vertex.attribute(net2, "test1")
)

expect_equal(
  network::get.vertex.attribute(net2, "test1")[attr(ans, "map10")],
  network::get.vertex.attribute(net1, "test1")
)

expect_equal(
  network::get.vertex.attribute(net1, "test2")[attr(ans, "map01")],
  network::get.vertex.attribute(net2, "test2")
)

expect_equal(
  network::get.vertex.attribute(net2, "test2")[attr(ans, "map10")],
  network::get.vertex.attribute(net1, "test2")
)

