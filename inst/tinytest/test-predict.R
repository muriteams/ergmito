data(fivenets)

# bernoulli graph
fit <- ergmito(fivenets ~ edges) 

# all ties have the same likelihood
ans0 <- mean(nedges(fivenets)/(nvertex(fivenets)*(nvertex(fivenets) - 1)))  
ans1 <- predict(fit)
ans1 <- lapply(ans1, `diag<-`, NA)
ans1 <- sapply(ans1, as.vector)

ans1 <- ans0 - ans1[!is.na(ans1)]
expect_equivalent(ans1, rep(0, length(ans1)), tol = 1e-5)

# bernoulli graph
fit <- ergmito(fivenets ~ edges + nodematch("female"))

net1 <- fivenets[[1]]
net1_list_net <- fivenets[1]

ans0 <- predict(fit)[1]
ans1 <- predict(fit, newdata = net1)
ans2 <- predict(fit, newdata = net1_list_net)

expect_equal(ans0, ans1)
expect_equal(ans0, ans2)

fit <- ergmito(fivenets ~ edges + ttriad)

ans0 <- predict(fit)[1]
ans1 <- predict(fit, newdata = net1)
ans2 <- predict(fit, newdata = net1_list_net)
ans3 <- predict(fit, newdata = as_adjmat(net1))
ans4 <- predict(fit, newdata = as_adjmat(net1_list_net))

expect_equal(ans0, ans1)
expect_equal(ans0, ans2)
expect_equal(ans0, ans3)
expect_equal(ans0, ans4)


