
# Directed powersets -----------------------------------------------------------
pset2 <- powerset(2, directed = TRUE)
pset3 <- powerset(3, directed = TRUE)
pset4 <- powerset(4, directed = TRUE)
pset5 <- powerset(5, directed = TRUE)

expect_length(pset2, 4)
expect_length(pset3, 2^6)
expect_length(pset4, 2^(4*3))
expect_length(pset5, 2^(5*4))

net <- network::network(rbernoulli(3), directed = TRUE)

# Computing using powersets
pset <- powerset(3, directed = TRUE, force = TRUE)
ans1 <- lapply(pset, function(i) {
  summary(i ~ edges + mutual)
})

ans1 <- do.call(rbind, ans1)
ans1 <- table(paste0(ans1[,1], ans1[,2]))
ans1 <- structure(as.vector(ans1), names = names(ans1))
ans1 <- ans1[order(names(ans1))]

# Baseline
ans0 <- ergm::ergm.allstats(net ~ edges + mutual, zeroobs = FALSE)
ans0 <- with(
  ans0,
  structure(weights, names = paste0(statmat[,1], statmat[,2]))
)
ans0 <- ans0[order(names(ans0))]

expect_equal(ans0, ans1)  

# Undirected powersets ---------------------------------------------------------
pset2 <- powerset(2, directed = FALSE)
pset3 <- powerset(3, directed = FALSE)
pset4 <- powerset(4, directed = FALSE)
pset5 <- powerset(5, directed = FALSE)

expect_length(pset2, 2)
expect_length(pset3, 2^(6/2))
expect_length(pset4, 2^(4*3/2))
expect_length(pset5, 2^(5*4/2)) 

net <- network::network(rbernoulli(4), directed = FALSE)

# Computing using powersets
pset <- powerset(4, directed = FALSE, force = TRUE)
pset <- lapply(pset, network::network, directed = FALSE)
ans1 <- lapply(pset, function(i) {
  summary(i ~ edges + cycle(3))
})

ans1 <- do.call(rbind, ans1)
ans1 <- table(paste0(ans1[,1], ans1[,2]))
ans1 <- structure(as.vector(ans1), names = names(ans1))
ans1 <- ans1[order(names(ans1))]

# Baseline
ans0 <- ergm::ergm.allstats(net ~ edges + cycle(3), zeroobs = FALSE)
ans0 <- with(
  ans0,
  structure(weights, names = paste0(statmat[,1], statmat[,2]))
  )
ans0 <- ans0[order(names(ans0))]

expect_equal(ans0, ans1)  

