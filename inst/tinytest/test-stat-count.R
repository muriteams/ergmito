# Generating data
x <- powerset(5)
ans0 <- count_stats(x[1:100], c("mutual", "edges"))

# Calculating using summary_formula
fm <- x[[i]] ~ mutual + edges
ans1 <- lapply(1:100, function(i) {
  environment(fm) <- environment()
  ergm::summary_formula(fm)
})
ans1 <- do.call(rbind, ans1)

expect_equivalent(ans0, ans1)

# Matches the formula syntax
ans2 <- count_stats(x[1:100] ~ mutual + edges)
expect_equivalent(ans0, ans2)
  

data(fivenets)
ans0 <- geodesita(fivenets)
ans1 <- lapply(fivenets, sna::geodist, inf.replace = NA_integer_)
ans1 <- lapply(ans1, "[[", "gdist")

# On the fivenets
expect_equal(ans0, ans1)

# On powerset of 4
pset <- powerset(4)
ans0 <- geodesita(pset)
ans1 <- lapply(pset, sna::geodist, inf.replace = NA_integer_)
ans1 <- lapply(ans1, "[[", "gdist")
expect_equal(ans0, ans1)

# Small matrices
ans0 <- geodesita(pset[[2]])
ans1 <- sna::geodist(pset[[2]], inf.replace = NA_integer_)$gdist
ans2 <- geodesita(network::as.network(pset[[2]]))
expect_equal(ans0[[1]], ans1)
expect_equal(ans0, ans2)


# Bug in ergm?????? ----------------------------------------------------------
set.seed(1)
pset3 <- powerset(3) # rbernoulli(rep(6, 100))
ans0 <- count_stats(pset3, c("ttriad"))[,1]
ans1 <- unname(sapply(pset3, function(p) summary(p ~ ttriad)))

expect_equivalent(ans0, ans1)
#-----------------------------------------------------------------------------

estats <- c(
  `mutual`          = "mutual",
  `edges`           = "edges",
  `ctriad`          = "ctriad",
  `ttriad`          = "ttriad",
  `triangle`        = "triangle",
  `balance`         = "balance",
  `triadcensus(15)` = "t300",
  `triadcensus(2)`  = "t102",
  `idegree1.5`      = "idegree1.5",
  `odegree1.5`      = "odegree1.5"
)

for (s in seq_along(estats)) {
  ans0 <- count_stats(pset3, c(estats[s]))[,1]
  ans1 <- unname(sapply(pset3, function(p) {
    fm <- as.formula(paste0("p ~ ", names(estats)[s]))
    environment(fm) <- environment()
    summary(fm)
    }))
  
  expect_equivalent(ans0, ans1)
}

pset5 <- powerset(5)
pset5 <- pset5[order(sapply(pset5, sum))]
pset5 <- c(head(pset5, 25), tail(pset5, 25))
for (star in c("istar", "ostar"))
for (i in 1:4) {
  ans0 <- count_stats(pset5, paste0(star, i))[,1]
  fm <- as.formula(paste0("p ~ ",star, "(", i, ")"))
  ans1 <- unname(sapply(pset5, function(p) {
    environment(fm) <- environment()
    summary(fm)
    }))
  
  expect_equivalent(ans0, ans1)
}


# Attribute based ------------------------------------------------------------
set.seed(44)
age <- lapply(nvertex(pset5), rpois, lambda=4)

ans0 <- count_stats(pset5, "nodeocov", age)
ans1 <- unname(sapply(seq_along(pset5), function(i) {
  summary(network::network(pset5[[i]], list(age = age[[i]]), "age") ~ nodeocov("age"))
  }))

expect_equivalent(as.vector(ans0), as.vector(ans1))

ans0 <- count_stats(pset5, "nodeicov", age)
ans1 <- unname(sapply(seq_along(pset5), function(i) {
  summary(network::network(pset5[[i]], list(age = age[[i]]), "age") ~ nodeicov("age"))
}))

expect_equivalent(as.vector(ans0), as.vector(ans1))

ans0 <- count_stats(pset5, "absdiff", age)
ans1 <- unname(sapply(seq_along(pset5), function(i) {
  summary(network::network(pset5[[i]], list(age = age[[i]]), "age") ~ absdiff("age"))
}))

expect_equivalent(as.vector(ans0), as.vector(ans1))

for (star in c("istar", "ostar"))
  for (i in 1:4) {
    ans0 <- as.vector(count_stats(pset5, paste0(star, i), age))
    ans1 <- integer(length(pset5))
    for (g in seq_along(pset5)) {
      net <- network::network(pset5[[g]], list(age = age[[g]]), "age")
      fm <- as.formula(paste0("net ~ ", star, "(",i, ", attr=\"age\")"))

      ans1[g] <- summary(fm)
    }

    expect_equivalent(ans0, ans1)
  }




# network(pset3[[1]], list(age = age[[1]])) %>%
#           gplot(label = age[[1]])
  
data("fivenets")
expect_error(
  ergmito(fivenets ~ edges + edgecitos),
  "not found"
)

expect_error(analyze_formula())

# Error when analyzing undirected networks
expect_error(
  count_stats(
    network::network(rbernoulli(4), directed = FALSE) ~ edges
    ), "undirected"
  )

# Not available
expect_error(count_stats(rbernoulli(4) ~ missing_term), "not available")
