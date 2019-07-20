context("Counting stats")

test_that("Matches ERGM", {
  
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
  
})

test_that("Geodesic distances", {
  
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
  
})

test_that("Sufficient statistics", {
  # Bug in ergm?????? ----------------------------------------------------------
  set.seed(1)
  pset3 <- powerset(3) # rbernoulli(rep(6, 100))
  ans0 <- count_stats(pset3, c("ttriad"))[,1]
  ans1 <- unname(sapply(pset3, function(p) summary(p ~ ttriad)))
  
  expect_equivalent(ans0, ans1)
  #-----------------------------------------------------------------------------
  
  ans0 <- count_stats(pset3, c("mutual"))[,1]
  ans1 <- unname(sapply(pset3, function(p) summary(p ~ mutual)))
  
  expect_equivalent(ans0, ans1)
  
  ans0 <- count_stats(pset3, c("edges"))[,1]
  ans1 <- unname(sapply(pset3, function(p) summary(p ~ edges)))
  
  expect_equivalent(ans0, ans1)
  
  ans0 <- count_stats(pset3, c("ctriad"))[,1]
  ans1 <- unname(sapply(pset3, function(p) summary(p ~ ctriad)))
  
  expect_equivalent(ans0, ans1)
  
  ans0 <- count_stats(pset3, c("triangle"))[,1]
  ans1 <- unname(sapply(pset3, function(p) summary(p ~ triangle)))
  
  expect_equivalent(ans0, ans1)
  
  ans0 <- count_stats(pset3, c("balance"))[,1]
  ans1 <- unname(sapply(pset3, function(p) summary(p ~ balance)))
  
  expect_equivalent(ans0, ans1)
  
  ans0 <- count_stats(pset3, c("t300"))[,1]
  ans1 <- unname(sapply(pset3, function(p) summary(p ~ triadcensus(15))))
  
  expect_equivalent(ans0, ans1)
  
  ans0 <- count_stats(pset3, c("t102"))[,1]
  ans1 <- unname(sapply(pset3, function(p) summary(p ~ triadcensus(2))))
  
  expect_equivalent(ans0, ans1)
  
  # Nodecovar
  set.seed(44)
  age <- lapply(nvertex(pset3), rpois, lambda=4)
  
  ans0 <- count_stats(pset3, "nodeocov", age)
  ans1 <- unname(sapply(seq_along(pset3), function(i) {
    summary(network::network(pset3[[i]], list(age = age[[i]]), "age") ~ nodeocov("age"))
    }))
  
  expect_equivalent(ans0, ans1)
  
  # Nodecovar
  set.seed(44)
  age <- lapply(nvertex(pset3), rpois, lambda=4)
  
  ans0 <- count_stats(pset3, "nodeicov", age)
  ans1 <- unname(sapply(seq_along(pset3), function(i) {
    summary(network::network(pset3[[i]], list(age = age[[i]]), "age") ~ nodeicov("age"))
  }))
  
  expect_equivalent(ans0, ans1)
  
  
  # network(pset3[[1]], list(age = age[[1]])) %>%
  #           gplot(label = age[[1]])
  
})

test_that("checking existance using ergm", {
  data("fivenets")
  expect_error(
    analyze_formula(fivenets ~ edges + edgecitos, check_w_ergm = TRUE),
    "not found"
  )
  
})