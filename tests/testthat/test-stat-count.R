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
  
})