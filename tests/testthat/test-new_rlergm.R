context("Sampler")

test_that("Using the count_stats vs ergm::summary_formula is equal", {
  
  set.seed(45532)
  nets <- rbernoulli(c(4,3,4))
  
  Sys.setenv(LERGM_TEST = "")
  set.seed(1);ans0 <- new_rlergm(nets ~ edges + mutual, sizes = 2:3, mc.cores = 2L)
  Sys.setenv(LERGM_TEST = 1)
  set.seed(1);ans1 <- new_rlergm(nets ~ edges + mutual, sizes = 2:3, mc.cores = 2L)
  
  expect_equal(ans0$prob, ans1$prob)
  
  
  
})
