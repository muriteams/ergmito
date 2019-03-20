test_that("boot works", {
  
  data(fivenets, package="ergmito")
  ans0 <- ergmito(fivenets ~ edges + nodematch("female"))
  ans1 <- suppressWarnings(ergmito_boot(ans0, R = 200L, ncpus = 2L))
  
  expect_output(print(ans1), "200 replicates")
  
})