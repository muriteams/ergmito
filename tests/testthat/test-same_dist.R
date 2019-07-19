context("Comparing networks same dist")

test_that("class network works", {
  
  data("fivenets")
  expect_true(same_dist(fivenets[[1]], fivenets[[2]]))
  expect_false(same_dist(fivenets[[1]], fivenets[[2]], "female"))
  
  expect_error(same_dist(fivenets[[1]], as_adjmat(fivenets[[1]])))
  expect_error(same_dist(fivenets[[1]], fivenets[[2]], "females"))
  
})

test_that("class matrix works", {
  
  data("fivenets")
  
  nets <- rbernoulli(c(4, 5))
  expect_true(same_dist(nets[[1]], nets[[1]]))
  expect_false(same_dist(nets[[1]], nets[[2]]))
  
  expect_error(same_dist(nets[[1]], fivenets[[1]]))
  # expect_error(same_dist(nets[[1]], nets[[2]], "females"))
  
})