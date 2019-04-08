context("ERGMito GOF")

test_that("As expected", {
  
  set.seed(12344)
  RNGversion("3.4.0")
  net <- rbernoulli(3, .3)
  model <- ergmito(net ~ edges)
  gof0  <- gof_ergmito(model)
  
  # Alternative calculation
  dat <- model$formulae$stats
  prb <- exact_loglik(dat[[1]]$statmat, coef(model), dat[[1]]$weights, dat[[1]]$statmat)
  prb <- exp(prb)
  
  # Adds up to one
  expect_equal(sum(prb*dat[[1]]$weights), 1)
  min_max <- range(dat[[1]]$statmat)
  
  # Quantiles
  ord    <- order(dat[[1]]$statmat)
  cumprb <- cumsum(prb[ord]*dat[[1]]$weights[ord])
  
  # 90% CI
  expect_equivalent(gof0$ci[[1]][, c("lower-q", "upper-q")], dat[[1]]$statmat[ord][c(1, 4)])
  expect_equivalent(gof0$ci[[1]][, c("lower-p", "upper-p")], cumprb[c(1, 4)])
  RNGversion("3.5.0")
  
})