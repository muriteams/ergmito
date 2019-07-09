context("ERGMito GOF")

test_that("As expected 1", {
  
  set.seed(12344)
  # suppressWarnings(RNGversion("3.4.0"))
  net <- rbernoulli(4, .3)
  model <- ergmito(net ~ edges)
  gof0  <- gof_ergmito(model)
  
  # Alternative calculation
  prb <- exact_loglik(
    x             = model$formulae$stats.statmat[[1]],
    params        = coef(model),
    stats.weights = model$formulae$stats.weights[[1]],
    stats.statmat = model$formulae$stats.statmat[[1]]
    )
  prb <- exp(prb)
  
  # Adds up to one
  expect_equal(sum(prb*model$formulae$stats.weights[[1]]), 1)
  min_max <- range(model$formulae$stats.statmat[[1]])
  
  # Quantiles
  ord    <- order(model$formulae$stats.statmat[[1]])
  cumprb <- cumsum(prb[ord]*model$formulae$stats.weights[[1]][ord])
  
  # 90% CI
  idx <- c(max(which(cumprb < .05)), min(which(cumprb > .95)))
  expect_equivalent(
    gof0$ci[[1]][, c("lower-q", "upper-q")],
    model$formulae$stats.statmat[[1]][ord][idx]
    )
  expect_equivalent(gof0$ci[[1]][, c("lower-p", "upper-p")], cumprb[idx])
  
  
})
  

test_that("As expected 1", {
  
  set.seed(1244)
  # suppressWarnings(RNGversion("3.4.0"))
  net <- rbernoulli(4, .3)
  model <- ergmito(net ~ edges + ttriad)
  gof0  <- gof_ergmito(model)
  
  # Alternative calculation
  prb <- exact_loglik(
    x             = model$formulae$stats.statmat[[1]],
    params        = coef(model),
    stats.weights = model$formulae$stats.weights[[1]],
    stats.statmat = model$formulae$stats.statmat[[1]]
  )
  prb <- exp(prb)
  
  # Adds up to one
  expect_equal(sum(prb*model$formulae$stats.weights[[1]]), 1)
  min_max <- range(model$formulae$stats.statmat[[1]])
  
  # Quantiles
  for (k in 1L:2L) {
    
    ord <- order(model$formulae$stats.statmat[[1]][, k])
    S   <- sort(unique(model$formulae$stats.statmat[[1]][ord, k]))
    
    # Aggregating probs
    cumprb <- numeric(length(S))
    for (i in seq_along(S)) {
      
      ids       <- which(model$formulae$stats.statmat[[1]][,k] == S[i])
      cumprb[i] <- sum(model$formulae$stats.weights[[1]][ids]*prb[ids])
    }
    
    cumprb <- cumsum(cumprb)
    
    # 90% CI
    idx <- c(max(1, which(cumprb < .05)), min(which(cumprb > .95)))
    expect_equivalent(
      gof0$ci[[1]][, c("lower-q", "upper-q")][k,],
      S[idx]
    )
    expect_equivalent(
      gof0$ci[[1]][k, c("lower-p", "upper-p")],
      cumprb[idx])
  }
  
  
  
})
