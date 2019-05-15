context("ERGMito GOF")

test_that("As expected 1", {
  
  set.seed(12344)
  # suppressWarnings(RNGversion("3.4.0"))
  net <- rbernoulli(4, .3)
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
  idx <- c(max(which(cumprb < .05)), min(which(cumprb > .95)))
  expect_equivalent(
    gof0$ci[[1]][, c("lower-q", "upper-q")],
    dat[[1]]$statmat[ord][idx]
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
  dat <- model$formulae$stats
  prb <- exact_loglik(dat[[1]]$statmat, coef(model), dat[[1]]$weights, dat[[1]]$statmat)
  prb <- exp(prb)
  
  # Adds up to one
  expect_equal(sum(prb*dat[[1]]$weights), 1)
  min_max <- range(dat[[1]]$statmat)
  
  # Quantiles
  for (k in 1L:2L) {
    
    ord <- order(dat[[1]]$statmat[, k])
    S   <- sort(unique(dat[[1]]$statmat[ord, k]))
    
    # Aggregating probs
    cumprb <- numeric(length(S))
    for (i in seq_along(S)) {
      
      ids       <- which(dat[[1]]$statmat[,k] == S[i])
      cumprb[i] <- sum(dat[[1]]$weights[ids]*prb[ids])
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
