set.seed(1)
net <- rbernoulli(rep(5, 5))
ans <- ergmito(net ~ edges + mutual)

expect_silent(plot(ans))
  
