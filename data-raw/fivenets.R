library(ergmito)
library(network)

# Random bernoulli graphs with random attribute age
set.seed(12312)

nets <- replicate(5, {
  network(
    rbernoulli(5, .5),
    vertex.attr = list(floor(rnorm(5, 30, 10))),
    vertex.attrnames = "age")
}, simplify = FALSE)

# Generating the samplers
examplers <- lapply(nets, function(net) {
  
  new_rergmito(net ~ edges + balance, theta = c(-2, 2))
  
})

# Generating a random sample of these networks (same size)
fivenets2 <- lapply(examplers, function(i) {
  # debug(i$sample)
  i$sample(1L, nvertex(i$network0))[[1L]]
})
ans <- ergmito(fivenets2 ~ edges + balance) 
confint(ans)
nedges(ans)

usethis::use_data(fivenets)

fivesamplers <- examplers
usethis::use_data(fivesamplers)


\frac{\exp{\transpose{\theta}\stats{y,x}}}{\transpose{\Mat{W}}\left[\exp{\transpose{\theta}\Mat{W}}\right]}