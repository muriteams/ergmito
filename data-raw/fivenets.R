library(ergmito)
library(network)

# Random bernoulli graphs with random attribute age
set.seed(1878)
n <- 4
nets <- replicate(5, {
  network(
    rbernoulli(n, .5),
    vertex.attr = list(sample(c(0,1), n, TRUE)),
    vertex.attrnames = "female")
}, simplify = FALSE)

# Generating the samplers
examplers <- lapply(nets, function(net) {
  
  new_rergmito(net ~ edges + nodematch("female"), theta = c(-2, 2))
  
})

# Generating a random sample of these networks (same size)
fivenets <- lapply(examplers, function(i) {
  # debug(i$sample)
  i$sample(1L, nvertex(i$network0))[[1L]]
})
ans <- ergmito(fivenets ~ edges + nodematch("female")) 
confint(ans)
nedges(ans)

usethis::use_data(fivenets, overwrite = TRUE)

fivesamplers <- examplers
usethis::use_data(fivesamplers, overwrite = TRUE)


# \frac{\exp{\transpose{\theta}\stats{y,x}}}{\transpose{\Mat{W}}\left[\exp{\transpose{\theta}\Mat{W}}\right]}