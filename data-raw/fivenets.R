library(ergmito)
library(network)

# Random bernoulli graphs with random attribute age
set.seed(12312)

nets <- replicate(5, {
  network(
    rbernoulli(4, .5),
    vertex.attr = list(rpois(4, 20)),
    vertex.attrnames = "age")
}, simplify = FALSE)

# Generating the samplers
examplers <- lapply(nets, function(net) {
  
  new_rergmito(net ~ edges + nodeicov("age"), theta = c(-4, .2))
  
})

# Generating a random sample of these networks (same size)
fivenets <- lapply(examplers, function(i) i$sample(1L, nvertex(i$network0))[[1L]])
ans <- ergmito(fivenets ~ edges + nodeicov("age")) 
confint(ans)

usethis::use_data(fivenets)

fivesamplers <- examplers
usethis::use_data(fivesamplers)


