library(microbenchmark)
library(lergm)

nets <- rbernoulli(5)

mysummary <- function(model) {
  
  net <- eval(model[[2]], envir = environment(model))
  net <- network::as.network(net)
  
  summary(
    ergm::ergm_model(model, net, response = NULL, role = "target"),
    net,
    response = NULL
  )
  
}

cmp <- compiler::cmpfun(mysummary)

microbenchmark(
  ergm::summary.formula(nets ~ edges + mutual),
  ergm::summary_formula(nets ~ edges + mutual),
  mysummary(nets ~ edges + mutual),
  cmp(nets ~ edges + mutual),
  times = 1000, unit="relative"
)


