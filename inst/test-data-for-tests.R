# if (file.exists("test-data-for-tests.rda"))
#   skip("ajaja")
# 
# 
# set.seed(76183)
# net <- matrix(rbinom(16, 1, .5), ncol=4)
# diag(net) <- 0
# library(ergm)
# suppressMessages(
#   ans0 <- ergm(net ~ mutual + edges)
# )
# 
# save.image("test-data-for-tests.rda")