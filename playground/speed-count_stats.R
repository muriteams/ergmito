library(ergmito)
library(ergm)
library(network)
library(magrittr)

set.seed(121)

nets <- replicate(100, {
  network(rbernoulli(5), vertex.attr = list(a=rpois(5, 4)), vertex.attrnames = "a")
}, simplify = FALSE)

bench::mark(
  as.matrix(nets[[1]]),
  as.adjmat(nets[[1]]), check = FALSE
) %>% plot


count_stats_ergm <- function(m) {
  LHS <- eval(m[[2]], envir = environment(m))
  m. <- update.formula(m, LHS[[i]] ~ .)
  environment(m.) <- environment()
  
  ans <- NULL
  for (i in seq_along(LHS))
    ans <- rbind(ans, summary(m.))
  
  ans
}

ans0 <- count_stats(nets ~ edges + nodeicov("a"))
ans1 <- count_stats_ergm(nets ~ edges + nodeicov("a"))

b <- bench::mark(
  ERGMito = count_stats(nets ~ edges + nodeicov("a")),
  ergm    = count_stats_ergm(nets ~ edges + nodeicov("a")), check = FALSE
)

microbenchmark::microbenchmark(
  ERGMito = count_stats(nets ~ edges + nodeicov("a")),
  ergm    = count_stats_ergm(nets ~ edges + nodeicov("a"))
)

plot(b)

z <- new_rergmito(nets[[1]] ~ edges + mutual + nodeicov("a"))
