library(ergmito)
library(ergm)

set.seed(1133)
nets <- lapply(1:40, function(i) {
  net <- rbernoulli(4, .3)
  net <- network::network(net)
  network::set.vertex.attribute(net, "age", runif(4))
})


# Should be the same
tmp0 <- tempfile()
tmp1 <- tempfile()
Rprofmem(tmp0, threshold = 1024)
ans0 <- ergmito_formulae(nets ~ edges + nodeicov("age") + ttriad)
Rprofmem(NULL)
Rprofmem(tmp1, threshold = 1024)
ans1 <- ergmito_formulae(nets ~ edges + nodeicov("age") + ttriad, use_ptr = TRUE)
Rprofmem(NULL)

# Iterative calls
ncalls <- 500
v0 <- vector("list", ncalls)
v1 <- vector("list", ncalls)
for (i in 1:ncalls) {
  p <- runif(ncol(ans0$target.stats))
  v0[[i]]$l <- ans0$loglik(
    p,
    stats.weights = ans0$stats.weights,
    stats.statmat = ans0$stats.statmat,
    target.stats = ans0$target.stats
    )
  
  v0[[i]]$g <- ans0$grad(
    p,
    stats.weights = ans0$stats.weights,
    stats.statmat = ans0$stats.statmat,
    target.stats = ans0$target.stats
  )
  
  v1[[i]]$l <- ans1$loglik(p)
  
  v1[[i]]$g <- ans1$grad(p)
}

bench::mark(
  old = ans0$loglik(
    p,
    stats.weights = ans0$stats.weights,
    stats.statmat = ans0$stats.statmat,
    target.stats = ans0$target.stats
  ),
  new = ans1$loglik(p)
)

range(abs(sapply(v0, "[[", "l") - sapply(v1, "[[", "l")))
s0 <- t(do.call(cbind, lapply(v0, "[[", "g")))
s1 <- t(do.call(cbind, lapply(v1, "[[", "g")))
range(abs(s0 - s1))

# Process should be the same
ans0 <- ergmito(nets ~ edges, use_ptr = FALSE)
ans1 <- ergmito(nets ~ edges, use_ptr = TRUE)

p0 <- profvis::profvis(ans0 <- ergmito(nets ~ edges + absdiff("age") + ttriad, use_ptr = FALSE))
p1 <- profvis::profvis(ans1 <- ergmito(nets ~ edges + absdiff("age") + ttriad, use_ptr = TRUE))
p2 <- profvis::profvis(ans1 <- ergm_blockdiag(nets ~ edges + absdiff("age") + ttriad))

library(microbenchmark)
ans <- microbenchmark(
  old = coef(ergmito(nets ~ edges + absdiff("age") , use_ptr = FALSE)),
  new = coef(ergmito(nets ~ edges + absdiff("age") , use_ptr = TRUE)), times = 10
)

library(bench)
ans <-mark(
  old = ergmito(nets ~ edges + absdiff("age"), use_ptr = FALSE),
  new = ergmito(nets ~ edges + absdiff("age"), use_ptr = TRUE),
  iterations = 4, check = FALSE
)

experiment <- function() {
  net <- rbernoulli(sample(3:4, size = sample(c(10, 20, 40), 1), replace = TRUE), p = .15)
  suppressWarnings({
    ans0 <- ergmito(net ~ edges + ttriad + mutual, keep.stats = FALSE)
    ans1 <- ergmito(net ~ edges + ttriad + mutual, keep.stats = FALSE, use_ptr = TRUE)
  })
  
  data.frame(
    par_old1 = ans0$optim.out$par[1],
    par_old2 = ans0$optim.out$par[2],
    par_old3 = ans0$optim.out$par[3],
    par_new1 = ans1$optim.out$par[1],
    par_new2 = ans1$optim.out$par[2],
    par_new3 = ans1$optim.out$par[3],
    time_old = ans0$timer["total"],
    time_new = ans1$timer["total"],
    counts_f_old = ans0$optim.out$counts["function"],
    counts_f_new = ans1$optim.out$counts["function"]
  )
  
}

set.seed(131)
ans <- parallel::mclapply(1:1000, function(i) experiment(), mc.cores = 4L)
ans <- do.call(rbind, ans)

hist(ans$counts_f_old - ans$counts_f_new)
hist(as.numeric(ans$time_old - ans$time_new))
