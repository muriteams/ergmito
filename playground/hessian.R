library(ergmito)

data(fivenets)
ans <- ergmito(fivenets ~ edges + ttriad + nodematch("female"))

library(microbenchmark)

microbenchmark(
  exact = with(ans$formulae, {
    exact_hessian(
      x             = target.stats,
      params = coef(ans),
      stats.weights = stats.weights,
      stats.statmat = stats.statmat
    )
  }),
  optiom = with(ans$formulae, {
    optimHess(
      par           = coef(ans),
      fn            = loglik,
      gr            = grad,
      target.stats  = target.stats,
      stats.weights = stats.weights,
      stats.statmat = stats.statmat
    )
  })
)
