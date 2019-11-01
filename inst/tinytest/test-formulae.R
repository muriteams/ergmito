
data("fivenets")
net0 <- fivenets
net0[[4]] <- USArrests

expect_error(ergmito_formulae(net0 ~ edges), "not a matrix")
expect_error(ergmito_formulae(USArrests ~ edges), "networks or a single")
  

data("fivenets")
ans <- ergmito(fivenets ~ edges)
expect_error(exact_loglik(
  x             = ans$formulae$target.stats,
  params        = ans$coef,
  stats.weights = ans$formulae$stats.weights[-2],
  stats.statmat = ans$formulae$stats.statmat 
  ), "weights.+should match"
)

expect_error(exact_loglik(
  x             = ans$formulae$target.stats,
  params        = ans$coef,
  stats.weights = ans$formulae$stats.weights,
  stats.statmat = ans$formulae$stats.statmat[-2] 
), "stats[.]statmat.+should match"
)

expect_error(exact_loglik(
  x             = ans$formulae$target.stats[-(1:5),,drop=FALSE],
  params        = ans$coef,
  stats.weights = ans$formulae$stats.weights,
  stats.statmat = ans$formulae$stats.statmat 
), "observed stat"
)


ans0 <- ans$formulae$grad(
  target.stats  = ans$formulae$target.stats,
  params        = ans$coef,
  stats.weights = ans$formulae$stats.weights,
  stats.statmat = ans$formulae$stats.statmat   
)

expect_equal(ans0[1], 0.0, tol=1e-5)

