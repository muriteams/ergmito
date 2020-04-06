
data("fivenets")
net0 <- fivenets
net0[[4]] <- USArrests

expect_error(ergmito_formulae(net0 ~ edges), "not a matrix")
expect_error(ergmito_formulae(USArrests ~ edges), "networks or a single")
  

data("fivenets")
ans <- ergmito(fivenets ~ edges)
expect_error(exact_loglik(
  x             = ans$formulae$target_stats,
  params        = ans$coef,
  stats_weights = ans$formulae$stats_weights[-2],
  stats_statmat = ans$formulae$stats_statmat 
  ), "should be numeric"
)

expect_error(exact_loglik(
  x             = ans$formulae$target_stats,
  params        = ans$coef,
  stats_weights = ans$formulae$stats_weights,
  stats_statmat = ans$formulae$stats_statmat[-2] 
), "should be matrix"
)

expect_error(exact_loglik(
  x             = ans$formulae$target_stats[-(1:5),,drop=FALSE],
  params        = ans$coef,
  stats_weights = ans$formulae$stats_weights,
  stats_statmat = ans$formulae$stats_statmat 
), "greater than 0"
)


ans0 <- ans$formulae$grad(params = ans$coef)
expect_equal(ans0[1], 0.0, tol=1e-5)

data(fivenets)
ans <- ergmito(fivenets ~ edges + mutual + ttriad)

ans0 <- ans$formulae$grad(params = rep(0, length(coef(ans))))
ans1 <- exact_gradient(
  x = ans$formulae$target_stats,
  params        = rep(0, length(coef(ans))),
  stats_weights = ans$formulae$stats_weights,
  stats_statmat = ans$formulae$stats_statmat
)

expect_equal(ans0, ans1)

# Checking offset terms -----------------------------------
data(fivenets)

m0 <- ergmito_formulae(fivenets ~ edges + ttriad, model_update = ~ . + offset(edges))
m1 <- ergmito_formulae(fivenets ~ edges + ttriad, model_update = ~ . + I(edges))

expect_equal(
  m0$loglik(c(.5, .5)),
  m1$loglik(c(.5, .5, 1))
)

expect_equal(
  m0$grad(c(.5, .5)),
  m1$grad(c(.5, .5, 1))[-3,,drop=FALSE]
)
