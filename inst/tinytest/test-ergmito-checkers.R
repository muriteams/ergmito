
# Fully connected network
x <- rbernoulli(c(4,4,4), 1)
ans0 <- ergmito(x ~ edges + ttriad)

# Empty graph
x <- rbernoulli(c(4,4,4), 0)
ans1 <- ergmito(x ~ edges + ttriad)

# Very low density
x <- lapply(x, function(x.) {
  x.[2:3] <- 1L
  x.
})
ans2 <- ergmito(x ~ edges + ttriad)

expect_equivalent(coef(ans0), c(Inf, Inf))
expect_equivalent(coef(ans1), -c(Inf, Inf))
expect_equivalent(coef(ans2)[2], -Inf)
expect_equivalent(summary(ans0)$coefs[, "Pr(>|z|)"], c(0, 0), tol = 1e-2)
expect_equivalent(summary(ans1)$coefs[, "Pr(>|z|)"], c(0, 0), tol = 1e-2)
expect_equivalent(summary(ans2)$coefs[, "Pr(>|z|)"], c(0, 0), tol = 1e-2)
