
# Fully connected network
x <- rbernoulli(c(4,4,4), 1)
ans0 <- ergmito(x ~ edges + ttriad)

# Very high density
x <- lapply(x, function(x.) {
  x.[2] <- 0L
  x.
})
ans0b <- ergmito(x ~ edges + ttriad)
ans0b$formulae$loglik(c(1e3, coef(ans0b)[2]))

# Empty graph
x <- rbernoulli(c(4,4,4), 0)
ans1 <- ergmito(x ~ edges + ttriad)

# Very low density
x <- lapply(x, function(x.) {
  x.[2:3] <- 1L
  x.
})
ans2 <- ergmito(x ~ edges + ttriad)

# Theoretical limit
s <- with(ans0$formulae, {
  sapply(1:nnets(ans2), function(i) {
    
    # To ease things
    t. <- target_stats[i, ]
    th <- coef(ans0)[1]
    w <- stats_weights[[i]]
    s <- stats_statmat[[i]]
    
    # Which equal the suff stat that is inf
    loc <- which(s[,2] == t.[2])
    
    delta <- s - matrix(t., ncol = 2, nrow = nrow(s), byrow = TRUE)
    
    t_delta <- ifelse(delta[loc, 1] == 0, 0, th * delta[loc, 1])
    
    # First bit
    - sum(delta[loc, 1]^2 * w[loc] * exp(t_delta))/
      sum(w[loc] * exp(t_delta)) +
      # Second bit
      sum( delta[loc] * w[loc] * exp(t_delta))^2/
      sum(w[loc] * exp(t_delta))^2
    
  })
})
-1/sum(s)
vcov(ans2)
ans2$formulae$hess(c(coef(ans2)[1], -1e10))

expect_true(all(is.infinite(coef(ans0))))
expect_true(all(coef(ans1) < 0) & all(is.infinite(coef(ans1))))
expect_true(coef(ans2)[2] < 0 & is.infinite(coef(ans2)[2]))
expect_equal(summary(ans0)$coefs[, "Pr(>|z|)"], c(0, 0), tol = 1e-2)
expect_equal(summary(ans1)$coefs[, "Pr(>|z|)"], c(0, 0), tol = 1e-2)
expect_equal(summary(ans2)$coefs[, "Pr(>|z|)"], c(0, 0), tol = 1e-2)

pretty_printer <- function(x) {
  paste(
    "\n",
    paste(capture.output(eval(x)), collapse = "\n"),
    "\n"
  )
}

message(
  "\nOn some OSs these are not easy to test. So here is the expected results:\n",
  pretty_printer(summary(ans0)), "\n----- should be c(Inf, Inf) \n----- with pval 0.\n",
  pretty_printer(summary(ans1)), "\n----- should be -c(Inf, Inf), \n----- with pval 0.\n",
  pretty_printer(summary(ans2)), "\n----- should be -Inf (2nd) \n----- with pval 0."
  )
  

