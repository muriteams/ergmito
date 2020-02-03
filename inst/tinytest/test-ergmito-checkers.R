
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
  

