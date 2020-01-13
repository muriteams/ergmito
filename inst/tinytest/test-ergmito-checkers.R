
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

if (.Platform$OS.type == "unix") {

  expect_equivalent(coef(ans0), c(Inf, Inf))
  expect_equivalent(coef(ans1), -c(Inf, Inf))
  expect_equivalent(coef(ans2)[2], -Inf)
  expect_equivalent(summary(ans0)$coefs[, "Pr(>|z|)"], c(0, 0), tol = 1e-2)
  expect_equivalent(summary(ans1)$coefs[, "Pr(>|z|)"], c(0, 0), tol = 1e-2)
  expect_equivalent(summary(ans2)$coefs[, "Pr(>|z|)"], c(0, 0), tol = 1e-2)

} else {
  
  pretty_printer <- function(x) {
    paste(
      "\n",
      paste(capture.output(eval(x)), collapse = "\n"),
      "\n"
    )
  }
  
  message(
    "\nOn windows these are not easy to test. So here is the expected results:\n",
    pretty_printer(summary(ans0)), "\n----- should be c(Inf, Inf) \n----- with pval 0.\n",
    pretty_printer(summary(ans1)), "\n----- should be -c(Inf, Inf), \n----- with pval 0.\n",
    pretty_printer(summary(ans2)), "\n----- should be -Inf (2nd) \n----- with pval 0."
    )
  
}
