library(ergmito)
library(ergm)

sampler <- new_rergmito(rbernoulli(4) ~ edges + ctriad, theta = c(-1, 1))

set.seed(1)
nsim <- 500
nets <- vector("list", nsim)
pars <- replicate(nsim, runif(2, -2, 2), simplify = FALSE)

for (i in seq_along(nets))
  nets[[i]] <- sampler$sample(30, s=4, theta = pars[[i]])

# Mode 0: Using 5 maxtries -----------------------------------------------------
ans_0 <- vector("list", nsim)
for (i in seq_along(ans_1)) {
  tmp <- tryCatch(ergmito(nets[[i]] ~ edges + ctriad, ntries = 1L), error = function(e) e)
  
  if (inherits(tmp, "error")) {
    ans_1[[i]] <- list(status=30)
    next
  }
  # if (tmp$status != 0) {
  #   message("An error!")
  # }
  
  ci <- confint(tmp)
  
  ans_0[[i]] <- list(
    coef   = coef(tmp),
    status = tmp$status ,
    best   = tmp$best_try,
    covered = (ci[,1] <= pars[[1]]) & (pars[[1]] <= ci[,2]),
    ll      = logLik(tmp),
    hist    = tmp$history
  )
}

# Mode 1: Using 5 maxtries -----------------------------------------------------
ans_1 <- vector("list", nsim)
for (i in seq_along(ans_1)) {
  tmp <- tryCatch(ergmito(nets[[i]] ~ edges + ctriad, ntries = 20L), error = function(e) e)
  
  if (inherits(tmp, "error")) {
    ans_1[[i]] <- list(status=30)
    next
  }
  # if (tmp$status != 0) {
  #   message("An error!")
  # }
  
  ci <- confint(tmp)
  
  ans_1[[i]] <- list(
    coef   = coef(tmp),
    status = tmp$status ,
    best   = tmp$best_try,
    covered = (ci[,1] <= pars[[1]]) & (pars[[1]] <= ci[,2]),
    ll      = logLik(tmp),
    hist    = tmp$history
  )
}

# Mode 2: A single attampt and using MPLE --------------------------------------
ans_2 <- vector("list", nsim)
for (i in seq_along(ans_1)) {
  
  net_tmp <- blockdiagonalize(nets[[i]])
  net_tmp <- suppressMessages(
    ergmMPLE(net_tmp ~ edges + ctriad, output = "fit")
    )
  
  net_tmp <- coef(net_tmp) 
  net_tmp[!is.finite(net_tmp)] <- sign(net_tmp[!is.finite(net_tmp)])*10
  
  tmp <- tryCatch(
    ergmito(nets[[i]] ~ edges + ctriad, init = net_tmp, ntries = 1L),
    error = function(e) e
    )
  
  if (inherits(tmp, "error")) {
    ans_2[[i]] <- list(status=30)
    next
  }
  
  ci <- confint(tmp)
  
  ans_2[[i]] <- list(
    coef   = coef(tmp),
    status = tmp$status ,
    best   = tmp$best_try,
    covered = (ci[,1] <= pars[[1]]) & (pars[[1]] <= ci[,2]),
    ll      = logLik(tmp),
    hist    = tmp$history
  )
}

# Mode 3: A single attampt and using MPLE --------------------------------------
ans_3 <- vector("list", nsim)
for (i in seq_along(ans_1)) {
  
  net_tmp <- blockdiagonalize(nets[[i]])
  net_tmp <- suppressMessages(
    ergmMPLE(net_tmp ~ edges + ctriad, output = "fit")
  )
  
  net_tmp <- coef(net_tmp) 
  net_tmp[!is.finite(net_tmp)] <- sign(net_tmp[!is.finite(net_tmp)])*10
  
  tmp <- tryCatch(
    ergmito(nets[[i]] ~ edges + ctriad, init = net_tmp, ntries = 20L),
    error = function(e) e
  )
  
  if (inherits(tmp, "error")) {
    ans_3[[i]] <- list(status=30)
    next
  }
  
  ci <- confint(tmp)
  
  ans_3[[i]] <- list(
    coef   = coef(tmp),
    status = tmp$status ,
    best   = tmp$best_try,
    covered = (ci[,1] <= pars[[1]]) & (pars[[1]] <= ci[,2]),
    ll      = logLik(tmp),
    hist    = tmp$history
  )
}


# Checking bias ----------------------------------------------------------------
true_coef <- do.call(rbind, pars)

good <- which(
  !sapply(lapply(ans_0, "[[", "coef"), is.null) &
    !sapply(lapply(ans_1, "[[", "coef"), is.null) & 
    !sapply(lapply(ans_2, "[[", "coef"), is.null) &
    !sapply(lapply(ans_3, "[[", "coef"), is.null)
)

bias_0 <- do.call(rbind, lapply(ans_1[good], "[[", "coef")) - true_coef[good, ]
bias_1 <- do.call(rbind, lapply(ans_1[good], "[[", "coef")) - true_coef[good, ]
bias_2 <- do.call(rbind, lapply(ans_2[good], "[[", "coef")) - true_coef[good, ]
bias_3 <- do.call(rbind, lapply(ans_3[good], "[[", "coef")) - true_coef[good, ]

summary(bias_0)
summary(bias_1)
summary(bias_2)
summary(bias_3)

mean(abs(bias_1) > abs(bias_2))
mean(abs(bias_1) > abs(bias_3))
mean(abs(bias_3) > abs(bias_2))
mean(abs(bias_0) > abs(bias_2))

# Differences in biases
finite <- is.finite(bias_0) & is.finite(bias_2)
t.test(as.vector(abs(bias_1))[finite], as.vector(abs(bias_2))[finite])

# Checking coverage ------------------------------------------------------------
colMeans(do.call(rbind, lapply(ans_1, "[[", "covered")), na.rm = TRUE)
colMeans(do.call(rbind, lapply(ans_2, "[[", "covered")), na.rm = TRUE)
colMeans(do.call(rbind, lapply(ans_3, "[[", "covered")), na.rm = TRUE)

# Checking how many steps ------------------------------------------------------
print(cumsum(prop.table(table(unlist(sapply(ans_1, "[[", "best"))))), digits=2)
print(cumsum(prop.table(table(unlist(sapply(ans_2, "[[", "best"))))), digits=2)
print(cumsum(prop.table(table(unlist(sapply(ans_3, "[[", "best"))))), digits=2)

# Checking likelihoods
mean(sapply(ans_3[good], "[[", "ll") > sapply(ans_2[good], "[[", "ll"), na.rm = TRUE)
mean(sapply(ans_1[good], "[[", "ll") > sapply(ans_2[good], "[[", "ll"), na.rm = TRUE)


# Are we doing better? ---------------------------------------------------------

