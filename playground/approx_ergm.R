library(ergmito)
library(network)
library(ergm)

# set.seed(1213)
# mat <- as.network(rbernoulli(50, .02))
# 
# sim1 <- ergm_MCMC_sample(
#   mat,
#   model    = ergm_model(mat ~ edges + ctriad),
#   proposal = ergm_proposal(
#     ~ .,
#     arguments = ergm:::InitErgmProposal.TNT(),
#     nw = mat
#     ),
#   control = control.ergm(), verbose = TRUE, eta = c(-2/10, 2/10)
#   )
# 
# library(sna)
# set.seed(1)
# gplot(sim1$networks[[1]])

# net <- as.network(mat)
# ans_ergm <- ergm(net ~ edges + mutual, control = control.ergm(seed = 991))

ergm_fracapprox <- function(formula, R = 500L, ...) {
  
  # Dividing the data
  LHS <- eval(formula[[2]])
  
  n <- nvertex(LHS)
  n <- n - (n %% 5)
  IDX <- replicate(R, sample.int(nvertex(LHS), n), FALSE)
  
  grps <- sort(rep(1:(n %/% 5), 5))
  IDX <- lapply(IDX, split, f = grps)
  
  estimates <- vector("list", R)
  formula <- stats::update.formula(formula, lhs ~ .)
  environment(formula) <- environment()
  for (i in seq_len(R)) {
    lhs <- induced_submat(LHS, IDX[[i]])
    estimates[[i]] <- ergmito(formula, ...)
  }

  estimates
}

# set.seed(12)
# ans_ergmito <- ergm_fracapprox(mat ~ edges + mutual)

# Trying with other ergm functions

data(sampson)

model <- samplk3 ~ edges + mutual

ans_ergmito  <- ergm_fracapprox(model, R = 200)
ans_ergmito_boot  <- ergmito_boot(ans_ergmito[[1]], R = 1000)
ans_ergm     <- ergm(model)
ans_ergmMPLE <- ergmMPLE(model)
ans_ergmMPLE <- with(ans_ergmMPLE, glm(
  response ~ 0 + predictor, family = binomial("logit"),
  weights = weights))

summary(ans_ergmito)$coef
summary(ans_ergmito_boot)$coef
summary(ans_ergm)$coefficients
summary(ans_ergmMPLE)$coef
