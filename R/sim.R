new_lergmsim <- function(model, theta = NULL, env = parent.frame(), ...) {
  
  # Getting the estimates
  if (!length(theta))
    theta <- coef(lergm(model))
  
  # Obtaining the network(s) object
  net <- eval(model[[2]])
  
  # Checking stats0
  dots <- list(...)
  if (length(dots$zeroobs) && (dots$zeroobs & nnets(net) > 1L)) {
    warning(
      "The option `zeroobs` was set to FALSE. `zeroobs = TRUE` invalidates the pooled ERGM.",
      call. = FALSE)
    
    dots$zeroobs <- FALSE
  } 

  # Generating powersets
  PSETS <- lapply(2:5, powerset)
  
  # Computing probabilities
  model. <- stats::update.formula(model, pset ~ .)
  P <- lapply(PSETS, function(psets) {
    
    # Getting the corresponding powerset
    ans <- vector("double", length(psets))
    environment(model.) <- environment()
    
    S <- NULL
    for (i in seq_along(psets)) {
      
      # Updating the model
      pset   <- psets[[i]]
      
      # In the first iteration we need to compute the statsmat
      if (i == 1L)
        S <- do.call(ergm::ergm.allstats, c(list(formula = model.), dots))
      
      ans[i] <- exact_loglik(
        params = theta,
        stat0  = summary(model.),
        stats  = S
        )
        # ergm::ergm.exact(
        # eta     = theta,
        # formula = model.,
        # statmat = S$statmat,
        # weights = S$weights
        # )
      
    }
    
    ans
    
  })
  
  P
  
}

