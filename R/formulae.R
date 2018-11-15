lergm_formulae <- function(model, env = parent.frame()) {
  
  # Capturing model
  model <- sys.call()[[2]]

  # What is the first component
  LHS <- eval(model[[2]], envir = env)

  if (inherits(LHS, "list")) {
    
    env. <- env
    
    # Creating one model per model
    F. <- lapply(model[[2]][-1], function(lhs) {
      
      model[[2]] <- lhs
      m <- as.formula(model)
      expr <- parse(
        text = sprintf("lergm_formulae(model = %s, env = env.)", deparse(m))
        )
      eval(expr)
      
    })
    
    # Additive loglike function
    function(params, weights, stats) {
      
      ll <- 0
      for (i in seq_along(f.)) {
        ll <- ll + f.[[i]](params, weights[[i]], stats[[i]])
        
      }
      
      ll
      
    }
    
    
  } else if (inherits(LHS, "matrix")) {
    
    # Network size
    n <- nrow(LHS)
    
      function(params., weights., stats.) {
        exact_loglik(params., weights., stats.)
      }

  } else 
    stop("One of the components is not a matrix `", deparse(model[[2]]),
         "` is of class ", class(LHS), ".", call. = FALSE)
  
  
}

# library(sna)
# 
# set.seed(1)
# x <- rgraph(4)
# y <- rgraph(5)
# 
# L <- lergm_formulae(list(x, y) ~ edges) 
# 
# library(ergm)
# ergm.allstats(x ~ edges)
