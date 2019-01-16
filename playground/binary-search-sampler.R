iterator <- function(n) {
  
  # Creating new environment
  env   <- new.env()
  env$n <- n
  
  # State of the bounds
  env$ub <- n - 1L
  env$lb <- 0L
  
  # Returns true if x is even
  env$even <- function(x) {
    !(x %% 2)
  }
  
  # First position
  env$pos <- if (env$even(n)) {
    n/2L - 1L
  } else {
    (n + 1L)/2L - 1L
  }
  
  # To call when updating
  env$tell <- function() {
    # 
    # with(
    #   env,
    #   message(
    #     "The iterator is now at position: ", pos, ". ub: ", ub, " and lb: ", lb)
    #   )
    # 
  }
  
  env$tell()
  
  env$correct <- function() {
    
    if (env$pos <= 0L)
      env$pos <- 1L
    else if (env$pos >= n)
      env$pos <- env$n - 1L
    
    invisible()
    
  }
  
  # Moving one position down
  env$left <- function() {
    
    # Updating upper bound
    current_pos <- env$pos
    env$ub      <- env$pos
    
    # Computing available steps
    n <- env$pos - env$lb
    env$pos <- env$pos - max(1L, if (env$even(n)) {
      n/2L - 1L
    } else {
      (n + 1L)/2L - 1L
    })
    
    env$correct()
    env$tell()
    
    invisible(env$pos != current_pos)
    
  }
  
  # Moing one position up
  env$right <- function() {
    
    # Updating upper bound
    current_pos <- env$pos
    env$lb      <- env$pos
    
    n <- env$ub - env$pos
    
    env$pos <- env$pos + max(1L, if (env$even(n)) {
      n/2L - 1L
    } else {
      (n + 1L)/2L - 1L
    })
    
    env$correct()
    env$tell()
    
    invisible(env$pos != current_pos)
    
  }
  
  env
  
}

sampler <- function(n, w) {
  
  # Cumulative sum weights
  W <- cumsum(w)
  W <- W/(W[length(W)] + 1)
  
  # Sampling
  R   <- runif(n)
  ans <- vector("double", n)
  
  # Matching
  for (i in seq_along(R)) {
    
    # Initializing iterator
    iter <- iterator(length(W))
    
    while (TRUE) {
      
        # Is it in between
      diff_right <- W[iter$pos + 1L] - R[i]
      diff_left  <- R[i] - W[iter$pos]
      
      # Score!
      if (diff_left > 0 & diff_right > 0) {
        
        if (diff_left < diff_right)
          ans[i] <- iter$pos - 1L
        else
          ans[i] <- iter$pos
        break
        
      }  else if (diff_left > 0 & diff_right < 0) {
        # message("moving to the right. r:", R[i], " W[iter$pos + 1]:", W[iter$pos + 1])
        change <- iter$right()
      } else if (diff_left < 0 & diff_right > 0) {
        # message("moving to the left. r:", R[i], " W[iter$pos + 1]:", W[iter$pos + 1])
        change <- iter$left()
      } 
      
      if (!change) {
        ans[i] <- if (abs(diff_left) > abs(diff_right)) 
          iter$pos 
        else
          iter$pos - 1L
        break
      } 
      
      
      
    }
    
  }
  
  list(
    x = ans,
    r = R
  )
  
}

x <- iterator(11)
for (i in 1:5)
  x$left()

x <- iterator(11)
for (i in 1:5)
  x$right()

set.seed(1397123)
sol <- sampler(10000, c(1,1,1,1,1,1,1))

do.call(cbind, sol)[which(sol$r > .5 & sol$x < 2),]


table(sol$x)

set.seed(1)
do.call(cbind, sol<-sampler(10, c(1,1,1,1)))
table(sol$x)
