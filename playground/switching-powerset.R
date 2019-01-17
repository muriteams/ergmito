

newenv <- function(n) {
  m <- 2^(n*(n-1))
  e <- list2env(list(
    x = vector("list", m),
    m = m
  ))
}

pset <- function(
  n, e = NULL, i = 1, pos = 1,
  lb = 1,
  ub = 2^(n*(n-1))
  ) {

  # En el ultimo?
  if (i == n^2) #  | pos > ub
    return(NULL)
  
  # Not in the diagonal?
  if ((i %% (n+1)) != 1)  {
    
    newub <- (ub - lb + 1)/2 + lb - 1
    
    cat(sprintf("ub: %04.1f, lb: %04.1f, newub: %04.1f, pos: %04.1f\n", ub, lb, newub, pos))

    if (pos > ub) {
      return(NULL)  
    }
    
    # Don't add anything
    if (length(e$x[[pos]]) > 0)
      e$x[[newub + 1]] <- e$x[[pos + 1]]

    # Add the i-th edge
    e$x[[pos + 1]] <- c(e$x[[pos]], i)

  } else
    return(pset(n, e, i + 1, pos, lb = lb, ub = ub))
    
  # Continue in the recursion
  c(
    pset(n = n, e = e, i = i + 1, pos = pos + 1, lb = lb, ub = newub),
    pset(n = n, e = e, i = i + 1, pos = newub + 1, lb = newub + 1, ub = ub)
    )
  
}

e <- newenv(3)
z <- pset(3, e)
# e$x
e$x <- lapply(e$x, unique)
length(unique(e$x))


f <- function(n, i = 1, x = NULL) {

  if (i == n)
    return(list(x))
  
  if ((i %% (sqrt(n)+1)) != 1)
    c(f(n, i + 1, x), f(n, i + 1, c(x, i)))
  else
    f(n, i + 1, x)
  
}

z <- f(9)
