powerset <- function(n) {
  
  set   <- 1:(n*(n-1))
  
  sets <- NULL
  for (i in set) {
    for (s in sets)
      sets <- c(sets, list(c(s, i)))
    sets <- c(sets, list(i))
    # print(sets)
  }
  
  c(sets, list(NULL))
  
}

ans <- powerset(4)
