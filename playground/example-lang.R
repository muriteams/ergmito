library(ggplot2)
library(dplyr)

data("diamonds")

draw_points <- function(dat, var) {
  
  mycall <- match.call()
  
  colnames(dat)[colnames(dat) == var] <- "var"
  
  ggplot(dat, aes(x = cut, y = var)) +
    ylab(var) + 
    geom_point()
  
}

draw_points(diamonds, var="color")
