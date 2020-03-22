# Testing how to handle formulas -----------------------------------------------

formulae2 <- function(x) {
  
  # Catch complex conditionals
  strsplit(deparse(x), split = "[|]\\s*[~]")
  
}

# library(ergmito)
data("fivenets")
net <- fivenets

# Simple 
f0 <- net ~ edges
formulae2(f0)

# Bit more complex
f1 <- net ~ edges + nodematch("female")
formulae2(f1)

# Hard 1
f2 <- net ~ edges + nodematch("female") | edges + nodematch("female") + log(nvertex(net))
formulae2(f2)

# Hard 2
f2 <- net ~ edges + nodematch("female") + log(nvertex(net)) |~ edge
formulae2(f2)

f2 <- ergmito_formulae(net ~ edges + ttriad, net ~ edges + log(edges) + ttriad)
