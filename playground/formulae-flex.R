# Testing how to handle formulas -----------------------------------------------

formulae2 <- function(x) {
  
  # First step, capture all ergm terms, these should be one of the following:
  # - name
  # - name.w.dot
  # - name_w_bar
  # - name_w_num
  # Names can also include an argument (or multiple arguments), so a function call
  
  x <- net ~ edges + nodematch('asda') + ndoe1match('1', 1231) - I(1/nnodes(net)) + offset(edges)
  trsm <- terms(x)
  as.formula(trsm)
  
  # But before that, we need to identify special characters, like 
  # 'offset' and 'I'.
  x_special_removed <- gsub("\\s*(+|-)?\\s*(offset|I)[(].+[)]\\s*(+|-)?\\s*", "", x)
  
  
  gsub(
    "[a-zA-Z_\\.][a-zA-Z0-9_\\.]*([(][^()]*[)])?",
    "yes",
    "edges + nodematch('asda') + ndoe1match('1', 1231) - I(1/nnodes(net))",
    perl = TRUE
    )
  
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
