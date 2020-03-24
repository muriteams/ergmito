# Testing how to handle formulas -----------------------------------------------

available_terms <- function(...) {
  
  # Reading the ERGM manual
  ergm_terms <- tools::Rd_db("ergm", ...)$"ergm-terms.Rd"
  ergm_terms <- as.character(ergm_terms)
  
  # Extracting list of terms (by the alist) and removing 
  # which terms matter.
  setdiff(
    c("edges", ergm_terms[which(grepl("[\\]alias", ergm_terms)) + 2]),
    c("ergm-terms", "ergm.terms", "terms.ergm", "InitErgmTerm")
  )
}

ERGM_TERMS <- available_terms()



# formulae2 <- function(x) {
  
  x <- net ~ edges + nodematch('asda') + ndoe1match('1', 1231) - I(1/nnodes(net)) + offset(edges)

  # First step, capture all ergm terms, these should be one of the following:
  # - name
  # - name.w.dot
  # - name_w_bar
  # - name_w_num
  # Names can also include an argument (or multiple arguments), so a function call
  
  # But before that, we need to identify special characters, like 
  # 'offset', 'I', 'factor' (for fixed effect), and re (for random effects).

  # Step 1: Chunk it off
  x_terms <- rownames(attr(terms(x), "factors"))
  
  # Step 2: Get the offset terms -----------------------------------------------
  
  # Which is offset
  x_offset <- grep("^offset[(]", x_terms)
  
  if (length(x_offset) > 1)
    stop(
      "There can be only one offset term: ",
      paste(x_terms[x_offset], collapse = ", "),
      ".", call. = FALSE
    )
  
  # Is it part of the ergm terms
  offset_ergm_term <- gsub("^offset[(]([a-zA-Z0-9_.]+).+[)]$", "\\1", x_terms[x_offset])
  if (offset_ergm_term %in% ERGM_TERMS)
    x_ergm_terms <- x_offset
  else
    x_ergm_terms <- NULL
    
  
  # Step 2: Get the ergm e.q to get stats counts ---------------------------------
  x_ergm_terms <- which(gsub("[(].+", "", x_terms) %in% ERGM_TERMS)
  
  x_ergm <- update.formula(
    x, paste(". ~", paste(x_terms[x_ergm_terms], collapse = " + "))
    )
  
  # Pass this to the old-fashion way to generate the counts
  # ... #
  
  
  # Step 3: 
  
  x_special_removed <- gsub("\\s*(+|-)?\\s*(offset|I)[(].+[)]\\s*(+|-)?\\s*", "", x)
  
  
  gsub(
    "[a-zA-Z_\\.][a-zA-Z0-9_\\.]*([(][^()]*[)])?",
    "yes",
    "edges + nodematch('asda') + ndoe1match('1', 1231) - I(1/nnodes(net))",
    perl = TRUE
    )
  
# }
# 
# # library(ergmito)
# data("fivenets")
# net <- fivenets
# 
# # Simple 
# f0 <- net ~ edges
# formulae2(f0)
# 
# # Bit more complex
# f1 <- net ~ edges + nodematch("female")
# formulae2(f1)
# 
# # Hard 1
# f2 <- net ~ edges + nodematch("female") | edges + nodematch("female") + log(nvertex(net))
# formulae2(f2)
# 
# # Hard 2
# f2 <- net ~ edges + nodematch("female") + log(nvertex(net)) |~ edge
# formulae2(f2)
# 
# f2 <- ergmito_formulae(net ~ edges + ttriad, net ~ edges + log(edges) + ttriad)
