#' Count Network Statistics
#' 
#' This function is similar to what [ergm::summary_formula] does, but it provides
#' a fast wrapper suited for matrix class objects.
#' @param X List of square matrices. (networks)
#' @param terms Character vector with the names of the statistics to calculate.
#' Currently, the only available statistics are: '\Sexpr{paste(ergmito::count_available(), collapse="', '")}'.
#' @export
#' @return A matrix of size `length(X) * length(terms)` with the corresponding
#' counts of statistics.
#' @examples 
#' # DGP 
#' x <- powerset(5)
#' ans0 <- count_stats(x[1:20], c("mutual", "edges"))
#' 
#' # Calculating using summary_formula
#' fm <- x[[i]] ~ mutual + edges
#' ans1 <- lapply(1:20, function(i) {
#'   environment(fm) <- environment()
#'   ergm::summary_formula(fm)
#' })
#' 
#' ans1 <- do.call(rbind, ans1)
#' 
#' # Comparing
#' all.equal(unname(ans0), unname(ans1))
#' 
count_stats <- function(X, ...) UseMethod("count_stats")

#' Parses a formula to the model parameters (and attributes if applicable)
#' @param x A formula.
#' @noRd
analyze_formula <- function(x, check_w_ergm = TRUE) {
  
  # Getting the parameters
  terms_passed <- attr(stats::terms(x), "term.labels")
  
  # Do these have attributes
  has_attr <- grepl("^[a-zA-Z0-9]+[(].+[)]", terms_passed)
  
  # Capturing attributes
  terms_attrs <- vector("list", length(terms_passed))
  for (term in which(has_attr)) {
    terms_attrs[[term]] <- gsub("\"?[)].*", "", gsub(".+[(]\"?", "", terms_passed[term]))
  }
  
  terms_names <- gsub("[(].+", "", terms_passed)
  
  # Do the terms exists?
  terms_exists <- sapply(terms_names, function(z) {
    out <- capture.output(ergm::search.ergmTerms(name = z))[1]
    !grepl("^No terms named", out, ignore.case = TRUE)
    })
  
  if (any(!terms_exists))
    stop("The following terms are not found in `ergm`: ", 
         paste(terms_names[!terms_exists], sep = ", "), ".", call. = FALSE)
  
  list(
    passed    = terms_passed,
    names     = terms_names,
    attrnames = terms_attrs,
    attrs     = vector("list", length(terms_names))
  )
  
  
}

#' @export
#' @rdname count_stats
count_stats.formula <- function(X, ...) {
  
  # Retrieving networks
  LHS <- eval(X[[2]], envir = environment(X))
  if (network::is.network(LHS))
    LHS <- list(LHS)
  
  # Analyzing the formula
  ergm_terms <- analyze_formula(X)
  
  # # Can we do it?
  # available <- which(!(ergm_terms$names %in% count_available()))
  # if (length(available)) 
  #   stop("The following terms cannot be computed")
  
  # Capturing attributes
  for (a in seq_along(ergm_terms$attrnames)) {
    ergm_terms$attrs[[a]] <- if (is.null(ergm_terms$attrnames[[a]]))
      numeric(0)
    else
      lapply(LHS, function(net) {
        network::get.vertex.attribute(net, attrname = ergm_terms$attrnames[[a]])
      })
  }
  
  # Coercing into the appropiate type
  if (network::is.network(LHS))
    LHS <- list(network::as.matrix(LHS))
  else if (is.list(LHS)) {
    is_net <- sapply(LHS, network::is.network)
    LHS[is_net] <- lapply(LHS[is_net], as.matrix)
  }
  
  out <- matrix(nrow = nnets(LHS), ncol = length(ergm_terms$names),
                dimnames = list(NULL, ergm_terms$passed))
  for (j in 1:ncol(out)) {
    
    out[, j] <- count_stats(
      X = LHS,
      ergm_terms$names[j],
      ergm_terms$attrs[[j]]
      )
    
  }
  
  out

  
}

#' @export
#' @rdname count_stats
count_stats.list <- function(X, terms, attrs = NULL) {
  
  chunks <- make_chunks(length(X), 2e5)
  
  if (!length(attrs))
    attrs <- replicate(length(X), numeric(0), simplify = FALSE)
  
  ans <- matrix(NA, nrow = length(X), ncol=length(terms))
  
  for (s in seq_along(chunks$from)) {
    
    i <- chunks$from[s]
    j <- chunks$to[s]
    
    ans[i:j,] <- count_stats.(X[i:j], terms, attrs)
    
  }
  
  ans
  
}
