#' Count Network Statistics
#' 
#' This function is similar to what [ergm::summary_formula] does, but it provides
#' a fast wrapper suited for matrix class objects.
#' @param X List of square matrices. (networks)
#' @param terms Character vector with the names of the statistics to calculate.
#' Currently, the only available statistics are: '\Sexpr{paste(ergmito::AVAILABLE_STATS(), collapse="', '")}'.
#' @param ... Passed to the method.
#' @param attrs A list of vectors. This is used when `term` has a nodal attribute
#' such as `nodeicov(attrname="")`.
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

#' @export
#' @rdname count_stats
AVAILABLE_STATS <- function() count_available()

#' Parses a formula to the model parameters (and attributes if applicable)
#' @param x A formula.
#' @noRd
analyze_formula <- function(x, check_w_ergm = FALSE) {
  
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
  if (check_w_ergm) {
    terms_exists <- sapply(terms_names, function(z) {
      out <- utils::capture.output(ergm::search.ergmTerms(name = z))[1]
      !grepl("^No terms named", out, ignore.case = TRUE)
      })
    
    if (any(!terms_exists))
      stop("The following terms are not found in `ergm`: ", 
           paste(terms_names[!terms_exists], sep = ", "), ".", call. = FALSE)
  }
  
  list(
    passed    = terms_passed,
    names     = terms_names,
    attrnames = terms_attrs,
    attrs     = replicate(length(terms_names), double(0))
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
  ergm_model <- analyze_formula(X)
  
  # # Can we do it?
  # available <- which(!(ergm_model$names %in% count_available()))
  # if (length(available)) 
  #   stop("The following terms cannot be computed")
  
  # Capturing attributes
  for (a in seq_along(ergm_model$attrnames)) {
    ergm_model$attrs[[a]] <- if (is.null(ergm_model$attrnames[[a]]))
      numeric(0)
    else
      lapply(LHS, function(net) {
        network::get.vertex.attribute(net, attrname = ergm_model$attrnames[[a]])
      })
  }
  
  # Coercing into the appropiate type
  if (network::is.network(LHS))
    LHS <- list(as.adjmat(LHS))
  else if (is.list(LHS)) {
    
    is_net <- sapply(LHS, network::is.network)
    
    # Coercing into a net
    for (i in which(is_net))
      LHS[[i]] <- as.adjmat(LHS[[i]])
    
  }
  
  out <- matrix(nrow = nnets(LHS), ncol = length(ergm_model$names),
                dimnames = list(NULL, ergm_model$passed))
  
  for (j in 1:ncol(out)) {
    
    out[, j] <- count_stats(
      X = LHS,
      ergm_model$names[j],
      ergm_model$attrs[[j]]
      )
    
  }
  
  out

  
}

#' @export
#' @rdname count_stats
count_stats.list <- function(X, terms, attrs = NULL, ...) {
  
  chunks <- make_chunks(length(X), 2e5)
  
  # if (!length(attrs))
  #   attrs <- replicate(length(X), numeric(0), simplify = FALSE)
  
  ans <- matrix(NA, nrow = length(X), ncol=length(terms))
  
  for (s in seq_along(chunks$from)) {
    
    i <- chunks$from[s]
    j <- chunks$to[s]
    
    for (k in seq_along(terms)) {
      if (!length(attrs))
        ans[i:j, k] <- count_stats.(X[i:j], terms[k], list(numeric(0)))
      else
        ans[i:j, k] <- count_stats.(X[i:j], terms[k], attrs[i:j])
    }
    
  }
  
  ans
  
}
