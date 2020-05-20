#' @param x List of networks
#' @param x_ref A single object of class `network` that provides the reference
#' values.
#' @param attr_type Type of attribute to replicate
#' @param attr_name Character vector of the attributes to replicate
#' @noRd
replicate_attr <- function(x, x_ref, attr_type, attr_name) {

  # Basic check
  if (!(attr_type %in% c("vertex", "edge", "network")))
    stop(
      "The -type- of network attribute \"", attr_type, "\" does not exist. ",
      "The only valid types are: \"vertex\", \"edge\", and \"network\".",
      call. = FALSE
      )
    
  # Getting the functions
  getter <- getFromNamespace(
    sprintf("get.%s.attribute", attr_type),
    ns = getNamespace("network")
    )
  
  setter <- getFromNamespace(
    sprintf("set.%s.attribute", attr_type),
    ns = getNamespace("network")
  )
  
  # Looping through the values
  for (a in seq_along(attr_name)) {
    
    # Getting the reference value
    value. <- getter(x_ref, attrname = attr_name[[a]])
    
    # Looping through the networks in x
    for (i in seq_along(x)) 
      setter(
        x        = x[[i]],
        attrname = attr_name[[a]],
        value    = value.
        )
  }
    
  x
}


#' Generate an ergmito model
#' @param object Either a model or an object of class [ergmito].
#' @param theta passed
#' @param ... Passed to [ergmito_formulae()].
#' @export
new_rergmito2 <- function(model, theta=NULL, ...) UseMethod("new_rergmito2")

#' @export
new_rergmito2.formula <- function(
  model, 
  theta = NULL,
  ...,
  env = parent.frame()
) {
  
  # Part 1: Get the list of networks -------------------------------------------
  LHS <- eval(model[[2]], envir = env)
  
  # If we have more than one graph, se recycle the formulae object (so they have
  # a shared model) but each object will have its own sampler.
  if (nnets(LHS) > 1L) {
    
    model.              <- stats::update.formula(model, LHS[[i]] ~ .)
    environment(model.) <- environment()
    ans <- vector("list", nnets(LHS))
    for (i in seq_along(ans))
      ans[[i]] <- do.call(
        new_rergmito2,
        c(list(model = model., env = environment()), list(...))
        )
    
    return(ans)
    
  }
  
  # Part 2: Generate the model -------------------------------------------------
  formulae <- ergmito_formulae(model = model, env = env, ...)
  
  # Part 3: Generating the powerset and adjusting for attributes ---------------
  pset <- powerset(nvertex(LHS), directed = is_directed(LHS))
  
  # Part 4: Counting the attributes on the powerset ----------------------------
  
  if (!is_directed(LHS) | !all(formulae$term_fun %in% AVAILABLE_STATS()) | (Sys.getenv("ERGMITO_TEST") != "")) {
    
    # We will need to check whether the model has attributes or not. For now, if
    # the user specified a model update, we assume that the model includes
    # attributes since these could be graph level attributes that are harder to
    # capture.
    if (nrow(formulae$vertex_attrs) | ( ("model_update" %in% names(list(...))) )) {
      
      # Turning into network class objects
      pset <- matrix_to_network(pset, directed = is_directed(LHS))
      
      # Passing the vertex and edge attributes
      pset <- replicate_attr(
        pset,
        x_ref     = LHS,
        attr_type = "vertex",
        attr_name = subset(formulae$vertex_attrs, type == "vertex")$attr
      )
      
      pset <- replicate_attr(
        pset,
        x_ref     = LHS,
        attr_type = "edge",
        attr_name = subset(formulae$vertex_attrs, type == "edge")$attr
      )
      
      # Network attributes (just in case)
      nattrs <- network::list.network.attributes(LHS)
      
      pset <- replicate_attr(
        pset,
        x_ref     = LHS,
        attr_type = "network",
        attr_name = nattrs
      )
      
    }
    
    # We now calculate all the statisics
    model_str <- deparse(stats::update.formula(model, p. ~ .))
    counts <- sapply(pset, function(p.) {
      ergm::summary_formula(as.formula(model_str))
    })
    
  } else {
    
    # Preparing counters
    attrs2pass <- NULL
    if (nrow(formulae$used_attrs)) {
      
      # Listing what are the attributes to be passed
      attrs2pass <- lapply( 
        X = formulae$term_attrs,
        function(a.) {
          
          if (!length(a.))
            return(double(0L))
          
          a. <- a.[names(a.) == "vertex"]
          
          if (!length(a.))
            return(double(0L))
          else if (length(a.) > 1L)
            stop(
              "For now, terms with more than one attribute are not supported on. ",
              "The current model you are trying to fit uses the following attributes ",
              "for a single term: ", a., call. = FALSE
            )
          
          network::get.vertex.attribute(LHS, a.[names(a.) == "vertex"])
      })
    }
    
    counts <- matrix(0, ncol = length(formulae$term_fun), nrow = length(pset)) 
    for (k in seq_along(formulae$term_fun))
      counts[, k] <- count_stats(
        X     = pset,
        terms = formulae$term_fun[k],
        attrs = list(attrs2pass[[k]])
      )
    
    
  }
  
  # Part 3: Functions to sample, re-compute, and subset ------------------------
  call_env <- environment()
  prob     <- NULL
  
  ans_calc <- function(theta. = NULL) {
    
    call_env$prob <- exp(exact_loglik(
      x             = call_env$counts,
      params        = if (is.null(theta.)) theta else theta.,
      stats_weights = call_env$formulae$stats_weights,
      stats_statmat = call_env$formulae$stats_statmat,
      target_offset = call_env$formulae$target_offset,
      stats_offset  = call_env$formulae$stats_offset 
      ))
    
    invisible()
    
  }
  
  call_env$ans_calc()
  return(
    list(counts = counts, pset = pset, calc = ans_calc, prob=prob)
  )
  
}