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
  getter <- utils::getFromNamespace(
    sprintf("get.%s.attribute", attr_type),
    ns = getNamespace("network")
    )
  
  setter <- utils::getFromNamespace(
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

# This avoids a warning
if (getRversion() >= "2.15.1")
  utils::globalVariables("type")

#' ERGMito sampler
#' 
#' Create a sampler object that allows you simulating streams of small networks
#' fast.
#' 
#' @param model A formula.
#' @param theta Named vector. Model parameters.
#' @param ... Further arguments passed to [ergmito_formulae()].
#' 
#' @details 
#' While the \CRANpkg{ergm} package is very efficient, it was not built to do some
#' of the computations required in the ergmito package. This translates in having
#' some of the functions of the package (ergm) with poor speed performance. This
#' led us to "reinvent the wheel" in some cases to speed things up, this includes
#' calculating observed statistics in a list of networks.
#' 
#' @return An environment with the following objects:
#' 
#' - `calc_prob` A function to calculate each graph's probability under the
#'   specified model.
#' - `call` A language object with the call.
#' - `counts` A list with 3 elements: `stats` the sufficient statistics of each network,
#'   `weights` and `statmat` the overall matrices of sufficient statistics used to
#'   compute the likelihood.
#' - `network` The baseline network used to either fit the model or obtain
#'   attributes.
#' - `networks` A list with the actual sample space of networks.
#' - `probabilities` A numeric vector with each graph's probability.
#' - `sample` A function to draw samples. `n` specifies the number of samples to
#'   draw and `theta` the parameter to use to calculate the likelihoods.
#' - `theta` Named numeric vector with the current values of the model parameters.
#' 
#' @export
#' @examples 
#' 
#' # We can generate a sampler from a random graph
#' set.seed(7131)
#' ans <- new_rergmito(rbernoulli(4) ~ edges, theta = -.5)
#' 
#' # Five samples
#' ans$sample(5)
#' 
#' # or we can use some nodal data:
#' data(fivenets)
#' ans <- new_rergmito(
#'   fivenets[[3]] ~ edges + nodematch("female"),
#'   theta = c(-1, 1)
#' )
#' 
#' # Five samples
#' ans$sample(5) # All these networks have a "female" vertex attr
#' 
new_rergmito <- function(model, theta, ...) UseMethod("new_rergmito")

#' @export
new_rergmito.formula <- function(
  model, 
  theta,
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
        new_rergmito,
        c(list(model = model., theta=theta, env = environment()), list(...))
        )
    
    return(ans)
    
  } else if (!inherits(LHS, "network") && is.list(LHS))
    LHS <- LHS[[1L]]
  
  # Initializing the object
  sampler <- new.env()
  
  # Part 2: Generate the model -------------------------------------------------
  sampler$formulae <- ergmito_formulae(model = model, env = env, ...)
  
  # This is a basic check
  if (length(theta) != sampler$formulae$npars)
    stop(
      sprintf(
        paste(
          "The length of -theta- (%i) does not match the",
          "number of parameters in the model (%i)."
        ),
        length(theta),
        sampler$formulae$npars
      ),
      call. = FALSE
      )
  
  # Part 3: Generating the powerset and adjusting for attributes ---------------
  sampler$networks <- powerset(nvertex(LHS), directed = is_directed(LHS))
  
  # Part 4: Counting the attributes on the powerset ----------------------------
  
  sampler_w_attributes <- FALSE
  if (!is_directed(LHS) | !all(sampler$formulae$term_fun %in% AVAILABLE_STATS()) | (Sys.getenv("ERGMITO_TEST") != "")) {
    
    # Turning into network class objects
    sampler$networks <- matrix_to_network(sampler$networks, directed = is_directed(LHS))
    
    # We will need to check whether the model has attributes or not. For now, if
    # the user specified a model update, we assume that the model includes
    # attributes since these could be graph level attributes that are harder to
    # capture.
    if (nrow(sampler$formulae$used_attrs) | ( ("model_update" %in% names(list(...))) )) {
      
      # Passing the vertex and edge attributes
      sampler$networks <- replicate_attr(
        sampler$networks,
        x_ref     = LHS,
        attr_type = "vertex",
        attr_name = subset(sampler$formulae$used_attrs, type == "vertex")$attr
      )
      
      sampler$networks <- replicate_attr(
        sampler$networks,
        x_ref     = LHS,
        attr_type = "edge",
        attr_name = subset(sampler$formulae$used_attrs, type == "edge")$attr
      )
      
      # Network attributes (just in case)
      nattrs <- network::list.network.attributes(LHS)
      
      sampler$networks <- replicate_attr(
        sampler$networks,
        x_ref     = LHS,
        attr_type = "network",
        attr_name = nattrs
      )
      
    } 
    
    # We now calculate all the statisics
    model_str <- deparse(stats::update.formula(model, p. ~ .))
    sampler$counts <- t(sapply(sampler$networks, function(p.) {
      ergm::summary_formula(as.formula(model_str))
    }))
    
  } else {
    
    # Preparing counters
    attrs2pass <- NULL
    if (nrow(sampler$formulae$used_attrs)) {
      
      # Listing what are the attributes to be passed
      attrs2pass <- lapply( 
        X = sampler$formulae$term_attrs,
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
      
      sampler_w_attributes <- TRUE
      
    } else {
      attrs2pass <- replicate(
        n        = length(sampler$formulae$term_attrs),
        expr     = double(0L),
        simplify = FALSE
        )
    }
    
    sampler$counts <- matrix(
      0, ncol = length(sampler$formulae$term_fun), nrow = length(sampler$networks)
      ) 
    for (k in seq_along(sampler$formulae$term_fun))
      sampler$counts[, k] <- count_stats(
        X     = sampler$networks,
        terms = sampler$formulae$term_fun[k],
        attrs = list(attrs2pass[[k]])
      )
    
  }
  
  # Making sure we have the names
  model. <- stats::update.formula(model, LHS ~ .)
  environment(model.) <- environment()
  colnames(sampler$counts) <- names(ergm::summary_formula(model.))
  
  # Figuring out
  # Are we updating the model? ------------------------------------------------
  graph_attrs <- graph_attributes_as_df(LHS)
  model_frame <- model_frame_ergmito(
    formula        = model,
    formula_update = sampler$formulae$model_update, 
    data           = sampler$counts,
    g_attrs.       = graph_attrs
  )
  
  sampler$counts        <- model_frame$stats
  sampler$target_offset <- model_frame$offsets
  
  # Part 5: Functions to sample, re-compute, and subset ------------------------
  sampler$probabilities <- NULL
  
  sampler$calc_probabilities <- function(theta. = NULL) {
    
    sampler$probabilities <- exact_loglik(
      x             = sampler$counts,
      params        = if (is.null(theta.)) theta else theta.,
      stats_weights = sampler$formulae$stats_weights,
      stats_statmat = sampler$formulae$stats_statmat,
      target_offset = sampler$target_offset,
      stats_offset  = sampler$formulae$stats_offset,
      as_prob       = TRUE
      )
    
    # Safeguard
    test_sum <- sum(sampler$probabilities)
    if (abs(test_sum - 1.0) > 1e-8)
      stop(
        "Something went wrong. The sum of the probabilities is different than 1. ",
        "This is an unexpected error. Check the passed parameters to new_rergmito(). ",
        "You can report bugs at https://github.com/muriteams/ergmito/issues",
        call. = FALSE
      )
    
    # Removing machine precision error.
    sampler$probabilities <- sampler$probabilities/test_sum
    
    invisible()
    
  }
  
  sampler$calc_probabilities()
  
  # A getter function ----------------------------------------------------------
  sampler$get_networks <- function(idx = NULL) {
    
    if (is.null(idx))
      idx <- seq_along(sampler$networks)
    
    nets <- sampler$networks[idx]
    
    # This is only needed if we are not using statnet's summary function.
    if (sampler_w_attributes) {
      
      # Turning into network class objects
      nets <- matrix_to_network(nets, directed = is_directed(LHS))
      
      # Passing the vertex and edge attributes
      nets <- replicate_attr(
        nets,
        x_ref     = LHS,
        attr_type = "vertex",
        attr_name = subset(sampler$formulae$used_attrs, type == "vertex")$attr
      )
      
      nets <- replicate_attr(
        nets,
        x_ref     = LHS,
        attr_type = "edge",
        attr_name = subset(sampler$formulae$used_attrs, type == "edge")$attr
      )
      
      # Network attributes (just in case)
      nattrs <- network::list.network.attributes(LHS)
      
      replicate_attr(
        nets,
        x_ref     = LHS,
        attr_type = "network",
        attr_name = nattrs
      )
      
    } else
      nets
    
  }
  
  # Sampling functions ---------------------------------------------------------
  sampler$sample <- function(n, theta = NULL, as_indexes = FALSE) {
    
    # If no new set of parameters is used, then 
    if (length(theta)) {
      oldp <- sampler$probabilities
      sampler$calc_probabilities(theta)
      on.exit(sampler$probabilities <- oldp)
    } 
    
    idx <- sample.int(
      n       = length(sampler$probabilities),
      size    = n,
      replace = TRUE,
      prob    = sampler$probabilities,
      useHash = FALSE
    )
    
    if (!as_indexes) 
      sampler$get_networks(idx)
    else 
      idx
    
  }
  
  # Call
  sampler$call     <- match.call()
  sampler$network  <- LHS
  sampler$theta    <- theta

  structure(
    sampler,
    class = "ergmito_sampler"
  )
  
}

#' @export
#' @param x An object of class `ergmito_sampler`.
#' @rdname new_rergmito
#' @param i `i` is an integer vector indicating the indexes of the networks to
#' draw.
#' @details The indexing method, `[.ergmito_sampler`, allows extracting networks
#' directly by passing indexes. `i` indicates the index of the networks to draw,
#' which go from 1 through `2^(n*(n-1))` if directed and `2^(n*(n-1)/2)` if
#' undirected . 
#' @return The indexing method `[.ergmito_sampler` returns a list of networks
`[.ergmito_sampler` <- function(x, i, ...) {
  
  # Sampling networks
  if (missing(i))
    i <- seq_along(x$networks)
  
  x$get_networks(i)
  

}

#' @export
print.ergmito_sampler <- function(x, ...) {
  
  cat("ERGMito simulator\n")
  cat("Call   :", deparse(x$call), "\n")
  cat("sample :", deparse(x$sample)[1], "\n")
  
  invisible(x)
  
}

#' @export
#' @importFrom utils ls.str
str.ergmito_sampler <- function(object, ...) {
  
  utils::ls.str(object, ...)
  
}
