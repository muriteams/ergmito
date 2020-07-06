#' Offset processer
#' @param x An object of class formula
#' @details
#' It will identify if the model has offset terms or not. In the case of
#' having offsets, the resulting object has an attribute indicating which
#' of the terms are offset terms.
#' @return A list of formulas.
#' @noRd
parse_offset <- function(x) {
  
  terms_x <- terms(x)
  
  if (!length(attr(terms_x, "term.labels")))
    stop(
      "Invalid model:\n",
      deparse(x), 
      "\nA offset-only model cannot be fitted.", call. = FALSE)
  
  # Are there any offset terms?
  offsets_x <- attr(terms_x, "offset")
  if (!length(offsets_x))
    return(list(x))
  
  # Updating offset terms
  variables <- rownames(attr(terms_x, "factors"))
  variables_free <- gsub("^offset[(]|[)]$", "", variables[offsets_x])
  
  # We need to memorize this later
  offsets_x <- structure(offsets_x, names = variables[offsets_x])
  variables[offsets_x] <- variables_free
  
  if (attr(terms_x, "response") > 0) 
    ans <- sprintf("%s ~ %s", variables[1L], variables[-1L])
  else 
    ans <- sprintf(" ~ %s", variables)
  
  ans <- lapply(ans, formula, env = environment(x))
  attr(ans, "ergmito_offset") <- offsets_x
  
  return(ans)
  
}

#' Checks the passed formula (puts it into the right form)
#' @noRd
model_update_parser <- function(x) {
  
  if (is.null(x))
    return(NULL)
  
  # parsing the formula object
  formula_update_d <- deparse(x, width.cutoff = 500L)
  
  if (length(formula_update_d) > 1L)
    formula_update_d <- paste(formula_update_d, collapse = " ")
  
  if (!grepl("^\\s*[~]", formula_update_d))
    stop(
      "The LHS of the update to the model (-model_update-) should be empty. ",
      "Right now, it has the following: ", gsub("[~].+", "", formula_update_d),
      call. = FALSE
    )
  if (!grepl("[~]\\s*[.]\\s*[+]", formula_update_d)) {
    formula_update_d <- gsub("[~]", "~ . +", formula_update_d)
  }
  
  as.formula(formula_update_d, env = environment(x))
}

#' 
#' @noRd
model_frame_ergmito <- function(formula, formula_update, data, g_attrs.) {
  
  # If nothing to do, then return nothing
  if (is.null(formula_update)) 
    return(list(model = formula, stats = data, offsets = double(nrow(data))))
  
  model_final <- stats::as.formula(
    sprintf("~ %s", paste(colnames(data), collapse = " + "))
  )
  
  # Building the new formula
  model_final <- stats::update.formula(model_final, formula_update)
  
  # Parsing offset terms. This removes the offset() around the term
  # names and then it adds an attribute that tags what are the offset
  # variables so it is easier to extract them from the model.
  model_final <- parse_offset(model_final)
  
  # Appending graph stats
  data <- as.data.frame(data)
  if (length(g_attrs.)) {
    data <- cbind(
      data, 
      if ((nrow(g_attrs.) == 1) && (nrow(data) > 1)) {
        do.call(rbind, replicate(
          nrow(data),
          g_attrs.,
          simplify = FALSE
        ))
      } else 
        g_attrs.
    )
  }
  
  # Computing stats
  data <- lapply(model_final, stats::model.frame, data = data)
  
  # Checking offset terms
  offset_terms <- attr(model_final, "ergmito_offset")
  if (length(offset_terms)) {
    
    # Extracting the offsets
    offsets. <- do.call(cbind, data[offset_terms])
    offsets. <- rowSums(as.matrix(offsets.))
    data     <- data[-offset_terms]
    data     <- as.matrix(do.call(cbind, data))
    
  } else {
    
    data     <- as.matrix(do.call(cbind, data))
    offsets. <- double(nrow(data))
    
  }
  
  list(
    model   = stats::update.formula(formula, formula_update),
    stats   = data,
    offsets = offsets.
  )
  
}

#' Parser of formulas, with an emphasis on checking terms against a list of
#' available terms.
#' 
#' @param x A [stats::formula].
#' @param graph A graph of reference, this is to compare names against
#' @noRd
analyze_formula <- function(x, graph = NULL) {
  
  if (!inherits(x, "formula"))
    stop("`x` must be a formula.", call. = FALSE)
  
  trms <- stats::terms(x)
  
  # Capture pattern
  pat    <- "(\"[^\"]+\")|(\'[^\']+\')"
  
  # Need to take care of the type of network that is beeing passed, if it is
  # empty, easy, if it is network, easy as well, but if it is something else,
  # then a bit complex
  obs_v_attrs <- NULL
  obs_e_attrs <- NULL
  obs_n_attrs <- NULL
  if (length(graph)) {
    
    if (!inherits(graph, "network")) {
      if (inherits(graph, "list"))
        if (inherits(graph[[1L]], "network")) {
          obs_v_attrs <- network::list.vertex.attributes(graph[[1L]])
          obs_e_attrs <- network::list.edge.attributes(graph[[1L]])
          obs_n_attrs <- network::list.network.attributes(graph[[1L]])
        }
        
    } else {
      obs_v_attrs <- network::list.vertex.attributes(graph)
      obs_e_attrs <- network::list.edge.attributes(graph)
      obs_n_attrs <- network::list.network.attributes(graph)
    }
  }
  
  # Parsing all labels
  attr_types <- NULL
  anames <- lapply(attr(trms, "term.labels"), function(t.) {
    
    # Capturing varnames
    m <- gregexpr(pat, t., perl = TRUE)
    v <- regmatches(t., m)[[1L]]
    
    # Removing quotes
    v <- gsub("^[\"\']|[\"\']$", "", v)
    
    if (length(obs_e_attrs)) {
      v_e <- intersect(v, obs_e_attrs)
      v_e <- structure(v_e, names = rep("edge", length(v_e)))
    } else
      v_e <- NULL
    
    if (length(obs_v_attrs)) {
      v_v <- intersect(v, obs_v_attrs)
      v_v <- structure(v_v, names = rep("vertex", length(v_v)))
    } else
      v_v <- NULL
    
    if (length(obs_n_attrs)) {
      v_n <- intersect(v, obs_n_attrs)
      v_n <- structure(v_n, names = rep("network", length(v_n)))
    } else
      v_n <- NULL
    
    v <- setdiff(v, union(union(v_v, v_e), v_n))
    v <- structure(v, names = rep("unknown", length(v)))
    
    c(v_v, v_e, v_n, v)
  })
  

  # Nice result
  all_attrs <- unlist(anames, use.names = TRUE)
  all_attrs <- data.frame(
    type = names(all_attrs),
    attr = unname(all_attrs),
    stringsAsFactors = FALSE
  )
  all_attrs <- unique(all_attrs)
  
  list(
    term_passed = attr(trms, "term.labels"),
    term_names  = gsub("[(].+", "", attr(trms, "term.labels")),
    term_attrs  = anames,
    attrs       = replicate(length(anames), double(0L), simplify = FALSE),
    nattrs      = sum(all_attrs$type != "unknown"),
    all_attrs   = all_attrs
  )
  
  
  
}

#' This function extracts networks attributes and returns a data
#' frame
#' @noRd
graph_attributes_as_df <- function(net) {
  
  if (inherits(net, "list") && nnets(net) > 1L) {
    net <- matrix_to_network(net)
    return(do.call(rbind, lapply(net, graph_attributes_as_df)))
  }
  
  if (!network::is.network(net))
    net <- matrix_to_network(net)
  
  if (is.list(net) && !is.network(net))
    net <- net[[1]]
  
  # Listing attributes, extracting and naming
  netattrs   <- network::list.network.attributes(net)
  ans        <- lapply(netattrs, network::get.network.attribute, x = net)
  names(ans) <- netattrs
  
  # Returning as data.frame
  as.data.frame(ans)
}

gmodel <- function(model, net) {
  
  netattrs   <- network::list.network.attributes(net)
  ans        <- lapply(netattrs, network::get.network.attribute, x = net)
  names(ans) <- netattrs
  
  ans <- stats::model.matrix(
    stats::update.formula(model, ~ 0 + .),
    as.data.frame(ans)
  )
  
  structure(
    ans[1, , drop=TRUE],
    names = colnames(ans)
  )
}
