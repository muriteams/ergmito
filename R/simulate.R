#' Draw samples from a fitted `ergmito` model
#' @param object An object of class [ergmito].
#' @param nsim Integer scalar. Number of samples to draw from the selected set 
#' of networks.
#' @param seed See [stats::simulate]
#' @param which_networks Integer vector. Specifies what networks to sample from.
#' It must be within 1 and `nnets(object)`.
#' @param theta,... Further arguments passed to [new_rergmito].
#' @export
#' @examples 
#' data(fivenets)
#' fit <- ergmito(fivenets ~ edges + nodematch("female"))
#' 
#' # Drawing 200 samples from networks 1 and 3 from the model
#' ans <- simulate(fit, nsim = 200, which_networks = c(1, 3))
#' @importFrom stats simulate
simulate.ergmito <- function(
  object,
  nsim           = 1,
  seed           = NULL,
  which_networks = 1L,
  theta          = NULL,
  ...
  ) {
  
  # Catching the model and preparing for local evaluation
  model  <- stats::formula(object)
  model. <- stats::update.formula(model, networks[[i]] ~ .)
  environment(model.) <- environment()
  
  if (!is.null(seed))
    set.seed(seed)
  
  # Retrieving the networks
  if (any((which_networks > nnets(object)) | (which_networks < 1)))
    stop("`which_networks` should be an integer within [1, nnets(object)].",
         call. = FALSE)
  
  networks <- eval(model[[2]], environment(model))[which_networks]

  # Step 2: Loop through the networks to generate the predictions:
  samples <- NULL # vector("list", length(which_networks) * nsim)
  dots <- list(...)
  
  if ("sizes" %in% names(dots))
    stop(
      "The parameter `sizes` cannot be passed to simulate. To specify different ",
      "sizes use the function `new_rergmito` directly. Otherwise, use the ",
      "argument `which_networks` to use one or more networks as reference.",
      call. = FALSE
      )
  
  for (i in seq_along(networks)) {
    
    # Generating sample, and later on, adding up matrices
    sampler <- new_rergmito(
      model = model.,
      theta = if (is.null(theta)) stats::coef(object) else theta,
      sizes = nvertex(networks[[i]]),
      ...
    )
    
    samples <- c(samples, sampler$sample(n = nsim, s = nvertex(networks[[i]])))
    
  }
  
  samples
  
}
