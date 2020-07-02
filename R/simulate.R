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
  model. <- stats::update.formula(model, networks ~ .)
  environment(model.) <- environment()
  
  if (!is.null(seed))
    set.seed(seed)
  
  # Retrieving the networks
  if (any((which_networks > nnets(object)) | (which_networks < 1)))
    stop("`which_networks` should be an integer within [1, nnets(object)].",
         call. = FALSE)
  
  networks <- eval(model[[2]], environment(model))[which_networks]

  # Step 2: Loop through the networks to generate the predictions:

  # Generating sample, and later on, adding up matrices
  sampler <- new_rergmito(
    model = model.,
    theta = if (is.null(theta)) stats::coef(object) else theta,
    ...
  )
  
  samples <- replicate(nsim, vector("list", length(which_networks)), simplify = FALSE)
  for (i in 1:length(which_networks)) {
    
    if (length(which_networks) == 1L)
      samples_tmp <- sampler$sample(n = nsim)
    else
      samples_tmp <- sampler[[i]]$sample(n = nsim)
    
    for (j in 1:nsim)
      samples[[j]][i] <- samples_tmp[j]
    
  }
  
  samples
  
}
