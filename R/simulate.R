#' Simulation method
#' @noRd
# simulate.ergmito <- function(
#   object,
#   nsims    = 1,
#   seed     = NULL,
#   mc_cores = 2L,
#   ...) {
#   
#   # Generating samples
#   n    <- network::network.size(object$network)
#   nets <- powerset(n)
#   
#   # Computing statistics
#   fm <- object$formula
#   fm[[2]] <- bquote(net0)
#   stats <- parallel::mclapply(nets, function(net0) {
#     environment(fm) <- environment()
#     summary(fm)
#   })
#   
#   stats0 <- summary(object$formula)
#   
#   stats <- do.call(rbind, stats) -
#     matrix(stats0, nrow = length(stats), ncol=length(stats0), byrow = TRUE)
#   
#   list(
#     net0 = object$network,
#     nets = nets,
#     prob = exp(stats %*% object$coef)*exp(object$loglikelihood)
#   )
#   
# }