#' Example of a group of small networks
#' 
#' This list of networks was generated using the [new_rergmito] sampler from a
#' set of 5 baseline networks with a random vector of `female`. The sufficient 
#' statistics that generate this data are `edges` and  `nodematch("female")` with parameters
#' -2.0 and 2.0 respectively.
#' 
#' @docType data
"fivenets"

# Five ERGMito samplers
# 
# This list contains five ERGMito samplers. Each one of these was built using
# a random Bernoulli graph with an attribute `female`. The parameters used for
# creating the sampler are `(edges = -2.0, nodematch("female") = 2.0)`. A example
# of a dataset generated with this is [fivenets].
# 
# 
# @docType data
# "fivesamplers"
