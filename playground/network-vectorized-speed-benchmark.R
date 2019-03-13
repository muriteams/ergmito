library(ergmito)
library(network)

set.seed(12)
A <- rbernoulli(50)
net0 <- network(A)

net1 <- ergmito:::matrix_to_network(list(A))

library(bench)

ans <- bench::mark(
  ergm = network::network(A),
  ergmito = ergmito::matrix_to_network(list(A)),
  check = FALSE,
  relative = TRUE,
  iterations = 100
  )

plot(ans)
ans

neta <- network::network(A)

# ans <- bench::mark(
#   ergm = w0 <- network::set.vertex.attribute(neta, "a", 1:50),
#   ergmito = w1 <- ergmito::add_vertex_attr(neta, list(1:50), "b"),
#   check = FALSE,
#   relative = TRUE
# )
# 
# plot(ans)
# ans

networks_list <- function(x) {
  res <- vector("list", length(x))
  for (i in seq_along(res))
    res[[i]] <- network(x[[i]])
  res
}

adjmats <- rbernoulli(rep(5, 2000))
nets <- matrix_to_network(adjmats)

(ans <- bench::mark(
  ergmito = matrix_to_network(adjmats),
  ergm    = networks_list(adjmats),
  check=FALSE
))
plot(ans)

# add_vattr_network <- function(x, attrvalue, attrname) {
#   
#   res <- vector("list", length(x))
#   for (i in seq_along(res))
#     res[[i]] <- set.vertex.attribute(x[[i]], attrname = attrname, value=attrvalue)
#   res
#   
# }
# 
# (ans <- bench::mark(
#   ergmito = ergmito::add_vertex_attr(nets, list(1:5), "a"),
#   ergm    = add_vattr_network(nets, list(1:5), "a"),
#   check=FALSE, iterations = 100
# ))
# plot(ans)
