library(ergmito)

data(fivenets)

# This is how you fit the model
res <- ergmito(
  fivenets ~ edges + nodeocov("female"),
  model_update = ~ . -nodeocov.female + offset(I(nodeocov.female * -1e5))
  )

summary(res)

# We can use the fitted value to create a sampler
sim <- new_rergmito(
  fivenets ~ edges + nodeocov("female"),
  theta = coef(res),
  model_update = ~ . -nodeocov.female + offset(I(nodeocov.female * -1e5))
)

# And with it, see if we get any network where females have outdegree different than 0
count_stats(sim[[1]]$sample(100) ~ edges + nodeicov("female") + nodeocov("female")) |>
  summary()

count_stats(sim[[2]]$sample(100) ~ edges + nodeicov("female") + nodeocov("female")) |>
  summary()


