data(fivenets)
bm <- benchmarkito(
  ergmito = count_stats(fivenets ~ edges + mutual),
  ergm    = lapply(fivenets, function(n) {
    ergm::summary_formula(n ~ edges + mutual)
  }), 
  times = 50
  )

expect_output(print(bm), "elapsed")
expect_silent(plot(bm))
