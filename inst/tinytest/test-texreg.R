data(fivenets)
suppressMessages(library(texreg))

ans <- ergmito(fivenets ~ edges + mutual)
expect_output(print(texreg(ans)), "tabular")
