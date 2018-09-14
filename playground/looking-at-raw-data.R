
library(tibble)
library(magrittr)

subjects <- read.table(
  "data-raw/MURI_Individual data set 090618.csv",
  sep="\t",
  header=TRUE
  ) %>% as_tibble

# Should this be different from the previous one?
surveys <- foreign::read.spss(
  "data-raw/MURI_AllSurveys - FINAL_050817_w clusters.sav"
  ) %>% as_tibble

