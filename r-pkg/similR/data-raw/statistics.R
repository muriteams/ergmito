
library(stringr)
library(magrittr)

getmetric <- readLines("src/similarity.cpp")
getmetric <- getmetric[
  which(grepl("void getmetric", getmetric)):
    (which(grepl("NumericMatrix similarity", getmetric)) - 2)
  ]

getmetric <- getmetric %>%
  str_extract_all("(?<=s\\s[=]{2}\\s.)[a-zA-Z0-9]+") %>%
  extract(sapply(., length)>0) %>%
  sapply("[", 1)

statistics <- list(
  similarity = getmetric[grepl("^s", getmetric)],
  distance   = getmetric[grepl("^d", getmetric)]
)

usethis::use_data(statistics)