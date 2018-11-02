library(readr)
library(magrittr)
library(similR)

#read in data of participant perceptions of advice ties among their other teammates 
advice_css <- readRDS("data/networks_advice_css.rds")
truth      <- readRDS("data/networks_truth.rds")[["3"]]$advice 

groups_ids <- names(advice_css) %>%
  unique

# Imputing missing preception
LAS        <- vector("list", length(advice_css))
names(LAS) <- groups_ids

for (g in groups_ids) {
  for (i in names(advice_css[[g]])) 
    advice_css[[g]][[i]][i,] <- truth[[g]][i,]
  
  # Creating the LAS matrix
  LAS[[g]] <- las(advice_css[[g]], rule = "i")
  
  dimnames(LAS[[g]]) <- dimnames(advice_css[[g]][[1]])
}
    
# Problem with the CSS network 31 B
# advice_css$`31B`
#     31A 31B 31C 31D
# 31A  NA   1   1   1
# 31B 999  NA 999 999
# 31C   0   0  NA   0
# 31D  99  99  99  NA

saveRDS(LAS, "data/networks_advice_las.rds")
