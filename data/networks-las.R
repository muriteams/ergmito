library(readr)
library(magrittr)
library(similR)

#read in data of participant perceptions of advice ties among their other teammates 
advice_css <- readRDS("data/networks-advice-css.rds")
truth      <- readRDS("data/networks-truth.rds") %>%
  lapply("[[", "s3q1advice")

groups_ids <- sapply(names(advice_css), gsub, pattern="[A-Z]+$", replacement="") %>%
  unique

MAT <- lapply(groups_ids, function(id) {
  m <- advice_css[grepl(sprintf("^%s", id), names(advice_css))]
  m <- lapply(m, function(z) {
    dimnames(z) <- list(
      gsub("^[0-9]+", "", rownames(z)),
      gsub("^[0-9]+", "", colnames(z))
    )
    z
  
  })
  names(m) <- gsub("^[0-9]+", "", names(m))
  
  m
})

names(MAT) <- groups_ids %>% as.integer %>% as.character
groups_ids <- groups_ids %>% as.integer %>% as.character

# Imputing missing preception
LAS        <- vector("list", length(MAT))
names(LAS) <- groups_ids

for (g in groups_ids) {
  for (i in names(MAT[[g]])) 
    MAT[[g]][[i]][i,] <- truth[[g]][i,]
  
  # Creating the LAS matrix
  LAS[[g]] <- las(MAT[[g]], rule = "i")
}
    
# Problem with the CSS network 31 B
# advice_css$`31B`
#     31A 31B 31C 31D
# 31A  NA   1   1   1
# 31B 999  NA 999 999
# 31C   0   0  NA   0
# 31D  99  99  99  NA

saveRDS(MAT, "data/networks-advice-las.rds")
