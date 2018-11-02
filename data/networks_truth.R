library(readr)
library(haven)
library(magrittr)
library(tidyr)
library(dplyr)
library(magrittr)

#read in data of participant perceptions of advice ties among their other teammates 

IDS   <- c(LETTERS[1:4], NA, "E")
IDS_E <- c(LETTERS[1:3], NA, "D")

groups_sizes <- readr::read_csv("data-raw/Study1_Group sizes.csv") %>%
  rename(group = Group, n = groupSize)


edgelist_pre <- haven::read_sav("data-raw/MURI_Survey_1.sav") %>%

  # Some weird cases with no upper case
  mutate(Survey1_PID = toupper(Survey1_PID)) %>%
  select(Survey1_PID, matches("^T[A-Z]_[0-9]_[0-9]")) %>%
  
  # Reshaping dataset as an edgelist, and keeping only values != NA
  gather("Code", "Value", -Survey1_PID) %>%
  filter(!is.na(Value)) %>%
  
  # Retrieving ego/alter
  transmute(
    group = as.integer(stringr::str_extract(Survey1_PID, "^[0-9]+")),
    ego   = stringr::str_extract(Survey1_PID, "[A-Z]$"),
    alter = as.integer(stringr::str_extract(Code, "[0-9]$")),
    alter = case_when(
      ego == "A" ~ IDS[-1][alter],
      ego == "B" ~ IDS[-2][alter],
      ego == "C" ~ IDS[-3][alter],
      ego == "D" ~ IDS[-4][alter],
      ego == "E" ~ IDS_E[alter],
      TRUE ~ "ERROR"
    ),
    survey = 1L,
    network = stringr::str_extract(Code, "(?<=[_])[0-9]+")
  ) 

edgelist_post <- haven::read_sav("data-raw/MURI_Survey_3.sav") %>%
  
  # Some weird cases with no upper case
  mutate(Survey3_PID = toupper(Survey3_PID)) %>%
  select(Survey3_PID, matches("^T[A-Z]_[0-9]post_[0-9]")) %>%
  
  # Reshaping dataset as an edgelist, and keeping only values != NA
  gather("Code", "Value", -Survey3_PID) %>%
  filter(!is.na(Value)) %>%
  
  # Retrieving ego/alter
  transmute(
    group = as.integer(stringr::str_extract(Survey3_PID, "^[0-9]+")),
    ego   = stringr::str_extract(Survey3_PID, "[A-Z]$"),
    alter = as.integer(stringr::str_extract(Code, "[0-9]$")),
    alter = case_when(
      ego == "A" ~ IDS[-1][alter],
      ego == "B" ~ IDS[-2][alter],
      ego == "C" ~ IDS[-3][alter],
      ego == "D" ~ IDS[-4][alter],
      ego == "E" ~ IDS_E[alter],
      TRUE ~ "ERROR"
    ),
    survey = 3L,
    network = stringr::str_extract(Code, "(?<=[_])[0-9]+")
  ) 

edgelists <- bind_rows(edgelist_post, edgelist_pre) %>%
  # Merging group size. We use full to make sure that we have all the data
  full_join(groups_sizes)

# saveRDS(edgelist, "data/edgelist_truth.rds")

# Now, creating the adjmat
networks <- vector("list", 2)
names(networks) <- c(1, 3)

for (s in c(1, 3)) {
  
  networks[[as.character(s)]] <- edgelists %>%
    
    filter(survey == s) %>%
    
    split(.$group) %>%
    lapply(function(g) {
      
      # Group ids
      n   <- g$n[1]
      ids <- LETTERS[1:n]
      
      # Looping at the group level
      M <- matrix(0, ncol = n, nrow = n, dimnames = list(ids, ids))
      M <- replicate(n, M, simplify = FALSE)
      names(M) <- ids
      
      # Keep complete cases
      g <- g[complete.cases(g),]
      
      if (nrow(g) == 0)
        return(M)
      
      message("Group ", g$group[1], "survey ", s, " done.")
      
      # Adding 1s
      for (i in 1:nrow(g))
        M[[g$viewer[i]]][g$ego[i], g$alter[i]] <- 1L
      
      M
      
    })
  
}


# saveRDS(networks, "data/networks_truth.rds")