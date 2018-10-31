library(readr)
library(haven)
library(magrittr)
library(tidyr)
library(dplyr)
library(magrittr)

#read in data of participant perceptions of advice ties among their other teammates 

IDS     <- c(LETTERS[1:4], NA, "E")

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
      ego == "E" ~ IDS[-5][alter],
      TRUE ~ "ERROR"
    ),
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
      ego == "E" ~ IDS[-5][alter],
      TRUE ~ "ERROR"
    ),
    network = stringr::str_extract(Code, "(?<=[_])[0-9]+")
  ) 

saveRDS(edgelist_post, "data/edgelist_post.rds")
saveRDS(edgelist_pre, "data/edgelist_pre.rds")
