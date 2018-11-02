library(readr)
library(haven)
library(magrittr)
library(tidyr)
library(dplyr)
library(magrittr)

#read in data of participant perceptions of advice ties among their other teammates 

IDS   <- c(LETTERS[1:4], NA, "E")
IDS_E <- c(LETTERS[1:3], NA, "D")

network_types <- c(
  "advice",
  "leader",
  "leadMotivate",
  "leadEmotion",
  "influence",
  "trust",
  "dislike"
)
ntypes <- length(network_types)

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
    network = as.integer(stringr::str_extract(Code, "(?<=[_])[0-9]+"))
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
    network = as.integer(stringr::str_extract(Code, "(?<=[_])[0-9]+"))
  ) 

edgelists <- bind_rows(edgelist_post, edgelist_pre) %>%
  
  # Individuals could have answered NA, so we have to remove these
  filter(!is.na(alter)) %>%
  
  # Merging group size. We use full to make sure that we have all the data
  full_join(groups_sizes) 


# Now, creating the adjmat
networks <- vector("list", 2)
names(networks) <- c(1, 3)

for (s in c(1, 3)) {
  
  # Nested at the type level
  networks[[as.character(s)]] <- vector("list", ntypes)
  names(networks[[as.character(s)]]) <- network_types
  
  for (nt in seq_along(network_types)) {
    
    # Nested at the group level
    networks[[as.character(s)]][[nt]] <- vector("list", nrow(groups_sizes))
    names(networks[[as.character(s)]][[nt]]) <- groups_sizes$group
    
    for (gid in as.character(groups_sizes$group)) {
      
      # Group features: Size, IDS and matrix
      n   <- filter(groups_sizes, group == gid)$n
      ids <- LETTERS[1:n]
      M   <- matrix(0L, ncol = n, nrow = n, dimnames = list(ids, ids))
  
      # Filtering the data
      g <- filter(edgelists, survey == s, network == nt, group == gid)
      
      if (!nrow(g)) {
        message("Group ", gid, " survey ", s, " network `",
                network_types[nt],"` done (but had no edges).")
        networks[[as.character(s)]][[network_types[nt]]][[gid]] <- M
        next
      }
      
      # Keep complete cases
      g <- g[complete.cases(g),]
      
      if (nrow(g) == 0) {
        networks[[as.character(s)]][[network_types[nt]]][[gid]] <- M
        next
      }
      
      # I've observed some errors in the data
      test <- which(!(g$alter %in% ids))
      if (length(test)) {
        
        # Throwing a warning
        warning(
          "Group ", g$group[1],
          " had reports on alters no included in their groups",
          " in wave ", s,
          sprintf(" (%i vs %i).", n, length(unique(g$alter))),
          call.=FALSE)
        
        # Resizing the dataset
        g <- filter(g, alter %in% ids)
        
      }
      
      # Is there any data left?
      if (!nrow(g)) {
        message("Group ", gid, " survey ", s, " network `",
                network_types[nt],"` done (but had no edges).")
        networks[[as.character(s)]][[network_types[nt]]][[gid]] <- M
        next
      }
      
      
      # Adding 1s
      for (i in 1:nrow(g))
        M[g$ego[i], g$alter[i]] <- 1L
      
      # Saving the result
      networks[[as.character(s)]][[network_types[nt]]][[gid]] <- M
      
      message("Group ", gid, " survey ", s, " network `",
              network_types[nt],"` done.")
      
      
    }
  }
}

# Saving the object
saveRDS(networks, "data/networks_truth.rds")

