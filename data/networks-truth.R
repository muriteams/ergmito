library(dplyr)
library(tidyselect)
library(tidyr)
library(magrittr)
library(igraph)

# Reading in Data
dat_nominations <- readxl::read_excel("data-raw/Study1_Dyad variables excel.xlsx")
dat_groupsize   <- readr::read_csv("data-raw/Study1_Group sizes.csv")

# Question, s1LeaderSumBin and s3LeaderSumBin?
nettypes <- dat_nominations %>% select(matches("^s[0-9]q[0-9]")) %>% colnames

# Reshaping long
dat_nominations <- dat_nominations %>%
  gather(tie, present, matches("^s[0-9]")) %>%
  filter(present == 1) %>%
  select(-present) %>%
  mutate(Ego = gsub("^[0-9]+", "", EgoPID))

# Creating list of nested networks
networks_truth <- vector("list", nrow(dat_groupsize))
names(networks_truth) <- dat_groupsize$Group

# Looping through groups. Notice that we use the group number as a name since
# 5 is missing.
for (g in dat_groupsize$Group) {
  
  # Creating the subvector to store all network types
  gchar <- as.character(g)
  
  networks_truth[[gchar]] <- vector("list", length = length(nettypes))
  names(networks_truth[[gchar]]) <- nettypes
  
  # Looping through network types
  for (nt in nettypes) {
    
    # From the original data, we get the network type and group name
    networks_truth[[gchar]][[nt]] <- dat_nominations %>%
      filter(tie == nt, GroupID == g) %>%
      select(Ego, Alter) %>%
      # Turn it into an igraph object passing vertices to make sure that we don't
      # miss any node.
      graph_from_data_frame(
        vertices = data.frame(
          id = LETTERS[1:dat_groupsize[dat_groupsize$Group == g, 2][[1]]]
        )
      ) %>%
      # And turn it into
      as_adj(sparse=FALSE)
    
  }
  
}

saveRDS(networks_truth, "data/networks-truth.rds")
