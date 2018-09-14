
########################################
######## CSS Analysis - MURI ###########
#######  march-april 2018 ###############


library(readr)

#open project "CSS Analyses"

# load("srdata/srnetworks.rda")

##########################################
### Step 1: create CSS networks ##########
##########################################

#read in data of participant perceptions of advice ties among their other teammates 
dat <- read_csv("data-raw/CSS_All.csv")

# Turning letters into factors
# so teammate letter (A,B, C etc) becomes a factor
dat$Letter <- as.factor(dat$Letter)


# first create list of matrices with complete info (missing)
MAT <- vector("list", nrow(dat))
for (i in 1:nrow(dat)) {
  # Step 1: create the matrix
  MAT[[i]] <- matrix(dat[i,-c(1:3)], ncol = 5, byrow = TRUE)
  
  # Step 2: set names
  groupid     <- gsub("[a-zA-Z]+", "",dat$PID[i])
  fancy_names <- sprintf("%s%s", groupid, LETTERS[1:5])
  dimnames(MAT[[i]])  <- list(fancy_names, fancy_names)
  
  # Step 3: Remove the 99
  MAT[[i]] <- MAT[[i]][
    1:dat$`Group Size`[[i]],
    1:dat$`Group Size`[[i]]
    ]
  
  diag(MAT[[i]]) <- NA
  
  ## at this point we have 178 complete matrices
  ## that contain rows of missing for the ego
  
  # Step 4: Remove individual
  # MAT[[i]] <- MAT[[i]][
  #   -as.integer(dat$Letter[[i]]),
  #   -as.integer(dat$Letter[[i]])
  #   ]
  
} 

names(MAT) <- dat$PID
