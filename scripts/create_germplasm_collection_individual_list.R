# cranberryHistoricalGBS
# 
# Prepare a file that contains the full sample names of only those individuals
# in the germplasm collection
# 

# Directories and packages
library(tidyverse)
library(readxl)

proj_dir <- here::here()
input_dir <- file.path(proj_dir, "input")

# Cranberry dir
cran_dir <- strsplit(proj_dir, "/")[[1]] %>%
  {.[seq_len(which(. == "CranberryLab"))]} %>%
  paste0(collapse = "/")


# Designate the key file name
keyfile_name <- file.path(input_dir, "cranberry_gbs_unique_keys_resolved_duplicates.txt")




# Create the file of individuals ------------------------------------------

# Read in the edit keyfile
keyfile <- read_tsv(file = keyfile_name) %>%
  # Filter out samples that start with CNJ or P
  filter(str_detect(FullSampleName, "^CNJ|^P[0-9]{1,}", negate = TRUE)) %>%
  distinct(FullSampleName) %>%
  arrange(FullSampleName)

# Save
write_tsv(x = keyfile, file.path(input_dir, "cranberry_gbs_germplasm_collection_individuals.txt"),
          col_names = FALSE)
