# cranberryHistoricalGBS
# 
# Consolidation of the keyfile
# 

# This script will edit the keyfile for use in the TASSEL5 pipeline. It will:
#   1. Remove blank samples
#   2. Give duplicated samples the same sample name.
#   

# Directories and packages
library(tidyverse)

proj_dir <- here::here()
input_dir <- file.path(proj_dir, "input")

# Read in the keyfile
key <- read_tsv(file = file.path(input_dir, "cranberry_gbs_all_keys.txt"))

# Remove blanks
key1 <- filter(key, !str_detect(toupper(DNASample), "BLANK")) %>%
  mutate(SeedLot = ifelse(is.na(SeedLot), DNASample, SeedLot),
         SeedLot = str_remove_all(SeedLot, " ")) %>%
  group_by(SeedLot) %>%
  mutate(FullSampleName = paste0(SeedLot, ":", LibraryPrepID[1])) %>%
  ungroup()

# Save
write_tsv(x = key1, path = file.path(input_dir, "cranberry_gbs_unique_keys.txt"))
