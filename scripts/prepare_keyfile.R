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
         SeedLot = str_remove_all(SeedLot, " "),
         Pedigree = ifelse(str_detect(SeedLot, "^P[0-9]{1,}"), "BG95xGH1", Pedigree)) %>%
  group_by(SeedLot) %>%
  mutate(FullSampleName = paste0(SeedLot, ":", LibraryPrepID[1])) %>%
  ungroup()

# Save
write_tsv(x = key1, path = file.path(input_dir, "cranberry_gbs_unique_keys.txt"))


# Vector of parent names
parents <- c("MullicaQueen", "CrimsonQueen", "Stevens", "GH1", "BG95")


# Create separate text lists for each bi-parental family
family_keys <- key1 %>% 
  filter(str_detect(SeedLot, "CNJ|^P[0-9]{1,}") | SeedLot %in% parents) %>%
  select(individual = SeedLot, pedigree = Pedigree, fullSampleName = FullSampleName) %>%
  # Correct GH1 pedigree
  mutate(pedigree = ifelse(pedigree == "wild", NA, pedigree)) %>%
  distinct() %>%
  mutate(family = case_when(
    str_detect(individual, "CNJ02") | individual %in% c("MullicaQueen", "CrimsonQueen") ~ "CNJ02",
    str_detect(individual, "CNJ04") | individual %in% c("Stevens") ~ "CNJ04",
    str_detect(individual, "^P") | individual %in% c("GH1", "BG95") ~ "GRYG")
  ) %>%
  bind_rows(., { subset(., individual == "MullicaQueen", c(individual, pedigree, fullSampleName)) }) %>% 
  mutate(family = ifelse(individual == "MullicaQueen" & is.na(family), "CNJ04", family),
         individual = as.factor(individual),
         individual = fct_relevel(individual, parents)) %>%
  arrange(individual)

# Write this
write_tsv(x = family_keys, path = file.path(input_dir, "cranberry_gbs_bp_family_metadata.txt"))


# Split by family and write files
family_keys_split <- family_keys %>% 
  split(.$family) %>%
  map(select, fullSampleName)

for (i in seq_along(family_keys_split)) {
  filename <- file.path(input_dir, paste0("family_", names(family_keys_split)[i], "_individuals.txt"))
  write_tsv(x = family_keys_split[[i]], path = filename, col_names = FALSE)
}
