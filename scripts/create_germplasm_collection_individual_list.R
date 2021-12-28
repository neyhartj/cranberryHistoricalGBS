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


# Read in population metadata
geno_pop_metadata <- read_excel(path = file.path(cran_dir, "Genotyping/MarkerDatabase/genotyped_germplasm_database.xlsx"))
pop_metadata <- read_excel(path = file.path(cran_dir, "Breeding/Germplasm/all_germplasm_metadata.xlsx"), 
                           col_types = "text") %>%
  mutate_all(parse_guess)



# Create the file of individuals ------------------------------------------
# Read in the edited keyfile
all_keys <- read_tsv(file = keyfile_name)

germplasm_keys <- all_keys %>%
  # Filter out samples that start with CNJ or P
  filter(str_detect(FullSampleName, "^CNJ|^P[0-9]{1,}", negate = TRUE)) %>%
  distinct(FullSampleName) %>%
  arrange(FullSampleName)

# Save
write_tsv(x = germplasm_keys, file.path(input_dir, "cranberry_gbs_germplasm_collection_individuals.txt"),
          col_names = FALSE)



# List of wild cranberry names --------------------------------------------

wild_germplasm_keys <- pop_metadata %>% 
  filter_at(vars(category, category2), any_vars(. == "Wild")) %>%
  select(marker_sample_name) %>%
  inner_join(., all_keys, by = c("marker_sample_name" = "SeedLot")) %>%
  distinct(FullSampleName) %>%
  arrange(FullSampleName)

# Save
write_tsv(x = wild_germplasm_keys, file.path(input_dir, "cranberry_gbs_wild_germplasm_individuals.txt"),
          col_names = FALSE)



# List of germplasm collection names minus wild ---------------------------

native_breeding_germplasm_keys <- pop_metadata %>% 
  filter(!is.na(marker_sample_name)) %>%
  filter_at(vars(category, category2), all_vars(! . %in% c("MP", "Wild"))) %>%
  select(marker_sample_name) %>%
  inner_join(., all_keys, by = c("marker_sample_name" = "SeedLot")) %>%
  distinct(FullSampleName) %>%
  arrange(FullSampleName)

write_tsv(x = native_breeding_germplasm_keys, file.path(input_dir, "cranberry_gbs_native_breeding_germplasm_individuals.txt"),
          col_names = FALSE)



