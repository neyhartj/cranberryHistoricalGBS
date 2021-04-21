## cranberryHistoricalGBS
## 
## Analyze the results of imputation
## 


# Directories and packages
library(tidyverse)
library(readxl)
library(vcfR)

proj_dir <- here::here()
pipeline_dir <- proj_dir
impute_dir <- file.path(proj_dir, "imputation/beagle_imputation/")

# Cranberry dir
cran_dir <- strsplit(proj_dir, "/")[[1]] %>%
  {.[seq_len(which(. == "CranberryLab"))]} %>%
  paste0(collapse = "/")

# Read in data ---------------------------------------------------------

# Read in the germplasm collection metadata
filename <- file.path(cran_dir, "Breeding/prior2021/Populations/GermplasmCollection/germplasmCollectionMetadata.xlsx")
metadata <- read_excel(path = filename, na = c("", "NA"))

# Subset those with marker data; select relevant columns
metadata1 <- metadata %>%
  filter(!is.na(marker_sample_name)) %>%
  mutate(genetic_id = str_split(genetic_id, ", ")) %>%
  select(marker_sample_name, individual, formatted_name, genetic_id, variety_designation, variety_selection_name,
         category, contains("origin"))



# Read in the VCF
filename <- file.path(impute_dir, "cranberryGBS_snps_to_impute_filtered_firstround.vcf.gz")
vcf_in <- read.vcfR(file = filename)

# Change column names
# Rename columns
colnames(vcf_in@gt) <- toupper(map_chr(str_split(colnames(vcf_in@gt), ":"), 1))

# Use the metadata to select columns
vcf_in1 <- vcf_in[, c("FORMAT", intersect(sort(metadata1$marker_sample_name), colnames(vcf_in@gt)))]

# Rename the columns with formatted sample names
samples_rename_df <- tibble(marker_sample_name = colnames(vcf_in1@gt)) %>%
  left_join(., select(metadata1, marker_sample_name, formatted_name)) %>%
  mutate(formatted_name = ifelse(marker_sample_name == "FORMAT", "FORMAT", formatted_name)) %>%
  # For samples belonging to the same marker sample name, pick the first
  group_by(marker_sample_name) %>%
  slice(1) %>%
  ungroup() %>%
  # sort
  mutate(marker_sample_name = factor(marker_sample_name, levels = colnames(vcf_in1@gt))) %>%
  arrange(marker_sample_name)

colnames(vcf_in1@gt) <- samples_rename_df$formatted_name


# Create genotype matrices / hapmap
# Pull out the genotype data
gt <- vcf_in1@gt[,-1]
# Add marker names as row
row.names(gt) <- vcf_in1@fix[,3]

# Transpose
gt2 <- t(gt)


# Create a genotype matrix
# Sum the number of reference alleles; subtract 1
geno_mat <- str_extract(string = gt2, pattern = "[01]/[01]") %>% 
  str_split(string = ., pattern = "/") %>% 
  map(as.numeric) %>% 
  map_dbl(sum) %>%
  matrix(data = ., nrow = nrow(gt2), ncol = ncol(gt2), byrow = FALSE, dimnames = dimnames(gt2))



# Determine sample duplication --------------------------------------------

# Calculate the proportion of identical genotypes
sample_geno_identical <- t(combn(x = row.names(geno_mat), m = 2)) %>%
  as.data.frame() %>%
  rename(individual = V1, individual2 = V2) %>%
  mutate(percent_identity = map2_dbl(individual, individual2, ~mean(geno_mat[.x,] == geno_mat[.y,], na.rm = TRUE)))
  

# Tidy
sample_geno_identical1 <- sample_geno_identical %>%
  # Calculate the pairwise number of missing markers
  mutate(n = map2_dbl(individual, individual2, ~sum(colSums(is.na(geno_mat[c(.x, .y),])) == 0)))

# Range in correlation
range(sample_geno_identical1$percent_identity)
hist(sample_geno_identical1$percent_identity)

# Subset pairs above 0.90
sample_geno_similar <- sample_geno_identical1 %>%
  filter(percent_identity > 0.95) %>%
  arrange(individual, individual2) %>%
  group_by_all() %>%
  do({
    row <- .
    
    row1 <- subset(metadata1, formatted_name == row$individual)
    row2 <- subset(metadata1, formatted_name == row$individual2)
    alias1 <- row1$genetic_id[[1]]
    alias2 <- row2$genetic_id[[1]]
  
    # Is individual 1 in alias 2 or vice versa?
    common_names <- union(intersect(row1$individual, alias2), intersect(row2$individual, alias1))
  
    # Pull out the variety selection names
    var1 <- row1$variety_selection_name
    var2 <- row2$variety_selection_name
    
    # Pull out origins
    origin1 <- row1$origin_name
    origin2 <- row2$origin_name
    
    # Get missingness for each individual
    missing1 <- mean(is.na(geno_mat[row1$formatted_name,]))
    missing2 <- mean(is.na(geno_mat[row2$formatted_name,]))
    
    
    # Return a tibble
    tibble(aliases = paste0(common_names, collapse = ", "),
           individual1_varname = var1, individual2_varname = var2,
           individual1_origin = origin1, individual2_origin = origin2,
           individual1_missing = missing1, individual2_missing = missing2)
    
  }) %>% ungroup() %>%
  # Designate the individual with the highest missingness to remove
  mutate(remove = ifelse(individual1_missing > individual2_missing, individual, individual2))

# Write to disk
write_csv(x = sample_geno_similar, path = file.path(proj_dir, "stats/sample_similarity_r90.csv"))




## Look at confirmed or suspected duplicates

# First create a list of duplicate groups
suspected_duplicate_group_list <- metadata1 %>%
  filter(map_lgl(genetic_id, ~!all(is.na(.)))) %>%
  mutate(genetic_id2 = map2(individual, genetic_id, union) %>% map(sort)) %>%
  pull(genetic_id2) %>%
  unique() %>%
  # Create list of data.frames
  map(~left_join(tibble(individual = .x), select(metadata1, formatted_name, individual), by = "individual") %>%
        mutate(genotyped = formatted_name %in% row.names(geno_mat)))

# For each group, create a pairwise similarity matrix
suspected_duplicate_similarity <- map(suspected_duplicate_group_list, ~{
  .x <- .x$individual
  # Change sample names
  sample_geno_identical_mat <- sample_geno_identical %>%
    left_join(., select(metadata1, formatted_name, individual), by = c("individual" = "formatted_name")) %>% 
    left_join(., select(metadata1, formatted_name, individual), by = c("individual2" = "formatted_name")) %>%
    select(individual1 = individual.y, individual2 = individual.y.y, percent_identity) %>%
    filter_at(vars(contains("individual")), all_vars(. %in% .x))
  
  # Build a similarity matrix
  indiv_names <- unique(unlist(select(sample_geno_identical_mat, contains("indiv"))))
  dist <- matrix(NA, nrow = length(indiv_names), ncol = length(indiv_names), 
                 dimnames = list(indiv_names, indiv_names)) %>% as.dist()
  dist[] <- sample_geno_identical_mat$percent_identity
  distmat <- as.matrix(dist)
  diag(distmat) <- 1
  
  # Return the similarity matrix
  return(distmat)
  
})

# Filter the null matrices
suspected_duplicate_similarity1 <- suspected_duplicate_similarity[!sapply(suspected_duplicate_similarity, is_empty)]



## High similar genotypes were inspected by hand; a decision was made as to which
## genotypes to discard and which to merge. 
## 

# Read in the excel sheet with duplicate resolution information.
duplicate_resolution <- read_excel(path = file.path(proj_dir, "stats/sample_similarity_r90.xlsx"))

# Create a vector of genotypes to remove from the marker matrix; those pairs of individuals
# designated with "merge" will not be included in this list
duplicate_resolution_merge <- duplicate_resolution %>% 
  filter(remove == "merge") 

individuals_to_merge <- duplicate_resolution_merge %>% 
  select(individual, individual2) %>%
  unlist() %>% 
  unique()

# Create a data.frame of the formatted name to merge, it's marker sample name,
# and the marker sample name of the final merged individual
individuals_to_merge_list <- list()
individuals_to_merge1 <- individuals_to_merge

while (length(individuals_to_merge1) > 0) {
  indiv_i <- individuals_to_merge1[1]
  # Find all mentions of this individual in the duplicate resolution
  indiv_i_connections <- indiv_i_connections1 <- duplicate_resolution_merge %>% 
    filter_at(vars(individual, individual2), any_vars(. == indiv_i)) %>% 
    select(individual, individual2) %>% 
    unlist() %>% 
    unique()
  
  # Find any connections with these
  change <- TRUE
  while (change) {
    indiv_i_connections <- duplicate_resolution_merge %>% 
      filter_at(vars(individual, individual2), any_vars(. %in% indiv_i_connections)) %>% 
      select(individual, individual2) %>% 
      unlist() %>% 
      unique()
    change <- any(!indiv_i_connections %in% indiv_i_connections1)
    indiv_i_connections1 <- indiv_i_connections
  }
  
  # add to the list
  individuals_to_merge_list[[indiv_i]] <- indiv_i_connections1
  
  # Remove the connections from the merge list
  individuals_to_merge1 <- setdiff(individuals_to_merge1, individuals_to_merge_list[[indiv_i]])
  
}

# The name of each list element will be the name of the merged sample
individuals_to_merge_list1 <- individuals_to_merge_list %>%
  imap(~tibble(formatted_name = .x, merged_name = .y)) %>%
  # add marker sample names
  map_df(~left_join(., select(metadata1, formatted_name, marker_sample_name)) %>%
           rename(former_marker_sample_name = marker_sample_name) %>%
           left_join(., select(metadata1, formatted_name, marker_sample_name), by = c("merged_name" = "formatted_name")) )



individuals_to_remove <- setdiff(unique(subset(duplicate_resolution, remove != "merge", remove, drop = TRUE)), individuals_to_merge)

# Use the individual names to lookup sample names for removal
individuals_to_remove_marker_names <- metadata1 %>%
  filter(formatted_name %in% individuals_to_remove) %>% 
  pull(marker_sample_name)


## Edit the key file to remove the individuals to remove and merge the individuals
## to merge using a common (first) name
keyfile <- read_tsv(file = file.path(proj_dir, "input/cranberry_gbs_unique_keys.txt"))

# Make sure the marker sample name get the library prep ID
individuals_to_merge_list2 <- individuals_to_merge_list1 %>%
  left_join(., distinct(keyfile, SeedLot, LibraryPrepID), by = c("marker_sample_name" = "SeedLot")) %>%
  rename(library_prep_ID = LibraryPrepID)


# Remove the individuals to remove
keyfile1 <- keyfile %>%
  filter(! SeedLot %in% individuals_to_remove_marker_names) %>%
  # Merge individuals using the first sample name
  left_join(., select(individuals_to_merge_list2, contains("sample_name"), library_prep_ID), by = c("SeedLot" = "former_marker_sample_name")) %>%
  # Edit the seedlot names
  mutate(SeedLot = ifelse(is.na(marker_sample_name), SeedLot, marker_sample_name),
         FullSampleName = ifelse(is.na(marker_sample_name), FullSampleName, paste0(SeedLot, ":", library_prep_ID))) %>%
  select(-marker_sample_name, -library_prep_ID)

# Save this new keyfile
write_tsv(x = keyfile1, path = file.path(proj_dir, "input/cranberry_gbs_unique_keys_resolved_duplicates.txt"))





