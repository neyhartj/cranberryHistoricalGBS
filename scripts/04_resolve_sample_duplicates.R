## cranberryHistoricalGBS
## 
## Analyze the results of imputation
## 


# Directories and packages
library(igraph)
library(tidyverse)
library(readxl)
library(vcfR)


proj_dir <- here::here()
pipeline_dir <- proj_dir
impute_dir <- file.path(proj_dir, "imputation/beagle_imputation/")
snps_dir <-  file.path(proj_dir, "/snps/")

# Cranberry dir
cran_dir <- strsplit(proj_dir, "/")[[1]] %>%
  {.[seq_len(which(. == "CranberryLab"))]} %>%
  paste0(collapse = "/")

# Minimum matching percent to merge
min_prop_matching <- 0.99


# A function that, given a list of individuals and a genotype matrix (n x m),
# will calculate the pairwise proportion of matching genotype calls
pairwise_geno_matching <- function(geno.mat, sample.names, table = TRUE) {
  # Generate pairs of sample names
  sample_name_pairs <- as.data.frame(t(combn(x = sample.names, m = 2)))
  names(sample_name_pairs) <- c("sample1", "sample2")
  sample_name_pairs <- as_tibble(sample_name_pairs)
  sample_name_pairs$prop_matching <- as.numeric(NA)
  sample_name_pairs$matching_table <- list(NULL)
  
  pb <- progress::progress_bar$new(total = nrow(sample_name_pairs))
  
  for (i in seq_len(nrow(sample_name_pairs))) {
    s1 <- sample_name_pairs$sample1[i]
    s2 <- sample_name_pairs$sample2[i]
    
    if (table) {
      matching_table <- table(geno.mat[s1,], geno.mat[s2,])
      # Add to the df
      sample_name_pairs$prop_matching[i] <- sum(diag(prop.table(matching_table)))
      sample_name_pairs$matching_table[[i]] <- matching_table
      
    } else {
      sample_name_pairs$prop_matching[i] <- mean(geno.mat[s1,] == geno.mat[s2,], na.rm = TRUE)
      
    }
    
    pb$tick()
    
  }
  
  # Return
  return(sample_name_pairs)
  
}



# Read in data ---------------------------------------------------------

# Read in the original keyfile
## Edit the key file to remove the individuals to remove and merge the individuals
## to merge using a common (first) name
keyfile <- read_tsv(file = file.path(proj_dir, "input/cranberry_gbs_unique_keys.txt"))


# Read in the VCF
vcf_in <- read.vcfR(file = file.path(snps_dir, "cranberryGBS_production_snps_allUniqueKeys_filtered.vcf.gz"))

# Read in germplasm metadata
germplasm_metadata <- read_excel(path = file.path(cran_dir, "Breeding/Germplasm/all_germplasm_metadata.xlsx"), col_types = "text") %>%
  mutate_all(parse_guess)

# Read in the genotyped germplasm metadata
genotyped_germplasm_metadata <- read_excel(path = file.path(cran_dir, "Genotyping/MarkerDatabase/genotyped_germplasm_database.xlsx")) %>%
  filter(marker_platform == "GBS")

# select relevant columns from the germplasm metadata
germplasm_metadata1 <- germplasm_metadata %>%
  select(individual, unique_id, individual_name_use, alias, scar_alias, variety_designation, variety_selection_name,
         maternal_parent, paternal_parent) %>%
  # Split up any columns with commas
  mutate(scar_alias = str_split(scar_alias, ", ")) %>%
  # Subset those that have been genotyped (i.e. have a marker sample name)
  inner_join(., genotyped_germplasm_metadata)
  


# Edit data ---------------------------------------------------------------

# Extract the sample names; cut out the full name
vcf_orig_sample_names <- map_chr(str_split(colnames(vcf_in@gt), ":"), 1)

# Match with the genotyped germplasm database - this will result in some duplication
sample_names_match <-  tibble(sample = vcf_orig_sample_names) %>%
  left_join(., genotyped_germplasm_metadata, by = c("sample" = "SeedLot")) %>%
  # If the individual name is NA, replace with the sample name
  mutate(individual = ifelse(is.na(individual), sample, individual)) %>%
  arrange(individual)

# Replace the column names with these names
vcf_in1 <- vcf_in
colnames(vcf_in@gt) <- vcf_orig_sample_names

# Create genotype matrices / hapmap
# Pull out the genotype data
gt <- vcf_in1@gt[,-1]
# Add marker names as row
row.names(gt) <- vcf_in1@fix[,3]

# Transpose
gt2 <- t(gt)

# Create a genotype matrix
# Sum the number of reference alleles
geno_mat <- str_extract(string = gt2, pattern = "[01]/[01]") %>% 
  str_split(string = ., pattern = "/") %>% 
  map(as.numeric) %>% 
  map_dbl(sum) %>%
  matrix(data = ., nrow = nrow(gt2), ncol = ncol(gt2), byrow = FALSE, dimnames = dimnames(gt2))




# 
# 

# 
# 
# 
# 
# 
# # Change column names in the VCF - split by colon and take the first in the split
# colnames(vcf_in@gt) <- toupper(map_chr(str_split(colnames(vcf_in@gt), ":"), 1))
# 
# # Use the metadata to select columns
# vcf_in1 <- vcf_in[, c("FORMAT", intersect(sort(metadata1$marker_sample_name), colnames(vcf_in@gt)))]
# 
# # Rename the columns with formatted sample names
# samples_rename_df <- tibble(marker_sample_name = colnames(vcf_in1@gt)) %>%
#   left_join(., select(metadata1, marker_sample_name, formatted_name)) %>%
#   mutate(formatted_name = ifelse(marker_sample_name == "FORMAT", "FORMAT", formatted_name)) %>%
#   # For samples belonging to the same marker sample name, pick the first
#   group_by(marker_sample_name) %>%
#   slice(1) %>%
#   ungroup() %>%
#   # sort
#   mutate(marker_sample_name = factor(marker_sample_name, levels = colnames(vcf_in1@gt))) %>%
#   arrange(marker_sample_name)
# 
# colnames(vcf_in1@gt) <- samples_rename_df$formatted_name



# Examine explicit sample duplication -------------------------------------

# These are samples with known alias (i.e. CrimsonQueen = NJS98-23)
explicit_duplicates <- sample_names_match %>%
  group_by(individual) %>%
  filter(n() > 1) %>%
  ungroup()

# Split by individual and compare
explicit_duplicates_similarity <- explicit_duplicates %>%
  split(.$individual) %>%
  map("sample") %>%
  map_df(~pairwise_geno_matching(geno.mat = geno_mat, sample.names = .x))
  
# All of these individuals had matching > 0.99; merge them in the resolved keyfile
explicit_duplicates_rename <- explicit_duplicates %>%
  select(original_sample_name = sample, new_sample_name = individual)



# Examine implicit duplicates ---------------------------------------------

# These are duplicates with low- or high-confidence aliases that could be merged
implicit_duplicates <- germplasm_metadata1 %>%
  filter(!map_lgl(scar_alias, ~all(is.na(.))) | ! is.na(alias) | !is.na(variety_designation) | !is.na(variety_selection_name))

## First look at aliases
implicit_duplicates_aliases <- implicit_duplicates %>%
  filter(!is.na(alias)) %>%
  select(individual, unique_id, individual_name_use, SeedLot, alias) %>%
  # Remove cases where the SeedLot is a different spelling of the alias
  mutate(SeedLot = toupper(str_replace_all(neyhart::str_add_space(SeedLot), " ", "_"))) %>%
  filter(SeedLot != alias)

# Find any aliases that are in the genotyped metadata
any(implicit_duplicates_aliases$alias %in% genotyped_germplasm_metadata$SeedLot)

# No aliases are individual names


## Next look at scar aliases (confidence level 1)
implicit_duplicates_scars <- implicit_duplicates %>%
  filter(!map_lgl(scar_alias, ~all(is.na(.)))) %>%
  select(individual, unique_id, individual_name_use, SeedLot, scar_alias) %>%
  # Add the individual to each scar list
  mutate(scar_alias = map2(scar_alias, individual, c))

# Merge groups of SCAR aliases if they share a member
scars_list <- implicit_duplicates_scars$scar_alias
scar_adjacency_matrix <- sapply(X = scars_list, FUN = function(x) sapply(scars_list, function(y) length(intersect(x,y)) > 0))
# Groups of scar aliases to merge
scar_alias_list_groups <- groups(components(graph_from_adjacency_matrix(scar_adjacency_matrix)))

implicit_duplicates_scars_merged <- lapply(X = scar_alias_list_groups, FUN = function(x) sort(unique(unlist(scars_list[x]))))
  

# Find groups of scar aliases
implicit_duplicates_scars_groups <- tibble(group = paste0("group", names(implicit_duplicates_scars_merged)), scar_alias = implicit_duplicates_scars_merged) %>%
  unnest(scar_alias) %>%
  rename(individual = scar_alias) %>%
  # Add seedlot names
  left_join(., select(genotyped_germplasm_metadata, individual, SeedLot), by = "individual") %>%
  filter(!is.na(SeedLot)) %>%
  # Remove groups of 1
  group_by(group) %>%
  filter(n_distinct(SeedLot) > 1) %>%
  ungroup() %>%
  # REmove any duplicates
  distinct() %>%
  # split
  split(.$group) 

# Calculate pairwise similarity
implicit_duplicates_scars_similarity <- implicit_duplicates_scars_groups %>%
  map("SeedLot") %>%
  map_df(~pairwise_geno_matching(geno.mat = geno_mat, sample.names = .x)) %>%
  # Add sample names
  left_join(., select(genotyped_germplasm_metadata, individual, SeedLot), by = c("sample1" = "SeedLot")) %>%
  left_join(., select(genotyped_germplasm_metadata, individual, SeedLot), by = c("sample2" = "SeedLot")) %>%
  select(individual1 = individual.x, SeedLot1 = sample1, individual2 = individual.y, SeedLot2 = sample2, names(.))

# Save this
implicit_duplicates_scars_similarity %>%
  select(-matching_table) %>%
  write_csv(x = ., file = "stats/sample_genotype_matching_scarSimilarityLists.csv")




## Next look at variety designations (confidence level 2)
implicit_duplicates_varietyDesignation <- implicit_duplicates %>%
  filter(!is.na(variety_designation)) %>%
  select(individual, unique_id, individual_name_use, SeedLot, variety_designation) %>%
  filter(!is.na(SeedLot)) %>%
  # Remove groups of 1
  group_by(variety_designation) %>%
  filter(n_distinct(SeedLot) > 1) %>%
  ungroup() %>% # split
  split(.$variety_designation) 

# Calculate pairwise similarity
implicit_duplicates_varietyDesignation_similarity <- implicit_duplicates_varietyDesignation %>%
  map("SeedLot") %>%
  map_df(~pairwise_geno_matching(geno.mat = geno_mat, sample.names = .x)) %>%
  # Add sample names
  left_join(., select(genotyped_germplasm_metadata, individual, SeedLot), by = c("sample1" = "SeedLot")) %>%
  left_join(., select(genotyped_germplasm_metadata, individual, SeedLot), by = c("sample2" = "SeedLot")) %>%
  # Add variety designation
  left_join(., select(germplasm_metadata1, individual.x = individual, variety_designation)) %>%
  select(variety_designation, individual1 = individual.x, SeedLot1 = sample1, individual2 = individual.y, 
         SeedLot2 = sample2, names(.))

# Save this
implicit_duplicates_varietyDesignation_similarity %>%
  select(-matching_table) %>%
  write_csv(x = ., file = "stats/sample_genotype_matching_varietyDesignationLists.csv")





## Next look at variety selection names (confidence level 3)
implicit_duplicates_varietySelectionName <- implicit_duplicates %>%
  filter(!map_lgl(variety_selection_name, ~all(is.na(.)))) %>%
  select(individual, unique_id, individual_name_use, SeedLot, variety_selection_name) %>%
  filter(!is.na(SeedLot), SeedLot %in% row.names(geno_mat)) %>%
  # Remove groups of 1
  group_by(variety_selection_name) %>%
  filter(n_distinct(SeedLot) > 1) %>%
  ungroup() %>% # split
  split(.$variety_selection_name) 

# Calculate pairwise similarity
implicit_duplicates_varietySelectionName_similarity <- implicit_duplicates_varietySelectionName %>%
  map("SeedLot") %>%
  map(unique) %>%
  map_df(~pairwise_geno_matching(geno.mat = geno_mat, sample.names = .x)) %>%
  # Add sample names
  left_join(., select(genotyped_germplasm_metadata, individual, SeedLot), by = c("sample1" = "SeedLot")) %>%
  left_join(., select(genotyped_germplasm_metadata, individual, SeedLot), by = c("sample2" = "SeedLot")) %>%
  # Add variety designation
  left_join(., select(germplasm_metadata1, individual.x = individual, variety_selection_name)) %>%
  select(variety_selection_name, individual1 = individual.x, SeedLot1 = sample1, individual2 = individual.y, 
         SeedLot2 = sample2, names(.)) %>%
  # Arrange by matching
  arrange(desc(prop_matching))

# Save this
implicit_duplicates_varietySelectionName_similarity %>%
  select(-matching_table) %>%
  write_csv(x = ., file = "stats/sample_genotype_matching_varietySelectionNameLists.csv")




# Examine all potential duplicates ----------------------------------------

# Calculate all pairwise within the germplasm collection
sample_names <- genotyped_germplasm_metadata %>%
  filter(str_detect(individual, "^CNJ02-1|^CNJ04-2|^P[0-9]", negate = TRUE)) %>%
  pull(SeedLot) %>%
  intersect(., row.names(geno_mat))

all_sample_pairwise_matching <- pairwise_geno_matching(geno.mat = geno_mat, sample.names = sample_names, table = FALSE) %>%
  left_join(., select(genotyped_germplasm_metadata, individual, SeedLot), by = c("sample1" = "SeedLot")) %>%
  left_join(., select(genotyped_germplasm_metadata, individual, SeedLot), by = c("sample2" = "SeedLot")) %>%
  select(individual1 = individual.x, SeedLot1 = sample1, individual2 = individual.y, SeedLot2 = sample2, prop_matching)

# Remove explicit duplicates
all_sample_pairwise_matching1 <- all_sample_pairwise_matching %>%
  mutate(group = map2_chr(SeedLot1, SeedLot2, ~paste0(sort(c(.x, .y)), collapse = ":"))) %>%
  filter(! group %in% map_chr(split(explicit_duplicates$sample, explicit_duplicates$individual), ~paste0(sort(.x), collapse = ":"))) %>%
  select(-group)


# Create a combined table of matching statistics --------------------------

complete_geno_matching_tables <- bind_rows(
  mutate(implicit_duplicates_scars_similarity, putative_matching_name = as.character(NA), reason_for_matching = "SCAR"),
  mutate(rename(implicit_duplicates_varietyDesignation_similarity, putative_matching_name = variety_designation), reason_for_matching = "varietyDesignation"),
  mutate(rename(implicit_duplicates_varietySelectionName_similarity, putative_matching_name = variety_selection_name), reason_for_matching = "varietySesignationName"),
  mutate(all_sample_pairwise_matching1, reason_for_matching = "None", putative_matching_name = as.character(NA))
) %>%
  select(putative_matching_name, reason_for_matching, names(.), -matching_table) %>%
  filter(prop_matching > min_prop_matching) %>%
  arrange(reason_for_matching, desc(prop_matching))

# Remove duplicates with the following priorties:
# SCAR > varietyDesignation > VarietySelectionName > None
complete_geno_matching_tables1 <- complete_geno_matching_tables %>%
  mutate(group = map2_chr(individual1, individual2, ~paste0(sort(c(.x, .y)), collapse = ":")),
         reason_for_matching = fct_inorder(reason_for_matching) %>% fct_relevel(., "None", after = Inf)) %>%
  split(.$group) %>%
  map_df(~arrange(., reason_for_matching) %>% head(1)) %>%
  select(-group) %>%
  arrange(reason_for_matching, desc(prop_matching))

# Save this
write_csv(x = complete_geno_matching_tables, file = "stats/sample_genotype_matching_combined.csv", na = "")




# Merge highly similar genotypes based on a priori evidence only ----------

matching_sample_to_merge <- complete_geno_matching_tables1 %>%
  filter(reason_for_matching != "None")

# Find common groups of SeedLots to merge
merge_lists <- map2(matching_sample_to_merge$SeedLot1, matching_sample_to_merge$SeedLot2, c)
merge_lists_adjacency_matrix <- sapply(X = merge_lists, FUN = function(x) sapply(merge_lists, function(y) length(intersect(x,y)) > 0))
# Groups of scar aliases to merge
merge_lists_groups <- groups(components(graph_from_adjacency_matrix(merge_lists_adjacency_matrix)))

seedLots_merged <- lapply(X = merge_lists_groups, FUN = function(x) sort(unique(unlist(merge_lists[x])))) %>%
  # Sort to prefer variety names
  lapply(sort, decreasing = T)

# Create a table of the name to use and previous names (with individual names)
seedLots_merged_df <- seedLots_merged %>%
  map_df(~tibble(new_sample_name = .x[1], original_sample_name = .x)) %>%
  left_join(., select(genotyped_germplasm_metadata, SeedLot, individual), by = c("original_sample_name" = "SeedLot")) %>%
  rename(original_sample_name_individual = individual)


# Combine this with the explicit duplicates
samples_to_merge <- bind_rows(explicit_duplicates_rename, seedLots_merged_df) %>%
  # Covert hyphens to underscores
  mutate(new_sample_name = str_replace_all(new_sample_name, "-", "_"))
  
# Save a table
write_csv(x = samples_to_merge, file = "input/sample_names_merged.csv", na = "")



# Rename samples in the keyfile -------------------------------------------

# Remove duplicates
samples_to_merge1 <- samples_to_merge %>%
  distinct(original_sample_name, new_sample_name)

# First separate entries that will not be renamed
keyfile_entries_unchanged <- keyfile %>%
  filter(! SeedLot %in% samples_to_merge1$original_sample_name)

# Next subset those  that will be renamed
keyfile_entries_changed <- keyfile %>%
  filter(SeedLot %in% samples_to_merge1$original_sample_name) %>%
  # Attach the new sample name
  left_join(., samples_to_merge1, by = c("SeedLot" = "original_sample_name")) %>%
  # Split by new sample name
  split(.$new_sample_name) %>%
  # Choose the full sample name the corresponds to the new sample name
  map_df(~{
    newFullSampleName <- str_subset(.x$FullSampleName, unique(.x$new_sample_name))
    if (length(newFullSampleName_search) == 0) {
      newFullSampleName <- paste0(unique(.x$new_sample_name), ":", unique(.x$LibraryPrepID)[1])
    } 
    mutate(.x, SeedLot = new_sample_name, FullSampleName = newFullSampleName)
  }) %>%
  select(-new_sample_name)

# Merge the keyfiles
new_keyfile <- bind_rows(keyfile_entries_unchanged, keyfile_entries_changed) %>%
  arrange(LibraryPlate, Flowcell, Lane, Row, Col)

# Save this new keyfile
write_tsv(x = new_keyfile, file = file.path(proj_dir, "input/cranberry_gbs_unique_keys_resolved_duplicates.txt"))


# Save a text file of germplasm collection individuals
# First remove mapping population individuals
mapping_population_seedLot_names <- genotyped_germplasm_metadata %>%
  filter(individual %in% subset(germplasm_metadata, category == "MP", individual, drop = TRUE)) %>%
  select(seed_lot = SeedLot)

new_keyfile %>%
  distinct(FullSampleName) %>%
  separate(FullSampleName, c("seed_lot", "pred_id"), sep = ":", remove = FALSE) %>%
  anti_join(., mapping_population_seedLot_names) %>%
  select(FullSampleName) %>%
  write_tsv(x = ., file = file.path(proj_dir, "input/cranberry_gbs_germplasm_collection_individuals.txt"), col_names = FALSE)
  

