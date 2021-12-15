## Cranberry Historical GBS
## 
## Prepare data for AlphaFamImpute
## 

# Packages
library(GWASTools) # For mendelian error checking
library(tidyverse)
library(vcfR)
library(readxl)
library(pedigree)
library(pedigreeTools)
library(kinship2)

proj_dir <- here::here()

# Cranberry dir
cran_dir <- strsplit(proj_dir, "/")[[1]] %>%
  {.[seq_len(which(. == "CranberryLab"))]} %>%
  paste0(collapse = "/")

# Set the keyfile name
keyfile <- file.path(proj_dir, "input/cranberry_gbs_unique_keys_resolved_duplicates_notForTassel.txt")
# Population metadata
population_metadata <- read_excel(path = file.path(cran_dir, "Breeding/Populations/population_inventory.xlsx"))

# Filepath of unphased genotype database
unphased_geno_db_file <- file.path(cran_dir, "Genotyping/MarkerDatabase/unphased_marker_genotype_db")

# Filepath of phased genotype database
phased_geno_db_file <- file.path(cran_dir, "Genotyping/MarkerDatabase/phased_marker_genotype_db")


# Read in the keyfile; create renaming dataframe --------------------------

sample_keys <- read_tsv(file = keyfile)

# Distinct sample names
distinct_sample_names <- sample_renames <- sample_keys %>% 
  distinct(FullSampleName, GenotypeName)

sample_renames <- setNames(object = distinct_sample_names$GenotypeName, 
                           nm = distinct_sample_names$FullSampleName)



# Read in the unphased GBS data -------------------------------------------

vcf_filename <- file.path(proj_dir, "snps/cranberryGBS_production_snps_filtered.vcf.gz")
unphased_vcf_in <- read.vcfR(file = vcf_filename)


# Rename columns
colmatch <- match(x = colnames(unphased_vcf_in@gt)[-1], table = names(sample_renames))
colnames(unphased_vcf_in@gt)[-1] <- sample_renames[colmatch]

# Edit sample names
sample_names_use <- colnames(unphased_vcf_in@gt)[-1] %>% 
  neyhart::str_add_space() %>% 
  str_replace_all(" ", "_") %>% 
  str_to_upper() %>%
  ifelse(. == "POTTERS_FAVORITE", "POTTER", .)

# Sample with underscores that need replacement
which_underscores_replace <- str_detect(string = sample_names_use, pattern = "[0-9]{1,}_[0-9]{1,}")
sample_names_use[which_underscores_replace] <- str_replace_all(string = sample_names_use[which_underscores_replace], pattern = "_", replacement = "-")

colnames(unphased_vcf_in@gt)[-1] <- sample_names_use

vcf_in1 <- unphased_vcf_in

# Extract the gt as numeric
genos <- extract.gt(x = vcf_in1, element = "GT", as.numeric = FALSE, return.alleles = FALSE, convertNA = TRUE)

# Convert to 0, 1, 2
# Iterate over samples
gt3 <- apply(X = genos, MARGIN = 2, FUN = function(snp) {
  sapply(str_split(string = snp, pattern = "/"), function(x) sum(as.numeric(x)))
})

# Rename
dimnames(gt3) <- dimnames(genos)


## For any duplicated columns, choose the one with the lowest missing data
dup_colnames <- colnames(gt3)[duplicated(colnames(gt3))]

# First remove the duplicated name
gt4 <- gt3[,!colnames(gt3) %in% dup_colnames]
# Next select the lowest-missingness duplicate
dups_resolved <- lapply(X = dup_colnames, FUN = function(dup_name) {
  snps <- gt3[,colnames(gt3) %in% dup_name,drop = FALSE]
  snps[,which.min(colMeans(is.na(snps))),drop = FALSE]
})

# Merge with gt4
gt5 <- cbind(gt4, do.call("cbind", dups_resolved))

# Get the map information
map_df <- vcf_in1@fix %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename_all(tolower) %>%
  select(id, chrom, pos) %>%
  mutate_at(vars(chrom, pos), parse_number) %>%
  as.data.frame() %>%
  # Reorder - third column and more are ignored by the program
  select(chrom, pos, marker = id)




# Create pedigree information for the populations -------------------------

mp_pop_meta <- population_metadata %>%
  filter(population_name == "MP")

# Iterate over the materal/paternal parents and find their parents and their parents
# and so on
parents_i <- parents <- union(unique(mp_pop_meta$maternal_parent), unique(mp_pop_meta$paternal_parent))
all_parents_found <- FALSE

while(!all_parents_found) {
  parent_pedigree1 <- filter(population_metadata, individual %in% parents)
  parents_i <- union(parents, union(unique(parent_pedigree1$maternal_parent), unique(parent_pedigree1$paternal_parent)))
  all_parents_found <- all(parents_i %in% parents)
  parents <- parents_i
}

# Format pedigree
pedigree_df <- bind_rows(parent_pedigree1, mp_pop_meta) %>% 
  select(individual, maternal_parent, paternal_parent)

# Create a sorted pedigree
ordered_ped <- pedigreeTools::editPed(sire = pedigree_df$paternal_parent, dam = pedigree_df$maternal_parent,
                                      label = pedigree_df$individual)

## Save this pedigree in AlphaFamImpute format
afi_ped <- ordered_ped %>% 
  rownames_to_column("row") %>% 
  select(individual = label, sire, dam) %>% 
  mutate_at(vars(sire, dam), ~ifelse(is.na(.), "0", .))

# Save
write_delim(x = afi_ped, path = file.path(proj_dir, "imputation/alphaFamImpute/cranberry_gbs_mp_pedigree.txt"), 
            delim = " ", col_names = FALSE)




# Test genotypes for Mendenlian errors -------------------------------------

geno_mat <- t(gt5)

## Make SEARLES all NA (the assigned Searles is probably incorrect)
geno_mat["SEARLES",] <- NA

# Convert to data.frame; sort by the pedigree
geno_mat1 <- as.data.frame(geno_mat) %>%
  rownames_to_column("individual") %>%
  mutate(individual = case_when(
    individual == "US88-8" ~ "#35",
    individual == "US88-72" ~ "MCFARLIN",
    TRUE ~ individual)) %>%
  # Filter out some individuals
  subset(individual %in% ordered_ped$label) %>%
  mutate(individual = factor(individual, levels = ordered_ped$label)) %>%
  arrange(individual) %>%
  mutate(individual = as.character(individual)) %>%
  column_to_rownames("individual") %>%
  as.matrix()

# Create integer vector of SNP id
map_df1 <- mutate(map_df, snpID = seq_len(nrow(map_df)))

# Create integer vector of sample ID
sample_df <- tibble(individual = row.names(geno_mat1)) %>%
  mutate(scanID = seq_len(nrow(.)))


# Create a GenotypeData object
mgr <- MatrixGenotypeReader(genotype = t(geno_mat1), snpID = map_df1$snpID, chromosome = as.integer(map_df1$chrom), 
                            position = as.integer(map_df1$pos), scanID = sample_df$scanID)
genoData <- GenotypeData(data = mgr, scanAnnot = ScanAnnotationDataFrame(data = data.frame(scanID = sample_df$scanID, sex = "M")))

# Split the pedigrees by nuclear families; create lists with each extended family
afi_ped_nuclear_list <- afi_ped %>%
  filter_at(vars(sire, dam), all_vars(. != "0")) %>%
  group_by(sire, dam) %>%
  nest() %>%
  ungroup() %>%
  mutate(fam_id = paste0("fam", seq_len(nrow(.))),
         fam_id = as.numeric(as.factor(fam_id))) %>%
  split(.$fam_id) %>%
  map(~{
    df <- .
    df1 <- unnest(df, cols = c(data))
    
    # Get all ancestors
    ped_i <- distinct(df1, sire, dam)
    pars <- as.character(unlist(ped_i))
    while (any(pars != "0")) {
      ped_i <- subset(afi_ped, individual %in% pars)
      df1 <- bind_rows(ped_i, df1)
      pars <- as.character(unlist(distinct(ped_i, sire, dam)))
    }
    
    # Fill in the fam_id and scan ID
    df1 %>%
      mutate(fam_id = unique(na.omit(fam_id)),
             scanID = sample_df$scanID[match(x = individual, table = sample_df$individual)])
    
  })
  

# Empty list to store errors
mendel_errors_out <- list()
# Copy the genotype matrix
geno_mat2 <- geno_mat1

# Iterate over nuclear families
for (i in seq_along(afi_ped_nuclear_list)) {
  # Subset the pedigree
  afi_ped_nuclear <- afi_ped_nuclear_list[[i]]

  # Create mendel lists
  mendel_list <- mendelList(familyid = afi_ped_nuclear$fam_id, offspring = afi_ped_nuclear$individual, 
                            father = afi_ped_nuclear$sire, mother = afi_ped_nuclear$dam, sex = ifelse(afi_ped_nuclear$individual %in% afi_ped_nuclear$sire, "M", "F"),
                            scanID = afi_ped_nuclear$scanID)
  
  # Screen for mendel errors
  mendel_err_out <- mendelErr(genoData = genoData, mendel.list = mendel_list, error.by.snp = TRUE, 
                              error.by.snp.trio = TRUE, verbose = FALSE)
  
  # Set mendelian errors to NA
  mendel_err_snp <- mendel_err_out$snp
  indiv <- tibble(listname = names(mendel_err_snp)) %>%
    filter(str_detect(listname, paste0("^", i))) %>%
    mutate(individual = str_remove_all(listname, paste0("^", i, "\\.")),
           snps_NA = map(listname, ~which(mendel_err_snp[[.x]] != 0)))
  
  # Iterate over individuals and set missing
  for (j in seq_len(nrow(indiv))) {
    snp_i <- indiv$snps_NA[[j]]
    geno_mat2[indiv$individual[j], snp_i] <- NA
    
  }
  
  # Add the summary to the list
  mendel_errors_out[[i]] <- mendel_err_out$all.trios
  
}


# Filter SNPs
# Remove singleton alleles
# Filter for less than 80% missingness
# 
geno_mat3 <- snps::filter_snps(x = geno_mat2 - 1, r2.max = 1, maf.min = 1 / nrow(geno_mat2), indiv.miss.max = 0.80, 
                               snp.miss.max = 0.80)

# Add SEARLES back in as completely missing
geno_mat4 <- rbind(geno_mat3, SEARLES = NA) + 1



# Convert the data for AlphaPeel ------------------------------

geno_mat4[is.na(geno_mat4)] <- 9

# Convert back to a data.frame; sort by the pedigree
geno_df <- as.data.frame(geno_mat4) %>%
  rownames_to_column("individual")

# Save
write_delim(x = geno_df, path = file.path(proj_dir, "imputation/alphaFamImpute/cranberry_gbs_mp_genotypes.txt"),
            delim = " ", col_names = FALSE)


# Save a map file ---------------------------------------------------------

# Filter the map for SNPs in the final matrix
map_df1 <- map_df %>%
  filter(marker %in% names(geno_df))

# Save
write_delim(x = map_df1, path = file.path(proj_dir, "imputation/alphaFamImpute/cranberry_gbs_mp_map.txt"),
            delim = " ", col_names = FALSE)




# Run AlphaFamImpute using python -----------------------------------------

# 1. Open Anaconda prompt
# 2. cd C:\Users\jeffrey.neyhart\OneDrive - USDA\Documents\CranberryLab\Genotyping\pre2021\cranberryHistoricalGBS\imputation\alphaFamImpute
# 3. AlphaFamImpute -gbs -genotypes cranberry_gbs_mp_genotypes.txt -map cranberry_gbs_mp_map.txt -pedigree cranberry_gbs_mp_pedigree.txt -out cranberry_gbs_mp




# Create a spec file for AlphaPeel ----------------------------------------

## Edit the pedigree to replace individual names with numbers
pedigree_input <- file.path(proj_dir, "imputation/alphaPeel/cranberry_gbs_mp_pedigree.txt")
pedigree_df <- read_delim(file = file.path(proj_dir, "imputation/alphaFamImpute/cranberry_gbs_mp_pedigree.txt"), 
                          delim = " ", col_names = c("individual", "sire", "dam"))
# Create a key
indiv_key <- c("0" = "0", set_names(seq_along(pedigree_df$individual), pedigree_df$individual))

# Replace and save
pedigree_df %>%
  # mutate_all(~map_chr(., ~indiv_key[.])) %>%
  write_delim(x = ., path = pedigree_input, delim = " ", col_names = FALSE)

# Edit files for alphapeel
# Split the genotypes by chromosome
markers_by_chrom <- split(map_df1$marker, map_df1$chrom)

for (i in seq_along(markers_by_chrom)) {
  chrom <- names(markers_by_chrom)[i]
  geno_df_i <- geno_df[,c("individual", markers_by_chrom[[i]])] #%>% mutate(individual = map_chr(individual, ~indiv_key[.]))
    
  # Save the geno file
  geno_input_i <- file.path(proj_dir, "imputation/alphaPeel",paste0("cranberry_gbs_mp_genotypes_chrom", chrom, ".txt"))
  write_delim(x = geno_df_i, path = geno_input_i, delim = " ", col_names = FALSE)
  
  # Save the spec file
  alphaPeel_specs_df <- tribble(
    ~X1, ~X2,
    "nsnp", as.character(ncol(geno_df_i) - 1),
    "inputfilepath", basename(geno_input_i),
    "pedigree", basename(pedigree_input),
    "outputfilepath", paste0("cranberry_gbs_mp_alphaPeelOutput_chrom", chrom),
    "runtype", "multi",
    "nCycles", "10",
  )
  
  write_delim(x = alphaPeel_specs_df, delim = ", ", col_names = FALSE,
              path = file.path(proj_dir, "imputation/alphaPeel",paste0("cranberry_gbs_mp_AlphaPeelSpec_chrom", chrom, ".txt")))
  
}



# Run AlphaPeel using linux -----------------------------------------

# 1. Open Ubuntu
# 2. cran
# 3. cd Genotyping/pre2021/cranberryHistoricalGBS/imputation/alphaPeel/
# 4. Run this:
# 
# for specfile in $(find . -name "*AlphaPeelSpec*"); do ./AlphaPeel $specfile; done
# 





