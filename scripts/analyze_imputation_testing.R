## cranberryHistoricalGBS
## 
## Analyze the results of imputation
## 


# Directories and packages
library(tidyverse)
library(vcfR)

proj_dir <- here::here()
pipeline_dir <- proj_dir
impute_dir <- file.path(proj_dir, "imputation/beagle_testing/")
impute_test_missing_dir <- file.path(proj_dir, "imputation/beagle_testing/random_missing_files/")
impute_test_dir <- file.path(proj_dir, "imputation/beagle_testing/imputed_files/")


# Read in the VCF ---------------------------------------------------------

# Read in
filename <- file.path(impute_dir, "snps_for_imputation_indv_filtered.recode.vcf")
vcf_in <- read.vcfR(file = filename)


# Get the genotype data
gt <- vcf_in@gt
gt1 <- gt[,-1,drop = FALSE]

# Extract the genotype calls; calculate number of reference alleles
gt2 <- apply(X = gt1, MARGIN = 2, FUN = function(indv) {
  indv_gt <- sapply(str_split(string = indv, pattern = ":"), "[[", 1)
  2 - sapply(X = lapply(X = str_split(string = indv_gt, pattern = "/"), as.numeric), sum)
})

# Remove any genotypes from unknown chromosomes
gt3 <- gt2[vcf_in@fix[,"CHROM"] != "UNKNOWN",,drop = FALSE]
# Index of non-missing for this file
gt3_non_missing <- !is.na(gt3)



# Compare results of the whole sample -------------------------------------

# List the files
subdir_name <- "all_individuals"
missing_files <- list.files(path = file.path(impute_test_missing_dir, subdir_name),
                            pattern = ".vcf.gz", full.names = TRUE)

# Create an empty list
missing_genotype_calls_list <- list()

# Read in each
for (i in seq_along(missing_files)) {

  vcf_missing_i <- read.vcfR(file = missing_files[i], verbose = FALSE)
  
  # Extract the genotype calls; convert to number of reference alleles
  gt_missing_i <- vcf_missing_i@gt[vcf_missing_i@fix[,"CHROM"] != "UNKNOWN",-1, drop = FALSE]
  
  gt_missing_i2 <- apply(X = gt_missing_i, MARGIN = 2, FUN = function(indv) {
    indv_gt <- sapply(str_split(string = indv, pattern = ":"), "[[", 1)
    2 - sapply(X = lapply(X = str_split(string = indv_gt, pattern = "/"), as.numeric), sum)
  })

  # Create a matrix pointing to the missing data
  missing_location_mat <- is.na(gt_missing_i2)
  
  # Add this to the list; use the basename of the file
  name_i <- gsub(pattern = ".vcf.gz", replacement = "", x = basename(missing_files[i]))
  missing_genotype_calls_list[[name_i]] <- missing_location_mat
  
}




# List the imputed files
imputed_files <- list.files(path = impute_test_dir, pattern = ".vcf.gz", full.names = TRUE)


# Calculate imputation accuracy

# Create an empty list
accuracy_list <- list()

# Read in each file, compare with the original, calculate the error rate
for (i in seq_along(imputed_files)) {
  file_i <- imputed_files[i]
  
  ## Read in the imputed VCF ##
  vcf_in_i <- read.vcfR(file = file_i, verbose = FALSE)
  
  # Extract the genotype calls
  # Remove unknown chromosomes
  gt_i <- vcf_in_i@gt[vcf_in_i@fix[,"CHROM"] != "UNKNOWN",-1, drop = FALSE]
  # Replace "|" with "/"
  gt_convt <- gsub(pattern = "\\|", replacement = "/", x = gt_i)
  # calculate number of reference alleles
  gt_convt1 <- apply(X = gt_convt, MARGIN = 2, FUN = function(indv) {
    2 - sapply(X = lapply(X = str_split(string = indv, pattern = "/"), as.numeric), sum)
  })
  
  # Find the corresponding missing file name
  base_filename <- str_extract(string = basename(file_i), pattern = "test_beagle_missing_[0-9]{1,}")
  missing_index_i <- missing_genotype_calls_list[[base_filename]]
   
  # Compare all original genotypes with imputed
  gt_compare_all <- gt3 == gt_convt1
  
  # Compare only missing original genotypes with imputed
  gt_compare_missing <- gt3[missing_index_i & gt3_non_missing] == gt_convt1[missing_index_i & gt3_non_missing]
  
  # Extract metadata from the file
  ne <- parse_number(str_extract(string = basename(file_i), pattern = "Ne[0-9]{1,}"))
  
  # Return both accuracies
  accuracy_list[[i]] <- tibble(file = base_filename,
                               Ne = ne,
                               accuracy_compare_all = mean(gt_compare_all, na.rm = TRUE), 
                               accuracy_compare_missing = mean(gt_compare_missing, na.rm = TRUE))

}


# Bind the list
accuracy_df_all_indiv <- bind_rows(accuracy_list)

accuracy_df %>%
  group_by(Ne) %>% 
  summarize_at(vars(contains("accuracy")), mean) %>%
  as.data.frame()




# Compare results of just the germplasm collection -------------------------------------

# List the files
subdir_name <- "germplasm_collection"
missing_files <- list.files(path = file.path(impute_test_missing_dir, subdir_name),
                            pattern = ".vcf.gz", full.names = TRUE)

# Create an empty list
missing_genotype_calls_list <- list()

# Read in each
for (i in seq_along(missing_files)) {
  
  vcf_missing_i <- read.vcfR(file = missing_files[i], verbose = FALSE)
  
  # Extract the genotype calls; convert to number of reference alleles
  gt_missing_i <- vcf_missing_i@gt[vcf_missing_i@fix[,"CHROM"] != "UNKNOWN",-1, drop = FALSE]
  
  gt_missing_i2 <- apply(X = gt_missing_i, MARGIN = 2, FUN = function(indv) {
    indv_gt <- sapply(str_split(string = indv, pattern = ":"), "[[", 1)
    2 - sapply(X = lapply(X = str_split(string = indv_gt, pattern = "/"), as.numeric), sum)
  })
  
  # Create a matrix pointing to the missing data
  missing_location_mat <- is.na(gt_missing_i2)
  
  # Add this to the list; use the basename of the file
  name_i <- gsub(pattern = ".vcf.gz", replacement = "", x = basename(missing_files[i]))
  missing_genotype_calls_list[[name_i]] <- missing_location_mat
  
}




# List the imputed files
subdir_name <- "germplasm_collection"
imputed_files <- list.files(path = file.path(impute_test_dir, subdir_name),
                            pattern = ".vcf.gz", full.names = TRUE)
 

# Calculate imputation accuracy

# Create an empty list
accuracy_list <- list()

# Read in each file, compare with the original, calculate the error rate
for (i in seq_along(imputed_files)) {
  file_i <- imputed_files[i]
  
  ## Read in the imputed VCF ##
  vcf_in_i <- read.vcfR(file = file_i, verbose = FALSE)
  
  # Extract the genotype calls
  # Remove unknown chromosomes
  gt_i <- vcf_in_i@gt[vcf_in_i@fix[,"CHROM"] != "UNKNOWN",-1, drop = FALSE]
  # Replace "|" with "/"
  gt_convt <- gsub(pattern = "\\|", replacement = "/", x = gt_i)
  # calculate number of reference alleles
  gt_convt1 <- apply(X = gt_convt, MARGIN = 2, FUN = function(indv) {
    2 - sapply(X = lapply(X = str_split(string = indv, pattern = "/"), as.numeric), sum)
  })
  
  # Vector of sample names
  gt_sample_names <- colnames(gt_convt1)
  
  # Find the corresponding missing file name
  base_filename <- str_extract(string = basename(file_i), pattern = "test_beagle_missing_[0-9]{1,}")
  missing_index_i <- missing_genotype_calls_list[[base_filename]]
  missing_index_i2 <- missing_index_i & gt3_non_missing[, gt_sample_names, drop = FALSE]
  
  # Compare all original genotypes with imputed
  gt_compare_all <- gt3[,gt_sample_names,drop = FALSE] == gt_convt1
  
  # Compare only missing original genotypes with imputed
  gt_compare_missing <- gt3[,gt_sample_names][missing_index_i2] == gt_convt1[missing_index_i2]
  
  # Extract metadata from the file
  ne <- parse_number(str_extract(string = basename(file_i), pattern = "Ne[0-9]{1,}"))
  
  # Return both accuracies
  accuracy_list[[i]] <- tibble(file = base_filename,
                               Ne = ne,
                               accuracy_compare_all = mean(gt_compare_all, na.rm = TRUE), 
                               accuracy_compare_missing = mean(gt_compare_missing, na.rm = TRUE))
  
}


# Bind the list
accuracy_df_germcol <- bind_rows(accuracy_list)

accuracy_df_germcol %>%
  group_by(Ne) %>% 
  summarize_at(vars(contains("accuracy")), mean) %>%
  as.data.frame()

