## cranberryHistoricalGBS
## 
## Generate random missing data samples
## 


# Directories and packages
library(tidyverse)
library(vcfR)

proj_dir <- here::here()
pipeline_dir <- proj_dir
impute_dir <- file.path(proj_dir, "imputation/beagle_testing/")
impute_test_dir <- file.path(proj_dir, "imputation/beagle_testing/random_missing_files/")

# Read in the VCF ---------------------------------------------------------

# Read in
filename <- file.path(impute_dir, "snps_for_imputation_indv_filtered.recode.vcf")
vcf_in <- read.vcfR(file = filename)


# Get the genotype data
gt <- vcf_in@gt
gt1 <- gt[,-1,drop = FALSE]

# Create a randomization matrix for assigning missing data
# This proportion is based on the missing data amount observed in the data
pMiss <- 0.15

# Number of replications
nRep <- 10

# Iterate
for (i in seq_len(nRep)) {
  
  # Apply over rows
  gt_miss <- apply(X = gt, MARGIN = 1, FUN = function(row) {
    # Pull genotypes
    geno <- row[-1]
    
    # Randomly assign missing
    missing_i <- as.logical(rbinom(n = length(geno), size = 1, prob = pMiss))
    
    # Subset these calls; split each string
    calls_missing_i <- str_split(string = geno[missing_i], pattern = ":")
    # Convert the first element to missing
    calls_missing_i1 <- map(calls_missing_i, ~modify_at(.x, .at = 1, .f = ~"./.")) %>%
      map_chr(~paste0(., collapse = ":"))
    
    # Replace these
    geno[missing_i] <- calls_missing_i1
    
    # Return the row
    row[-1] <- geno
    return(row)
    
  })
  
  # Transpose
  gt_miss1 <- t(gt_miss)
  
  # Copy the VCF
  vcf1 <- vcf_in
  
  # Replace the gt
  vcf1@gt <- gt_miss1
  
  # Write the file
  new_filename <- file.path(impute_test_dir, 
                            paste0("test_beagle_missing_", formatC(x = i, width = nchar(nRep), flag = "0"), ".vcf.gz"))
  
  write.vcf(x = vcf1, file = new_filename)
  
}










