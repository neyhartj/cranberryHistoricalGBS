## cranberryHistoricalGBS
## 
## Analyze the results of imputation
## 


# Directories and packages
library(igraph)
library(tidyverse)
library(readxl)
library(vcfR)
library(neyhart)


proj_dir <- here::here()
pipeline_dir <- proj_dir
snps_dir <-  file.path(proj_dir, "/snps/")
fig_dir <- file.path(proj_dir, "/figures/")

# Cranberry dir
cran_dir <- find_dir(pattern = "CranberryLab")



# Analyze depth -----------------------------------------------------------

# Read in the depth data
mean_indiv_depth <- read_tsv(file = file.path(snps_dir, "cranberryGBS_production_snps_resolvedDuplicates_filtered_meanIndivDepth.idepth")) %>%
  rename_all(tolower)
mean_site_depth <- read_tsv(file = file.path(snps_dir, "cranberryGBS_production_snps_resolvedDuplicates_filtered_meanSiteDepth.ldepth.mean")) %>%
  rename_all(tolower)
geno_depth <- read_tsv(file = file.path(snps_dir, "cranberryGBS_production_snps_resolvedDuplicates_filtered_genoDepth.gdepth"))

# Convert geno_depth to matrix
geno_depth_mat <- geno_depth %>%
  unite("site", c("CHROM", "POS")) %>%
  as.data.frame() %>%
  column_to_rownames("site") %>%
  as.matrix()


# Histogram of average individual depth
hist(mean_indiv_depth$mean_depth, breaks = 30, main = "Mean Individual Depth", xlab = "Depth")

# Histogram of average site depth
hist(mean_site_depth$mean_depth, breaks = 30, main = "Mean Site Depth", xlab = "Depth")

### 
### Potential maximum genotype depths
### 

# Twice the average depth plus 1
maxDP <- round((mean(geno_depth_mat) * 2) + 1)

# How many genotype calls would this erase?
geno_depth_mat_filter1 <- geno_depth_mat
geno_depth_mat_filter1[geno_depth_mat > maxDP] <- NA

mean(is.na(geno_depth_mat_filter1))
hist(colMeans(is.na(geno_depth_mat_filter1)), main = "Mean Individual Missingness")
hist(rowMeans(is.na(geno_depth_mat_filter1)), main = "Mean Site Missingness")









