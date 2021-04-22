## cranberryHistoricalGBS
## 
## Save SNPs to a common SNP database (housed as a VCF file)
## 
## 


# Packages
library(tidyverse)
library(vcfR)

proj_dir <- here::here()

# Cranberry dir
cran_dir <- strsplit(proj_dir, "/")[[1]] %>%
  {.[seq_len(which(. == "CranberryLab"))]} %>%
  paste0(collapse = "/")

# Set the keyfile name
keyfile <- file.path(proj_dir, "input/cranberry_gbs_unique_keys_resolved_duplicates.txt")

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



# Organize the unphased VCF -----------------------------------------------

vcf_filename <- file.path(proj_dir, "snps/cranberryGBS_production_snps_filtered.vcf.gz")

unphased_vcf_in <- read.vcfR(file = vcf_filename)


# Rename columns
colmatch <- match(x = colnames(unphased_vcf_in@gt)[-1], table = names(sample_renames))
colnames(unphased_vcf_in@gt)[-1] <- sample_renames[colmatch]

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

# Gather SNP metadata
snp_info <- vcf_in1@fix %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate_if(is.character, parse_guess) %>%
  rename_all(tolower) %>%
  unite(alleles, ref, alt, sep = "/") %>%
  select(marker = id, chrom, pos, alleles)

# Combine into hapmap version
geno_hmp <- cbind(snp_info, gt3) %>%
  remove_rownames() %>%
  as_tibble()


# Write a VCF
write.vcf(x = vcf_in1, file = paste0(unphased_geno_db_file, ".vcf"))
# Write a hmp file
write_tsv(x = geno_hmp, path = paste0(unphased_geno_db_file, "_hmp.txt"))



# Organize the phased VCF ---------------------------------------------------

# Read in the VCF file
filename <- file.path(proj_dir, "/imputation/beagle_imputation/cranberryGBS_germplasm_imputed_snps.vcf.gz")
vcf_in <- read.vcfR(file = filename)

# Rename columns
colmatch <- match(x = colnames(vcf_in@gt)[-1], table = names(sample_renames))
colnames(vcf_in@gt)[-1] <- sample_renames[colmatch]

vcf_in1 <- vcf_in

# Extract the gt as numeric
genos <- extract.gt(x = vcf_in1, element = "GT", as.numeric = FALSE, return.alleles = FALSE, convertNA = TRUE)

# Create a haplotype array (2 x m x n)
haplo_array <- array(NA, dim = c(2, dim(genos)), dimnames = list(NULL, rownames(genos), colnames(genos)))

# Iterate over samples
for (i in seq_len(ncol(genos))) {
  entry <- genos[,i]
  haplo_array[,,i] <- t(apply(X = do.call("rbind", str_split(entry, "\\|")), MARGIN = 2, FUN = as.numeric))
}


# Create a genotype matrix
# Sum the number of reference alleles; subtract 1
geno_mat <- t(apply(X = haplo_array, MARGIN = c(2,3), sum))

# Gather SNP metadata
snp_info <- vcf_in1@fix %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate_if(is.character, parse_guess) %>%
  rename_all(tolower) %>%
  unite(alleles, ref, alt, sep = "/") %>%
  select(marker = id, chrom, pos, alleles)

# Combine into hapmap version
geno_hmp <- cbind(snp_info, t(geno_mat)) %>%
  remove_rownames() %>%
  as_tibble()


# Write a VCF
write.vcf(x = vcf_in1, file = paste0(phased_geno_db_file, ".vcf"))
# Write a hmp file
write_tsv(x = geno_hmp, path = paste0(phased_geno_db_file, "_hmp.txt"))

# Write an RData file
phased_geno_hmp <- geno_hmp
phased_geno_mat <- geno_mat
phased_geno_haplo_array <- haplo_array
save("phased_geno_hmp", "phased_geno_mat", "phased_geno_haplo_array",
     file = paste0(phased_geno_db_file, ".RData"))

