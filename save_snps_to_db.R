## cranberryHistoricalGBS
## 
## Save SNPs to a common SNP database (housed as a VCF file)
## 
## 


# Packages
library(tidyverse)
library(vcfR)
library(snps)

proj_dir <- here::here()

# Cranberry dir
cran_dir <- strsplit(proj_dir, "/")[[1]] %>%
  {.[seq_len(which(. == "CranberryLab"))]} %>%
  paste0(collapse = "/")

# Read in the recombination rate map
cMMb <- read_csv(file = file.path(cran_dir, "Genotyping/GeneticMaps/genetic_physical_rate_map.csv"))


# Output directory
marker_dir <- file.path(cran_dir, "Genotyping/MarkerDatabase/")

# Set the keyfile name
keyfile <- file.path(proj_dir, "input/cranberry_gbs_unique_keys_resolved_duplicates_notForTassel.txt")

# Filepath of unphased genotype database
unphased_geno_db_file <- file.path(cran_dir, "Genotyping/MarkerDatabase/unphased_marker_genotype_db")
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

# Edit sample names
sample_names_use <- colnames(genos) %>% 
  neyhart::str_add_space() %>% 
  str_replace_all(" ", "_") %>% 
  str_to_upper() %>%
  ifelse(. == "POTTERS_FAVORITE", "POTTER", .)

colnames(genos) <- sample_names_use


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


# Gather SNP metadata
snp_info <- vcf_in1@fix %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate_if(is.character, parse_guess) %>%
  rename_all(tolower) %>%
  unite(alleles, ref, alt, sep = "/") %>%
  select(marker = id, chrom, pos, alleles) %>%
  mutate(chrom = parse_number(chrom)) %>%
  # Predict genetic positions
  interp_gen_position(map = ., cMMb = cMMb) %>%
  mutate(chrom = str_pad(chrom, 2, pad = "0"))

# Combine into hapmap version
geno_hmp <- cbind(snp_info, gt5) %>%
  remove_rownames() %>%
  as_tibble()


# Write a VCF
write.vcf(x = vcf_in1, file = paste0(unphased_geno_db_file, ".vcf"))
# Write a hmp file
write_tsv(x = geno_hmp, path = paste0(unphased_geno_db_file, "_hmp.txt"))



# Organize the phased genotypes ---------------------------------------------------


# Beagle ------------------------------------------------------------------



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


## For any duplicated columns, choose the first one (all will have the same
## level of missingness; 0%)
#
# Remove the duplicated name
haplo_array1 <- haplo_array[,,!duplicated(dimnames(haplo_array)[[3]])]

# Gather SNP metadata
snp_info_beagle <- snp_info %>%
  filter(marker %in% colnames(haplo_array1))


# Save the haplotype array
phased_geno_haplo_array <- haplo_array1
save("phased_geno_haplo_array", "snp_info_beagle", file = file.path(marker_dir, "beagle_phased_marker_genotypes.RData"))


# AlphaFamImpute ----------------------------------------------------------


## Read in results from AlphaFamImpute
afi_genos <- read_delim(file = file.path(proj_dir, "imputation/alphaFamImpute/cranberry_gbs_mp.genotypes"), 
                        delim = " ", col_names = FALSE)
afi_phase <- read_delim(file = file.path(proj_dir, "imputation/alphaFamImpute/cranberry_gbs_mp.phase"), 
                        delim = " ", col_names = FALSE)
afi_dosage <- read_delim(file = file.path(proj_dir, "imputation/alphaFamImpute/cranberry_gbs_mp.dosages"),
                         delim = " ", col_names = FALSE)
# Read in the map used for this imputation
afi_map_in <- read_delim(file = file.path(proj_dir, "imputation/alphaFamImpute/cranberry_gbs_mp_map.txt"),
                         delim = " ", col_names = c("chrom", "pos", "marker"))

# Genotype matrix
afi_genos_mat <- afi_genos %>%
  as.data.frame() %>%
  column_to_rownames("X1") %>%
  as.matrix() %>%
  `colnames<-`(., afi_map_in$marker)

# Dosage matrix
afi_dosage_mat <- afi_dosage %>%
  as.data.frame() %>%
  column_to_rownames("X1") %>%
  as.matrix() %>%
  `colnames<-`(., afi_map_in$marker)

## Create another haplotype array for the AlphaFamImpute output
afi_phase_mat <- `names<-`(afi_phase, c("individual", afi_map_in$marker)) %>%
  mutate(individual = paste0(individual, c("_mat", "_pat"))) %>% # add _mat and _pat to haplotype name
  column_to_rownames("individual") %>%
  as.matrix()

afi_phase_mat_split <- asplit(x = afi_phase_mat, MARGIN = 1)
afi_haplo_array <- sapply(X = seq(1, length(afi_phase_mat_split), by = 2), 
                          FUN = function(i) rbind(afi_phase_mat_split[[i]], afi_phase_mat_split[[i+1]]), simplify = FALSE) %>%
  unlist() %>%
  array(data = ., dim = c(2, length(afi_phase_mat_split[[1]]), length(afi_phase_mat_split) / 2),
        dimnames = list(c("maternal", "paternal"), names(afi_phase_mat_split[[1]]), unique(gsub(pattern = "_mat|_pat", replacement = "", x = row.names(afi_phase_mat)))))

# Gather SNP metadata
snp_info_afi <- snp_info %>%
  filter(marker %in% colnames(afi_haplo_array))



# Save the haplotype array
save("afi_haplo_array", "afi_dosage_mat", "afi_genos_mat", "snp_info_afi",
     file = file.path(marker_dir, "alphaFamImpute_phased_marker_genotypes.RData"))


# AlphaPeel ---------------------------------------------------------------

# List the output files
all_alphapeel_output_files <- list.files(path = file.path(proj_dir, "imputation/alphaPeel/"),
                                         pattern = "Output", full.names = TRUE)

# Combine the dosage files
ap_dosage_mat <- str_subset(string = all_alphapeel_output_files, pattern = ".dosages") %>%
  map(read_delim, delim = " ", col_names = FALSE) %>%
  map(as.data.frame) %>%
  map(column_to_rownames, "X1") %>%
  bind_cols() %>%
  `names<-`(afi_map_in$marker) %>%
  as.matrix()

ap_dosage_mat1 <- ap_dosage_mat %>%
  str_trim() %>% 
  parse_number() %>%
  matrix(data = ., nrow = nrow(ap_dosage_mat), ncol = ncol(ap_dosage_mat), dimnames = dimnames(ap_dosage_mat))

# Round the dosage matrix to get a genotype matrix
ap_genos_mat <- round(ap_dosage_mat1)

# Combine the haps files
ap_haps_mat <- str_subset(string = all_alphapeel_output_files, pattern = ".haps") %>%
  map(read_delim, delim = " ", col_names = FALSE) %>%
  map(~mutate(., X1 = paste0(X1, c("_aa", "_aA", "_Aa", "_AA")))) %>%
  map(as.data.frame) %>%
  map(column_to_rownames, "X1") %>%
  bind_cols() %>%
  `names<-`(afi_map_in$marker) %>%
  as.matrix()

ap_haps_mat1 <- ap_haps_mat %>%
  str_trim() %>% 
  parse_number() %>%
  matrix(data = ., nrow = nrow(ap_haps_mat), ncol = ncol(ap_haps_mat), dimnames = dimnames(ap_haps_mat))

# Combine the segregation files
ap_segs_mat <- str_subset(string = all_alphapeel_output_files, pattern = ".seg") %>%
  map(read_delim, delim = " ", col_names = FALSE) %>%
  map(~mutate(., X1 = paste0(X1, c("_pafa_pama", "_pafa_mama", "_mafa_pama", "_mafa_mama")))) %>%
  map(as.data.frame) %>%
  map(column_to_rownames, "X1") %>%
  bind_cols() %>%
  `names<-`(afi_map_in$marker) %>%
  as.matrix()

ap_segs_mat1 <- ap_segs_mat %>%
  str_trim() %>% 
  parse_number() %>%
  matrix(data = ., nrow = nrow(ap_segs_mat), ncol = ncol(ap_segs_mat), dimnames = dimnames(ap_segs_mat))


# ## Plot #35 inheritance probability for the CNJ02 population
# rows <- intersect(str_which(string = row.names(ap_segs_mat1), pattern = "CNJ02"),
#                   str_which(string = row.names(ap_segs_mat1), pattern = "mafa")) # mafa = maternal grandparent from father
# 
# no35_rows <- ap_segs_mat1[rows,]
# no35_rows1 <- sapply(X = seq(1, length(rows), by = 2), FUN = function(i) colSums(no35_rows[c(i, i + 1),, drop = FALSE]),
#                      simplify = FALSE)
# 
# no35_rows2 <- do.call("rbind", no35_rows1)
# row.names(no35_rows2) <- unique(str_remove_all(string = row.names(no35_rows), pattern = "_mafa_mama|_mafa_pama"))
# 
# # Visualize
# image(t(no35_rows2[,subset(afi_map_in, chrom == 1, marker, drop = T)]))



# Rename
ap_haps_mat <- ap_haps_mat1
ap_dosage_mat <- ap_dosage_mat1
ap_segs_mat <- ap_segs_mat1

# Gather SNP metadata
snp_info_ap <- snp_info %>%
  filter(marker %in% colnames(ap_genos_mat))


# Save the results
save("ap_haps_mat", "ap_dosage_mat", "ap_segs_mat", "ap_genos_mat", "snp_info_ap",
     file = file.path(marker_dir, "alphaPeel_phased_marker_genotypes.RData"))



# Combine genotype matrices ------------------------------------------------

# Edit the sample names in phased_geno_haplo_array
sample_names <- dimnames(phased_geno_haplo_array)[[3]]
sample_names1 <- sample_names %>%
  neyhart::str_add_space() %>% 
  str_replace_all(" ", "_") %>% 
  str_to_upper() %>%
  ifelse(. == "POTTERS_FAVORITE", "POTTER", .)

dimnames(phased_geno_haplo_array)[[3]] <- sample_names1


## Merge the haplotype arrays
# Remove samples from phased_geno_haplo_array that are in afi_haplo_array
phased_geno_haplo_array1 <- phased_geno_haplo_array[,,!dimnames(phased_geno_haplo_array)[[3]] %in% dimnames(afi_haplo_array)[[3]]]
# Markers in afi_haplo_array that are missing in phased_geno_haplo_array get NA
missing1_snps <- setdiff(dimnames(afi_haplo_array)[[2]], dimnames(phased_geno_haplo_array1)[[2]])
missing1_snps_array <- array(data = NA, dim = c(nrow(phased_geno_haplo_array1), length(missing1_snps), dim(phased_geno_haplo_array1)[3]),
                             dimnames = list(NULL, missing1_snps, dimnames(phased_geno_haplo_array1)[[3]]))

phased_geno_haplo_array2 <- abind::abind(phased_geno_haplo_array1, missing1_snps_array, along = 2)
phased_geno_haplo_array2 <- phased_geno_haplo_array2[,sort(dimnames(phased_geno_haplo_array2)[[2]]),]

# Repeat the same for the afi_haplo_array
# Markers in phased_geno_haplo_array that are missing in afi_haplo_array get NA
missing2_snps <- setdiff(dimnames(phased_geno_haplo_array2)[[2]], dimnames(afi_haplo_array)[[2]])
missing2_snps_array <- array(data = NA, dim = c(nrow(afi_haplo_array), length(missing2_snps), dim(afi_haplo_array)[3]),
                             dimnames = list(NULL, missing2_snps, dimnames(afi_haplo_array)[[3]]))

afi_haplo_array2 <- abind::abind(afi_haplo_array, missing2_snps_array, along = 2)
afi_haplo_array2 <- afi_haplo_array2[,sort(dimnames(afi_haplo_array2)[[2]]),]

# Combine
phased_geno_haplo_array_merged <- abind::abind(phased_geno_haplo_array2, afi_haplo_array2, along = 3)



# Create a genotype matrix
# Sum the number of reference alleles; subtract 1
geno_mat <- t(apply(X = phased_geno_haplo_array_merged, MARGIN = c(2,3), sum))

snp_info_phased <- filter(snp_info, marker %in% colnames(geno_mat))

# Combine into hapmap version
geno_hmp <- cbind(snp_info_phased, t(geno_mat)) %>%
  remove_rownames() %>%
  as_tibble()


# Write a hmp file
write_tsv(x = geno_hmp, path = paste0(phased_geno_db_file, "_hmp.txt"))

# Write an RData file
phased_geno_hmp <- geno_hmp
phased_geno_mat <- geno_mat
phased_geno_haplo_array <- phased_geno_haplo_array_merged

save("phased_geno_hmp", "phased_geno_mat", "phased_geno_haplo_array", "snp_info_phased",
     file = paste0(file.path(marker_dir, "phased_marker_genotype_db.RData")))



# Build a VCF
phased_vcf1 <- unphased_vcf_in
phased_vcf1@fix <- phased_vcf1@fix[phased_vcf1@fix[,"ID"] %in% geno_hmp$marker,]
gt <- apply(X = phased_geno_haplo_array, MARGIN = c(2,3), paste0, collapse = "|")

phased_vcf1@gt <- cbind(FORMAT = "GT", gt)


# Write a VCF
write.vcf(x = phased_vcf1, file = file.path(marker_dir, "phased_marker_genotype_db.vcf"))

