#!/bin/bash

## 
## PLINK stats
## 
## This script will use PLINK to do the following:
## 1. generate haplotype blocks in the germplasm collection
##    a. do this for both wild and wild+native selections
## 2. 
## 
## 


# Load the modules
module load plink/0.1.9

## Set variables
# Working directory
WD=/project/cranberrygbs/cranberryHistoricalGBS/
# Name of input directory
INPUT=$WD/imputation/beagle_imputation/
# Output directory
OUTPUT=$WD/haplotyping

# Name of input vcf
VCFIN=$INPUT/cranberryGBS_germplasm_imputed_snps.vcf.gz


# Change working directory
cd $WD 


## Filter

# Filter to retain the germplasm collection, then filter SNPs for excessive 
# missingness (> 50%) and then filter individuals for excessive missingness (> 50%)
vcftools --gzvcf $VCFIN --keep $WD/input/cranberry_gbs_germplasm_collection_individuals.txt \
--recode-INFO-all --recode --stdout | \
vcftools --vcf - --mac 15 --max-missing 0.5 --recode-INFO-all --recode --out $OUTPUT/cranberryGBS_snps_to_impute

# Calculate individual missingness
vcftools --vcf $OUTPUT/cranberryGBS_snps_to_impute.recode.vcf --missing-indv \
--out $OUTPUT/cranberryGBS_germplasm_missingness

# Filter out individuals with missingness > 0.50
awk 'NR > 1 { if ($5 > 0.50) print $1 }' $OUTPUT/cranberryGBS_germplasm_missingness.imiss > $OUTPUT/individual_remove_imputation.txt

vcftools --vcf $OUTPUT/cranberryGBS_snps_to_impute.recode.vcf --remove $OUTPUT/individual_remove_imputation.txt \
--recode-INFO-all --recode --out $OUTPUT/cranberryGBS_snps_to_impute_filtered


## Run beagle

# Run the imputation software
java -Xmx32g -jar /software/7/apps/beagle-geno/5.0/beagle.16May19.351.jar \
gt=$OUTPUT/cranberryGBS_snps_to_impute_filtered.recode.vcf \
out=$OUTPUT/cranberryGBS_germplasm_imputed_snps \
ne=1000 burnin=6 iterations=25

