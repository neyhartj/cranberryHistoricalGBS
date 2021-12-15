#!/bin/bash

## 
## Beagle imputation testing preparation
## 
## This script will prepare files for beagle imputation
## testing
## 
## 


# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load vcftools

## Set variables

# Working directory
WD=/project/cranberrygbs/cranberryHistoricalGBS/

# Name of input directory
INPUT=$WD/snps/

# Output directory
OUTPUT=$WD/imputation/beagle_testing

# Change working directory
cd $WD

# First create a set of high-confidence SNPs by filtering:
vcftools --vcf $INPUT/cranberryGBS_production_snps_filtered.vcf --minDP 25 --mac 25 \
--max-missing 0.9 --recode-INFO-all --recode --out $OUTPUT/snps_for_imputation_testing

# Calculated individual missingness
vcftools --vcf $OUTPUT/snps_for_imputation_testing.recode.vcf --missing-indv --out $OUTPUT/snps_for_imputation_testing

# Filtered individuals with <= 5% missing data
awk 'NR > 1 { if ($5 > 0.05) print $1 }' $OUTPUT/snps_for_imputation_testing.imiss > $OUTPUT/individual_remove_imputation_test.txt

vcftools --vcf $OUTPUT/snps_for_imputation_testing.recode.vcf --remove $OUTPUT/individual_remove_imputation_test.txt \
--recode-INFO-all --recode --out snps_for_imputation_testing_indv_filtered
