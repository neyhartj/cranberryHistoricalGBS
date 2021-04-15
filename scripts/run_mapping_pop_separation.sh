#!/bin/bash

## 
## Subset individuals and parents for each biparental mapping population
## 
## Subset individual and filter SNPs based on minor allele count
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
INPUT=$WD/input/
# Output directory
OUTPUT=$INPUT
# SNP dir
SNPS=$WD/snps/

# Name of the VCF file
VCFIN=$SNPS/cranberryGBS_production_snps_filtered.vcf.gz

# Change working directory
cd $WD

##

# Find the files listing individuals for each population
files=$(find $INPUT -name "family*")

# Iterate over the list of files
for inputfile in $files; do

  # Create the output filename
  outfile=$SNPS/cranberryGBS_$(basename $inputfile | grep -o "family_.*_individuals")
  
  # Run vcftools 
  vcftools --gzvcf $VCFIN --mac 10 --keep $inputfile --recode --recode-INFO-all --out $outfile
  
done



