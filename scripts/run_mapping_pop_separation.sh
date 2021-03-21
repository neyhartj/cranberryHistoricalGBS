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
INPUT=$WD/snps/

# Output directory
OUTPUT=$INPUT

# Change working directory
cd $WD

