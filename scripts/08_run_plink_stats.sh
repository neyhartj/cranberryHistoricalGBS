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
module load plink/1.9
module load vcftools

## Set variables
# Working directory
WD=/project/gifvl_vaccinium/cranberrygbs/cranberryHistoricalGBS/
# Name of input directory
INPUT=$WD/imputation/beagle_imputation/
# Output directory
OUTPUT=$WD/haplotyping

# Name of input vcf
VCFIN=$INPUT/cranberryGBS_germplasm_imputed_snps.vcf.gz

# Change working directory
cd $WD 



### PLINK HAPLOTYPING ###


# First haplotyping of all germplasm
plink --vcf $VCFIN --blocks no-pheno-req --blocks-max-kb 1000 --double-id \
--out haplotyping/cranberryGBS_germplasm_haplotyping

## Next haplotyping of wild accessions
# First subset the VCF for wild germplasm and write to a plink format file
vcftools --gzvcf $VCFIN --keep $WD/input/cranberry_gbs_wild_germplasm_individuals.txt --plink \
--out cranberryGBS_wild_germplasm_plink

# Now pass to plink
plink --file cranberryGBS_wild_germplasm_plink --blocks no-pheno-req --blocks-min-maf 0.05 \
--blocks-max-kb 1000 --out haplotyping/cranberryGBS_wild_germplasm_haplotypes

# Remove the .ped files
rm cranberryGBS_wild_germplasm_plink.*

