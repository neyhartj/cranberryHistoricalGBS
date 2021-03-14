#!/bin/bash

# Variant filtration

##### Make changes below this line #####

# Working directory
WD=/project/cranberrygbs/cranberryHistoricalGBS/
# Name of input directory
INPUT=$WD/input

# The complete filepath to the VCF file
VCFIN=$WD/snps/cranberryGBS_production_snps.vcf

# The desired name of the VCF file (including the .vcf extension)
## (e.g. '2row_GBS_filtered_snps.vcf')
VCFOUT=$WD/snps/cranberryGBS_production_snps_filtered.vcf

# The keyfile produced in the demultiplexing step. Note: this file should
## look like the following: Samples+barcode+ID_${PROJECT}.txt
## One may provide multiple keyfiles, separated by a comma (',').
KEYFILE=$INPUT/cranberry_gbs_all_keys.txt

# You may change the filtering parameters below:
# Note: by default, this filtering script will extract only biallelic 
## SNP variants.

# Default parameters are assigned already 
# Any variables left blank will be assigned default values 

# Minimum mapping quality score (phred scaled) # to include that variant
MinQ=0
# Minimum genotype quality score (phred scaled) to include that genotype
MinGQ=40
# The minimum average read depth for a variant to include that variant
MinMeanDP=10
# The minimum allele depth for any genotype
MinDP=7
# The minimum minor allele frequency to retain a site
MinMAF=0
# The minimum non-missing data proportion to retain a site (where 0 allows 
## completely missing data and 1 retricts to no missing data)
MaxMISS=0

# Other parameters can be changed as well. You may edit the script after the line beginning
## with "bcftools reheader...". Please see https://vcftools.github.io/man_latest.html for
## a list of adjustable parameters.


####################################################
##### USE CAUTION WHEN EDITING BELOW THIS LINE #####
####################################################

set -e
set -u
set -o pipefail

# Load modules
module load bcftools
module load vcftools


# Change working directory
cd $WD

## Run pre-filtering stats
bcftools stats $VCFIN > $WD/snps/cranberryGBS_production_snps_prefilter_stats.out


## Remove blank samples using the keyfile

# First, find the appropriate column name in the keyfile
COLS=1
for i in $(head $KEYFILE -n 1) ; do
      if [ $i == "FullSampleName" ] ; then
        break ;
    else
        COLS=$(( $COLS + 1 ))
    fi
done

# Cut the file for those sample names
BLANKS=$(cut -f $COLS $KEYFILE | grep "BLANK\|blank" -)

# Capture the BLANK samples in a variable
# THIS NEEDS TO BE CHANGED
SAMPLESREMOVE=$(echo $BLANKS | sed 's/ / --remove-indv /g')

# Use bcftools reheader to change the sample names, then pipe the output to vcftools for filtering
echo -e "\nStarting VCF filtering."
vcftools --vcf $VCFIN \
	--remove-indels \
	--remove-filtered-all \
	--remove-indv ${SAMPLESREMOVE} \
	--min-alleles 2 \
	--max-alleles 2 \
	--min-meanDP $MinMeanDP \
	--minDP $MinDP \
	--minGQ $MinGQ \
	--maf $MinMAF \
	--max-missing $MaxMISS \
	--recode \
	--recode-INFO-all \
	--out $VCFOUT \
&& echo -e "\nFiltering complete and output file created."

# Rename the file
mv $VCFOUT.recode.vcf $VCFOUT

## Run post-filtering stats
bcftools stats $VCFOUT > $WD/snps/cranberryGBS_production_snps_postfilter_stats.out
