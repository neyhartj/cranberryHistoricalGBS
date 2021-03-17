#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=01:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=32G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="beagle_imputation_testing"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## Beagle imputation testing
## 
## This script will run Beagle on a set of VCF with random missing
## genotypes using various Beagle parameters. The output will be a set
## of imputed and phased VCF files for comparison with the original.
## 
## 


# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load beagle-geno
module load parallel

## Set variables

# Working directory
WD=/project/cranberrygbs/cranberryHistoricalGBS/imputation/beagle_testing

# Name of input directory
INPUT=$WD/random_missing_files/

# Output directory
OUTPUT=$WD/imputed_files

# Change working directory
cd $WD


## Define a function to be executed in parallel
function runBeagle () { 
  vcfin=$1;
  Ne=$2;
  outdir=$3;
  vcfout=$outdir/$(basename $vcfin | sed 's,.vcf.gz,_Ne'"$Ne"'_imputed,g');
  java -Xmx4g -jar /software/7/apps/beagle-geno/5.0/beagle.16May19.351.jar \
  gt=$vcfin out=$vcfout ne=$Ne burnin=6 iterations=25
};

# Export the function
export -f runBeagle

# List the files
FILES=( $(find $INPUT -name "*.vcf.gz") )


## Beagle Ne 25
parallel -j 8 runBeagle {} 25 $OUTPUT ::: ${FILES[@]}

## Beagle Ne 100
parallel -j 8 runBeagle {} 100 $OUTPUT ::: ${FILES[@]}

## Beagle Ne 1000
parallel -j 8 runBeagle {} 1000 $OUTPUT ::: ${FILES[@]}


