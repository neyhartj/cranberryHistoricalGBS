#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 1 processor core(s) per node X 2 threads per core
#SBATCH --mem=6400M   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="ReadAlignment"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## TASSEL5 pipeline for cranberry GBS
## 
## Read Alignment
## 
## 


# Set error handling options
set -e
set -u
set -o pipefail

# Load bowtie
module load bowtie2

## Set variables
# Working directory
WD=/project/cranberrygbs/cranberryHistoricalGBS/testing/
# Name of tag fasta
TAGFASTA=$WD/tags/gbs_tags_for_alignment.fa.gz
# Name of output sam file
SAMOUT=$WD/tags/gbs_tags_aligned.sam
# Basename of reference index
REFIND=/KEEP/cranberrygbs/genome_assemblies/Vm_BenLear_v2_bowtie_index/Vaccinium_macrocarpon_BenLear_v2

# Change working directory
cd $WD

# Run the alignment
bowtie2 \
-p 8 \ 
-f \
--sensitive \
-x $REFIND \ 
-U $TAGFASTA \
-S $SAMOUT









