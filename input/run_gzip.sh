#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=08:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=2   # 1 processor core(s) per node X 2 threads per core
#SBATCH --mem=6400M   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="gzip fastq"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## gzip FASTQ files
##


# Set error handling options
set -e
set -u
set -o pipefail

# Directory containing fastq files
FASTQDIR=/project/cranberrygbs/cranberryHistoricalGBS/input/fastq_files

# Change into that directory
cd $FASTQDIR

# Search for fastq files that are not compressed
files=$(find . -name "*_fastq")

# Loop over files
for fi in $files; do
  # gzip
  gzip $fi
  
done




