#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 1 processor core(s) per node X 2 threads per core
#SBATCH --mem=48G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="tassel5_production_pipeline"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## TASSEL5 pipeline for cranberry GBS
## 
## ProductionSNPCallerPluginV2
## 
## 


# Set error handling options
set -e
set -u
set -o pipefail

# Load the TASSEL module
module load tassel5

## Set variables
# Working directory
WD=/project/cranberrygbs/cranberryHistoricalGBS/

# Name of input directory
INPUT=$WD/input
# Name of input directory with FASTQ
FASTQDIR=$INPUT/fastq_files
# Name of database
DBNAME=$WD/database/cranberry_gbs_discovery.db
# Name of keyfile
# KEY=$INPUT/cranberry_gbs_unique_keys.txt
KEY=$INPUT/cranberry_gbs_unique_keys_resolved_duplicates.txt


# Output file name
OUTFILE=$WD/snps/cranberryGBS_production_snps.vcf


# Change working directory
cd $WD


# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -ProductionSNPCallerPluginV2 \
-db $DBNAME \
-e EcoT22I \
-i $FASTQDIR \
-k $KEY \
-kmerLength 64 \
-do true \
-o $OUTFILE \
-endPlugin -runfork1

