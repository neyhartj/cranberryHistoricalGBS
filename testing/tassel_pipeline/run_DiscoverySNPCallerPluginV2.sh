#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=2   # 1 processor core(s) per node X 2 threads per core
#SBATCH --mem=6400M   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="DiscoverySNPCallerPluginV2"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## TASSEL5 pipeline for cranberry GBS
## 
## DiscoverySNPCallerPluginV2
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
WD=/project/cranberrygbs/cranberryHistoricalGBS/testing/
# GBS database created in step 1
DBNAME=$WD/database/cranberry_gbs_discovery.db
# FASTA reference genome
REF=/KEEP/cranberrygbs/genome_assemblies/Vaccinium_macrocarpon_BenLear_v2.fasta


# Change working directory
cd $WD


# Execute the plugin
run_pipeline.pl -Xms1G -Xmx4G -fork1 -DiscoverySNPCallerPluginV2 \
-db $DBNAME \
-maxTagsCutSite 1000 \
-mnLCov 0.1 \
-mnMAF 0.00001 \
-ref $REF \
-deleteOldData true \
-endPlugin -runfork1

