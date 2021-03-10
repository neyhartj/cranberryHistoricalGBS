#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=2   # 1 processor core(s) per node X 2 threads per core
#SBATCH --mem=6400M   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="GBSSeqToTagDBPlugin_testing"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## TASSEL5 pipeline for cranberry GBS
## 
## TagExportToFastqPlugin
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
# Name of database
DBNAME=$WD/database/cranberry_gbs_discovery.db



# Change working directory
cd $WD


# Execute the plugin
run_pipeline.pl -Xms1G -Xmx4G -fork1 -TagExportToFastqPlugin \
-c 1 \
-db $DBNAME \
-o $WD/gbs_tags_for_alignment.fa.gz \
-endPlugin -runfork1

