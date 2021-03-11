#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=48G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="Tassel5_pipeline_steps1-2"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## TASSEL5 pipeline for cranberry GBS
## 
## GBSSeqToTagDBPlugin
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
# KEY=$INPUT/cranberry_gbs_all_keys.txt
KEY=$INPUT/cranberry_gbs_all_keys.txt
# Directory to store fa file
OUTDIR=$WD/tags



# Change working directory
cd $WD

## GBSSeqToTagDBPlugin
# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -GBSSeqToTagDBPlugin \
-c 10 \
-db $DBNAME \
-i $INPUT \
-k $KEY \
-e EcoT22I \
-kmerLength 64 \
-minKmerL 20 \
-mnQS 20 \
-mxKmerNum 100000000 \
-endPlugin -runfork1

## TagExportToFastqPlugin 
# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -TagExportToFastqPlugin \
-c 1 \
-db $DBNAME \
-o $OUTDIR/gbs_tags_for_alignment.fa.gz \
-endPlugin -runfork1

# Unzip the fasta file
gunzip $OUTDIR/gbs_tags_for_alignment.fa.gz


