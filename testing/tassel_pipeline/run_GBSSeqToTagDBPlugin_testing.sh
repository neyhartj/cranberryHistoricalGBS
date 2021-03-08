#!/bin/bash

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
WD=/project/cranberrygbs/cranberryHistoricalGBS/testing/
# Name of input directory with fastqs
INPUT=$WD/input/fastq_files
# Name of database
DBNAME=$WD/database/cranberry_gbs_discovery.db
# Name of keyfile
KEY=$INPUT/cranberry_gbs_all_keys.txt



# Execute the plugin
./run_pipeline.pl -Xms1G -Xmx4G -fork1 -GBSSeqToTagDBPlugin \
-c 10 \
-db $DBNAME \
-i $INPUT \
-k $KEY \
-kmerLength 64 \
-minKmerL 20 \
-mnQS 20 \
-mxKmerNum 100000000 \
-endPlugin -runfork1

