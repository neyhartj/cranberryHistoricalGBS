#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=48G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="tassel5_pipeline_set2"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## TASSEL5 pipeline for cranberry GBS
## 
## Read Alignment
## SAMToGBSdbPlugin
## DiscoverySNPCallerPluginV2
## SNPQualityProfilerPlugin
## 
## 


# Set error handling options
set -e
set -u
set -o pipefail

# Load tassel
module load tassel5
# Load bowtie
module load bowtie2

## Set variables
# Working directory
WD=/project/cranberrygbs/cranberryHistoricalGBS/
# Name of tag fasta
TAGFASTA=$WD/tags/gbs_tags_for_alignment.fa
# Name of output sam file
SAMOUT=$WD/alignment/gbs_tags_aligned_BenLearv2.sam
# SAMOUT1=$WD/alignment/gbs_tags_aligned_Stevensv1.sam
# Basename of reference index
REFIND=/KEEP/cranberrygbs/genome_assemblies/Vm_BenLear_v2_bowtie_index/Vaccinium_macrocarpon_BenLear_v2
# REFIND1=/KEEP/cranberrygbs/genome_assemblies/Vm_Stevens_v1_bowtie_index/Vaccinium_macrocarpon_Stevens_v1

# DB name
DBNAME=$WD/database/cranberry_gbs_discovery.db

# FASTA reference genome
REF=/KEEP/cranberrygbs/genome_assemblies/Vaccinium_macrocarpon_BenLear_v2.fasta

# Output stat file
OUTFILE=$WD/stats/snpStats.txt


# Change working directory
cd $WD

## Read alignment

# Run the alignment
bowtie2 -p 8 --sensitive -x $REFIND -U $TAGFASTA -S $SAMOUT

# # Alternative
# bowtie2 -p 8 --sensitive -x $REFIND1 -U $TAGFASTA -S $SAMOUT1



## SAMToGBSdbPlugin

# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -SAMToGBSdbPlugin  \
-i $SAMOUT \
-db $DBNAME \
-aProp 0.0 \
-aLen 0 \
-minMAPQ 20 \
-endPlugin -runfork1


## DiscoverySNPCallerPluginV2

# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -DiscoverySNPCallerPluginV2 \
-db $DBNAME \
-maxTagsCutSite 1000 \
-mnLCov 0.1 \
-mnMAF 0.00001 \
-ref $REF \
-deleteOldData true \
-endPlugin -runfork1


## SNPQualityProfilerPlugin

# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -SNPQualityProfilerPlugin \
-db $DBNAME \
-statFile $OUTFILE \
-deleteOldData true \
-endPlugin -runfork1







