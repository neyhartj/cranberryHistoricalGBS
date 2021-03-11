# Script for processing the statistics from the SNP SNPQualityProfilerPlugin
# 
# 

# Load packages and set the directory
library(tidyverse)

proj_dir <- here::here()
pipeline_dir <- file.path(proj_dir, "testing")
stat_dir <- file.path(pipeline_dir, "stats")


# Read in the data --------------------------------------------------------

# Load the data
snp_stats_df <- read_tsv(file = file.path(stat_dir, "snpStats.txt"))

# Tidy
snp_stats_df1 <- snp_stats_df %>%
  arrange(Chromosome, PositionID)


# Read in chromosome lengths
chromLens <- read_tsv(file = "../../ReferenceGenomes/Vaccinium_macrocarpon_BenLear_v2_chromlens.txt",
                      col_names = c("Chromosome", "Length")) %>%
  arrange(Chromosome) %>%
  mutate(Chromosome = str_remove(Chromosome, "chr") %>% toupper)


# Summaries ---------------------------------------------------------------

# SNPs per chrom
snp_stats_df1 %>%
  group_by(Chromosome) %>%
  summarize(nSnps = n()) %>%
  ggplot(aes(x = Chromosome, y = nSnps)) +
  geom_col()

# Plot snp positions
ggplot(data = chromLens, aes(y = Chromosome, yend = Chromosome)) +
  geom_segment(aes(x = 0, xend = Length), lwd = 1, color = "grey85") + 
  # geom_point(data = snp_stats_df1, aes(x = PositionID), size = 0.25) + 
  geom_segment(data = snp_stats_df1, 
               aes(x = PositionID, xend = PositionID, y = as.numeric(Chromosome) - 0.2, yend = as.numeric(Chromosome) + 0.2), 
               lwd = 0.01) +
  theme_classic()

# Remove the unknown chrom
snp_stats_df2 <- filter(snp_stats_df1, Chromosome != "UNKNOWN")


# Average read depth per site
snp_stats_df2 %>%
  ggplot(aes(x = avgDepth)) +
  geom_histogram(bins = 50) +
  facet_wrap(~ Chromosome, scales = "free_x")

# Percentage of total depth consisting of the minor allele
snp_stats_df2 %>%
  ggplot(aes(x = minorDepthProp)) +
  geom_histogram(bins = 50) +
  facet_wrap(~ Chromosome, scales = "free_x")

# Number of multi-allelic SNPS
snp_stats_df2 %>%
  filter(minor2DepthProp != 0) %>%
  group_by(Chromosome) %>%
  summarize(nSnps = n())










