
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Cranberry GBS - Scripts

<!-- badges: start -->
<!-- badges: end -->

The subfolder contains the scripts used to process reads generated via
genotyping-by-sequencing and call SNP variants. The order and
description of the scripts is below:

## Read processing and variant calling

1.  `01_prepare_keyfile.R` - removes blank samples and assigns
    duplicated samples (based on ONLY information in the original
    keyfile) the same sample name.

2.  `02_run_tassel5_pipeline_allUniqueKeys.sh` - runs the entire TASSEL5
    variant-calling pipeline. This version uses the keyfile
    corresponding to all unique entries, which have been grouped based
    on ONLY information in the original keyfile (i.e. explicitly
    duplicated sample names.)

3.  `03_filter_variants_allUniqueKeys.sh` - filters the variants
    generated from Step 2 based on the following: MinGQ = 40, MinMeanDP
    = 10, MinDP = 7, and bi-allelic SNP variants only.

4.  `04_resolve_sample_duplicates.R` - uses the variants generated from
    Step 2 to identify accessions that are likely duplicates; also
    identifies and merges duplicates based on previous DNA
    fingerprinting using “SCAR” markers.

5.  `05_run_tassel5_pipeline_resolvedDuplicates.sh` - re-run the TASSEL5
    pipeline after merging duplicate samples in Step 4.

6.  `06_filter_variants_resolvedDuplicates.sh` - filters the variates
    discovered in Step 5 using the same criteria as in Step 3.

## Imputation and other post-processing

7.  `07_run_beagle_imputation.sh` - runs BEAGLE to impute variants.

8.  `08_run_plink_stats.sh` - generates some haplotype stats using
    PLINK.
