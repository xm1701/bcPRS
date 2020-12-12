# bcPRS

Bias-corrected genetic correlation estimators of cross-trait polygenic risk scores (PRS)

# Description

This tutorial details the steps of obtaining bias-corrected genetic correlation estimator using cross-trait polygenic risk scores. The preprint could be found at https://arxiv.org/abs/1903.01301. Briefly, the genetic correlation between two complex traits estimated by cross-trait PRS can be biased towarded zero. We quantify the asymptotic limit of this bias and propose the bias-corrected estimator of genetic correlation. 

# Overview

There are three steps to obtain the bias-corrected estimator of genetic correlation:

1. Construct the cross-trait polygnic risk scores. 

2. Obtain the raw estimator of genetic correlation.

3. Perform bias-correction using the bcPRS package. 


# Example with demo codes

## Step 0: preparation of input data
To generate cross-trait PRS and obtain the bias-corrected genetic correlation estimator between trait 1 and trait 2, we need the following data:

1) GWAS summary statistics of trait 1. 
Regular full GWAS summary statistics dataset available from GWAS data consortium, such as GWAS catalog https://www.ebi.ac.uk/gwas/downloads/summary-statistics.

2) Sample size of the above GWAS.
Such information is typically available along with the GWAS summary statistics dataset or can be found in the reference paper. 

3) Individual-level genetic data of trait 2.
This is the in-house individual-level GWAS dataset that you have access to. 

4) SNP heritability estimator of the two traits. 
For example, the heritability estimated using the GREML method https://cnsgenomics.com/software/gcta/#Overview. The GREML estimator of many complex traits have been made publicly available, such as from https://nealelab.github.io/UKBB_ldsc/ and https://atlas.ctglab.nl/. 

5) Number of independent genetic variants. This can be obtained by performing LD-based prunning or clumping via plink (https://www.cog-genomics.org/plink2/) on your individual-level genetic data. 

Demo code: 

~/plink --bfile your_plink_data --indep-pairwise 50 5 0.1 --out list_pruned
~/plink --bfile your_plink_data  --extract list_pruned.prune.in --make-bed  --out your_plink_data_pruned


## Step 1: construct the cross-trait polygnic risk scores.

We recommand the following two options to generate the cross-trait polygnic risk scores. 

###Option 1: Within UK Biobank 

If you are using UK Biobank data, we suggest the number of independent genetic variants to be 150k. 

We have tested two options for p 

We recommend to use genetic variants (i.e., SNPs) after linkage disequilibrium (LD)-based pruning (or clumping). No thresholding is required. That is, all pruned SNPs are use to constructe the polygenic risk score. 

## Step 2: construct the cross-trait polygnic risk scores.


## Step 3: Perform bias-correction using the bcPRS package. 



