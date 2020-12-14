# bcPRS

Bias-corrected genetic correlation estimators of cross-trait polygenic risk scores (PRS)

# Description

This tutorial details the steps of obtaining bias-corrected genetic correlation estimator using cross-trait polygenic risk scores. Briefly, the genetic correlation between two complex traits estimated by cross-trait PRS can be biased towarded zero. We quantify the asymptotic limit of this bias and propose the bias-corrected estimator of PRS-based genetic correlation. The preprint can be found at https://arxiv.org/abs/1903.01301. 

# Overview

There are three steps to obtain the bias-corrected estimator of genetic correlation:

1. Construct the cross-trait polygnic risk scores. 

2. Obtain the raw estimator of genetic correlation.

3. Perform bias-correction using the bcPRS package. 


# Example and Demo codes

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

Demo code of LD-based prunning with your genetic data in plink binary format (window size 50, step 5, R2 threshold 0.1): 

~/plink --bfile your_plink_data --indep-pairwise 50 5 0.1 --out list_pruned

~/plink --bfile your_plink_data  --extract list_pruned.prune.in --make-bed  --out your_plink_data_pruned


## Step 1: construct the cross-trait polygenic risk scores.

We recommand the following procedure to generate the cross-trait polygnic risk scores. 

###Case 1: Within one study (e.g., UK Biobank) 

If both of your training and testing daat come from one study (e.g., UK Biobank), you can use either unimputed genotyping genetic variants or imputed genetic variants to construct your PRS. 

Demo code of constructing PRS

~/plink --bfile your_plink_data_pruned   --score  your_gwas_summary_statistics.file  --out prs_scores

In your_gwas_summary_statistics.file, typically we have the following columns: snpid ,A1, A2, and Zscore.
More information about --score function can be found at https://www.cog-genomics.org/plink/1.9/score. 

Note: to perform bias-correction in Step 3, we perform LD-based pruning but no p-value thresholding is required in this step.  

###Case 2: Across two studies (e.g., UK Biobank and a non-UKB study) 

In this situation, we recommend to use imputed genetic variants to increase overlapping rate of genetic variants genotyped in two studies. 
We also need to remove ambiguous genetic variants (i.e. variant with complementary alleles, either C/G or A/T) before generating PRS.

Demo code of removing ambiguous genetic variants (suppose #5 #6 columns are your A1 and A2 data)

awk '!( ($5=="A" && $6=="T") || \
($5=="T" && $6=="A") || \
($5=="G" && $6=="C") || \
($5=="C" && $6=="G")) {print $2}' your_plink_data_pruned.bim > your_plink_data_pruned_atgc.list

## Step 2: Obtain the raw estimator of genetic correlation.

With the PRS of trait 1 generated in Step 1, we can evaludate the genetic correlation between trait 1 and trait 2 in the testing GWAS. 
Typically, we can fit a linear model to estimate the partial correlation between trait 2 and PRS of trait 1, while adjusting the effects of age, gender, and genetic principal components.  

## Step 3: Perform bias-correction using the bcPRS package. 

With this raw estimator of genetic correlation (r_g) estimated in Step 2, we can use the bcPRS package to obtain the bias-corrected estimator. 
For example, suppose the r_g-0.1, training GWAS sample size is 50k, number of indpendent genetic variants is 500k, hertability of the two traits are both 0.5, and the training and testing GWAS are independent (overlapping samples is 0), then we can apply the following function: 

bc_prs(raw_estimator=0.1,n_train=50000,p_indep=500000, h2_trait1=0.5,h2_trait2=0.5,n_overlap=0)

$corrected_estimator
[1] 0.6480741

The bias-corrected PRS-based geneticcorrelation estimator is 0.648. 


