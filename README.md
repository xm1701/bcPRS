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

1) GWAS summary statistics of trait 1;
2) Sample size of GWAS 
3) Individual-level genetic data of trait 2;
4) SNP heritability estimator of the two traits. 

## Step 1: construct the cross-trait polygnic risk scores.

We recommend to use genetic variants (i.e., SNPs) after linkage disequilibrium (LD)-based pruning (or clumping). No thresholding is required. That is, all pruned SNPs are use to constructe the polygenic risk score. 

## Step 2: construct the cross-trait polygnic risk scores.


## Step 3: Perform bias-correction using the bcPRS package. 

We have tested two options for p 

