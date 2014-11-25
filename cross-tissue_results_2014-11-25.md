---
title: Genetic architecture of cross-tissue and tissue-specific gene expression
#output: beamer_presentation
output: pdf_document
author: Heather E. Wheeler
date: 25 November 2014
---



Modeling cross-tissue expression
========================================================
Linear mixed effect model

```r
library(lme4)

fit <- lmer(expression ~ (1|SUBJID) + TISSUE 
+ GENDER + PEERs) 

#cross-tissue expression
fitranef <- ranef(fit) 

#tissue-specific expression
fitresid <- resid(fit) 
```

Estimating heritability with GCTA
========================================================

Tested several genetic relationship matrix (GRM) models for each expressed gene

- localGRM (SNPs within 1 Mb of gene)
- localGRM + globalGRM (all SNPs)
- localGRM + chrGRM + globalGRM

First pass: estimated h2 of cross-tissue expression and tissue-specific expression in the 7 tissues with the most samples

GCTA heritability: Y ~ localGRM h2
========================================================
![alt text](plots/hist_GTEx_localGRM_h2_2014-11-25.pdf)

GCTA heritability: Y ~ localGRM h2 **ZOOM**
========================================================
![alt text](plots/hist_GTEx_localGRM_h2_ylim3000_2014-11-25.pdf)

GCTA heritability: Y ~ localGRM p-values
========================================================
![alt text](plots/hist_GTEx_localGRM_p_2014-11-25.pdf)

GCTA heritability: Y ~ localGRM + globalGRM h2 
========================================================
![alt text](plots/scatter_GTEx_localGRM_globalGRM_h2_2014-11-25.pdf)

GCTA heritability: Y ~ localGRM + globalGRM h2 
========================================================
![alt text](plots/hist_GTEx_localGRM_globalGRM_h2_2014-11-25.pdf)

GCTA heritability: Y ~ localGRM + globalGRM SE 
========================================================
![alt text](plots/hist_GTEx_localGRM_globalGRM_SE_2014-11-25.pdf)

GCTA heritability: Y ~ localGRM + chrGRM + globalGRM h2 
======================================================== 
![alt text](plots/scatter_GTEx_localGRM_chrGRM_globalGRM_h2_2014-11-25.pdf)

GCTA heritability: Y ~ localGRM + chrGRM + globalGRM h2 
========================================================
![alt text](plots/hist_GTEx_localGRM_chrGRM_globalGRM_h2_2014-11-25.pdf)

GCTA heritability: Y ~ localGRM + chrGRM + globalGRM SE 
========================================================
![alt text](plots/hist_GTEx_localGRM_chrGRM_globalGRM_SE_2014-11-25.pdf)
