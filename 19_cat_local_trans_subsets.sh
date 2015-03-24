#!/bin/bash
### concatenate output from: 
### 16_est_cross-tissue_local_and_trans_h2.r
### 17_est_tissue-specific_local_and_trans_h2.r
### 18_est_tissue-wide_local_and_trans_h2.r

IFS=$'\n'

## cat tissue-specific
for tis in `cat tissue.list`; ##edit tissue.list to add more tissues
do
    echo $tis
    head -n 1 subsets/GTEx.resid.tissue-specific.h2_${tis}_marginal.local_subset00.2015-03-23.txt > header
    cat subsets/GTEx.resid.tissue-specific.h2_${tis}_marginal.local_subset*.2015-03-23.txt >o
    grep -v tissue o > p
    sort -t$'\t' -rgk5 p > o
    cat header o > GTEx.resid.tissue-specific.h2_${tis}_marginal.local_2015-03-23.txt 
    rm o p header

    head -n 1 subsets/GTEx.resid.tissue-specific.h2_${tis}_FHSfdr0.05.subset00_transForGene.2015-03-23.txt > header
    cat subsets/GTEx.resid.tissue-specific.h2_${tis}_FHSfdr0.05.subset*_transForGene.2015-03-23.txt > o
    grep -v tissue o > p
    sort -t$'\t' -rgk9 p > o
    cat header o > GTEx.resid.tissue-specific.h2_${tis}_FHSfdr0.05_transForGene.2015-03-23.txt
    rm o p header
done

## cat tissue-wide
for tis in `cat gtex-rnaseq/ind-tissues-from-nick/GTEx_PrediXmod.tissue.list`;
do
    echo $tis
    head -n 1 subsets/GTEx.resid.tissue-wide.h2_${tis}_marginal.local_subset00.2015-03-24.txt > header
    cat subsets/GTEx.resid.tissue-wide.h2_${tis}_marginal.local_subset*.2015-03-24.txt >o
    grep -v tissue o > p
    sort -t$'\t' -rgk5 p > o
    cat header o > GTEx.tissue-wide.h2_${tis}_marginal.local_2015-03-24.txt
    rm o p header

    head -n 1 subsets/GTEx.resid.tissue-wide.h2_${tis}_FHSfdr0.05.subset00_transForGene.2015-03-24.txt > header
    cat subsets/GTEx.resid.tissue-wide.h2_${tis}_FHSfdr0.05.subset*_transForGene.2015-03-24.txt > o
    grep -v tissue o > p
    sort -t$'\t' -rgk9 p > o
    cat header o > GTEx.tissue-wide.h2_${tis}_FHSfdr0.05_transForGene.2015-03-24.txt
    rm o p header
done

## cat cross-tissue
head -n 1 subsets/GTEx.cross-tissue.h2.marginal.local_subset00.2015-03-23.txt > header
cat subsets/GTEx.cross-tissue.h2.marginal.local_subset*.2015-03-23.txt >o
grep -v gene o > p
sort -t$'\t' -rgk5 p > o
cat header o > GTEx.ranef.cross-tissue.h2_marginal.local_2015-03-23.txt
rm o p header

head -n 1 subsets/cross-tissue.h2.all.models_FHSfdr0.05.subset00_transForGene.2015-03-23.txt > header
cat subsets/cross-tissue.h2.all.models_FHSfdr0.05.subset*_transForGene.2015-03-23.txt > o
grep -v gene o > p
sort -t$'\t' -rgk9 p > o
cat header o > GTEx.ranef.cross-tissue.h2_FHSfdr0.05_transForGene.2015-03-23.txt
rm o p header


