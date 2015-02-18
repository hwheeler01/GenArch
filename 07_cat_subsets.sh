#!/bin/bash
###concatenate output from 05_est_cross-tissue_h2.r and 06_est_tissue-specifc_h2.r

IFS=$'\n'

for tis in `cat tissue.list`; ##edit tissue.list to add more tissues

do
    echo $tis
    cat subsets/GTEx.resid.tissue-specific.h2_${tis}_all.models_subset*2-17* >o
    head -n 1 o > header
    sort -t$'\t' -rgk5 o > p
    cat header p > o
    head -n 17023 o > GTEx.resid.tissue-specific.h2_${tis}_all.models_2015-02-17.txt ##removes headers at bottom of file
    rm o p header
done


cat subsets/GTEx.cross-tissue.h2.all.models_subset*2-17* >o
head -n 1 o > header
sort -rgk5 o > p
cat header p  > o
head -n 17023 o > GTEx.cross-tissue.h2.all.models_2015-02-17.txt
rm o p header
