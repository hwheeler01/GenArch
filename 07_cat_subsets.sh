#!/bin/bash
###concatenate output from 05_est_cross-tissue_h2.r and 06_est_tissue-specifc_h2.r

IFS=$'\n'

for tis in `cat tissue.list`; ##edit tissue.list to add more tissues

do
    echo $tis
    cat GTEx.resid.tissue-specific.h2_${tis}_all.models_subset* >o
    head -n 1 o > header
    sort -t$'\t' -rgk5 o > p
    cat header p > o
    head -n 18737 o > GTEx.resid.tissue-specific.h2_${tis}_all.models_2014-11-24.txt ##removes headers at bottom of file
    rm o p header
done


cat GTEx.cross-tissue.h2.all.models_subset* >o
head -n 1 o > header
sort -rgk5 o > p
cat header p  > GTEx.cross-tissue.h2.all.models_2014-11-25.txt
rm o p header
