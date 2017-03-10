#!/bin/bash
### concatenate output from: 
### 04_est_h2_reml-no-constrain_globalAll.r

IFS=$'\n'

Nk=( 20 40 60 100)
for i in "${Nk[@]}"
do
    head -n 1 subsets/DGN-WB.h2_transPCA-unconstrained_40_PFs_Chr11_globalAll_reml-no-constrain.2017-03-09.txt > header
    cat subsets/DGN-WB.h2_transPCA-unconstrained_${i}_PFs_Chr*_globalAll_reml-no-constrain.2017-03-09.txt > o
    grep -v tissue o > p #rm header rows
    sort -t$'\t' -rgk13 p > o #sort by global h2
    cat header o > DGN-WB.h2_transPCA-unconstrained_${i}_PFs_Chr1-22_globalAll_reml-no-constrain.2017-03-09.txt
    rm o p header
done



