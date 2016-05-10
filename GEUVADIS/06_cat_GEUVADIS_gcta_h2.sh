#!/bin/bash
## cat pop-specific output from 05_est_gcta_h2.r

IFS=$'\n'

head -n 1 chr_results/GEUVADIS_FIN_exp_gcta_reml_local_h2_1KGP_geuMAF0.01_chr9_2016-04-28.txt > header

for pop in `cat poplist`; 
do
    echo $pop
    cat chr_results/GEUVADIS_${pop}_exp_gcta_reml_local_h2_1KGP_geuMAF0.01_chr*.txt > o
    grep -v ensid o > p
    cat header p > GEUVADIS_${pop}_exp_gcta_reml_local_h2_1KGP_geuMAF0.01_all_chr1-22_2016-05-02.txt
    rm o p
done

rm header 

