#!/bin/bash
## cat pop-specific output from 05_est_gcta_h2.r

IFS=$'\n'

head -n 1 geu_bslmm/GEUVADIS_ALL_BSLMM-s100K_iterations_ALL_chr10A_2016-05-03.txt > header

for pop in `cat poplist`; 
do
    echo $pop
    cat geu_bslmm/GEUVADIS_${pop}_BSLMM-s100K_iterations_${pop}_chr*.txt > o
    grep -v gene o > p
    cat header p > GEUVADIS_${pop}_BSLMM-s100K_iterations_all_chr1-22_2016-05-11.txt
    rm o p
done

rm header 

