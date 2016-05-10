#!/bin/bash
## cat pop-specific output from 01_GEUVADIS_CV_elasticNet.r

IFS=$'\n'

head -n 1 GEUVADIS_GBR_exp_10-foldCV_elasticNet_alpha0_1KGP_geuMAF0.01_chr10_2016-05-09.txt > header

head -n 1 geu_weights/GEUVADIS_GBR_elasticNet_alpha0_1KGP_geuMAF0.01_weights_chr10_2016-05-09.txt >weightheader

for pop in `cat poplist`; 
do
    echo $pop
    for i in 0 0.5 1;
    do
	cat GEUVADIS_${pop}_exp_10-foldCV_elasticNet_alpha${i}_1KGP_geuMAF0.01_chr*.txt > o
	grep -v gene o > p
	cat header p > GEUVADIS_${pop}_exp_10-foldCV_elasticNet_alpha${i}_1KGP_geuMAF0.01_all_chr1-22_2016-05-10.txt
	rm o p

	cat geu_weights/GEUVADIS_${pop}_elasticNet_alpha${i}_1KGP_geuMAF0.01_weights_chr*.txt > o
	grep -v gene o > p
	cat weightheader p > geu_weights/GEUVADIS_${pop}_elasticNet_alpha${i}_1KGP_geuMAF0.01_weights_all_chr1-22_2016-05-10.txt
	rm o p 
    done
done

rm header weightheader
