#!/bin/bash
## cat tissue-specific output from 26_GTEx_TW_CV_elasticNet.r

IFS=$'\n'

head -n 1 TW_Lung_exp_10-foldCV_elasticNet_alpha0.05_hapmapSnpsCEU_chr10_2015-08-20.txt > header

head -n 1 gtex-OTD-weights/TW_Lung_elasticNet_alpha0.05_hapmapSnpsCEU_weights_chr10_2015-08-20.txt >weightheader

for tis in `cat nine.tissue.list`; 
do
    echo $tis
    for i in 0.05 0.5 0.95 1;
    do
	cat TW_${tis}_exp_10-foldCV_elasticNet_alpha${i}_hapmapSnpsCEU_chr*.txt >o
	grep -v gene o > p
	cat header p > TW_${tis// /}_exp_10-foldCV_elasticNet_alpha${i}_hapmapSnpsCEU_all_chr1-22_2015-08-27.txt
	rm o p

	cat gtex-OTD-weights/TW_${tis}_elasticNet_alpha${i}_hapmapSnpsCEU_weights_chr*.txt >o
	grep -v gene o > p
	cat weightheader p > gtex-OTD-weights/TW_${tis// /}_elasticNet_alpha${i}_hapmapSnpsCEU_weights_all_chr1-22_2015-08-27.txt
	rm o p 
    done
done

rm header weightheader
