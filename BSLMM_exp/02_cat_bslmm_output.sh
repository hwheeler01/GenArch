#!/bin/bash
## cat output from 01_dgn_bslmm.r

head -n 1 working_DGN-WB_exp_BSLMM_1M_iterations_chr1A_2015-06-11.txt > header
cat working_DGN-WB_exp_BSLMM_1M_iterations_chr* > o
grep -v gene o > p
cat header p > working_DGN-WB_exp_BSLMM_1M_iterations_genes_finished_in_100hrs_2015-06-16.txt
rm o p header

head -n 1 DGN-WB_exp_BSLMM-s100K_iterations_chr22A_2015-06-14.txt > header
cat DGN-WB_exp_BSLMM-s100K_iterations_chr* > o
grep -v gene o > p
cat header p > DGN-WB_exp_BSLMM-s100K_iterations_all_genes_2015-06-14.txt
rm o p header
