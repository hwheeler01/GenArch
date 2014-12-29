#!/bin/bash
###concatenate output from 09_perms_est_cross-tissue_h2.r and 10_perms_est_tissue-specifc_h2.r

IFS=$'\n'

for tis in `cat tissue.list.prefix`; ##edit tissue.list to add more tissues

do
    echo $tis
    cat ${tis}_h2_globaljoint_subset*_100perms* > o
    sort -t$'\t' -rgk5 o > ${tis}_h2_globaljoint_allgenes_100perms_2014-12-19.txt

    cat ${tis}_h2_globalonly_subset*_100perms* > o
    sort -t$'\t' -rgk5 o > ${tis}_h2_globalonly_allgenes_100perms_2014-12-19.txt

    cat ${tis}_h2_localjoint_subset*_100perms* > o
    sort -t$'\t' -rgk5 o > ${tis}_h2_localjoint_allgenes_100perms_2014-12-19.txt

    cat ${tis}_h2_localonly_subset*_100perms* > o
    sort -t$'\t' -rgk5 o > ${tis}_h2_localonly_allgenes_100perms_2014-12-19.txt

done
