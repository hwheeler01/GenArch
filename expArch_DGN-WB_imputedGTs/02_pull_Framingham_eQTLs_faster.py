#!/usr/bin/env python

'''pull Framingham eQTL SNP lists at various thresholds: p<0.0001 (entire file), FDR<0.05
For now, only include hapmap2 SNPs because that's what we selected from the DGN imputation'''

import sys, subprocess

#framing = '/group/im-lab/nas40t2/haky/Data/Transcriptome/eQTL-results/Framingham/eqtl-gene-1000g-peer-validate_c1-23_c28.csv.gz'
framing = '10k.test.csv.gz'
hapmap = '/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/hapmapSnpsCEU.list' ##from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz


outfileAll = open('fast.Framingham_eqtl-gene_p0.0001_hapmapSnpsCEU.rsIDlist','w')
outfileFDR = open('fast.Framingham_eqtl-gene_fdr0.05_hapmapSnpsCEU.rsIDlist','w')

haplist = open(hapmap).read().splitlines()  ##put each line of file in list with no newline

##function to call zcat and return line from http://spyced.blogspot.com/2006/12/wow-gzip-module-kinda-sucks.html
def gziplines(fname):
  from subprocess import Popen, PIPE
  f = Popen(['zcat', fname], stdout=PIPE)
  for line in f.stdout:
      yield line

for line in gziplines(framing):
    linelist = line.split(',')
    snp = linelist[3]
    if snp == 'Rs_ID':
        header_line = line
    elif snp in haplist:
        outfileAll.write(snp+'\n')
        logfdr = float(linelist[22])
        fdr = 10**logfdr
        if fdr < 0.05:
            outfileFDR.write(snp+'\n')
        
outfileAll.close()
outfileFDR.close()
