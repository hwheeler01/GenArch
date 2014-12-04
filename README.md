Collaboration with Haky Im to explore the genetic architecture of cross-tissue and tissue-specific gene expression

## Plan November 2014



## Plan Summer 2014

Genetic Architecture of Transcriptome Regulation Based on Cross Validated Prediction Performance

We will investigate the cis and trans components of genetic regulation of gene expression traits using the prediction performance of gene expression levels

It has been shown that there is a large concentration of eQTLs around the transcription start site of genes. On the other hand, estimates of trans heritability in blood is higher than cis: 75/25% in adipose and 63/37% in blood. This large trans component can be utilized to improve the prediction of gene expression levels.

Superior performance of LASSO relative to polygenic score indicates that in cis, genetic model of gene expression traits is rather sparse.

Using gEUVADIS data only

We will need to do CV to avoid overfitting messing up results. Split into 10 subsets. No need to do this multiple times for now.
Compute cross validated LASSO, i.e. for each test set of people run LASSO using remaining individuals and save the SNPs selected by LASSO
Compute h2 for each gene using GCTA and globalGRM
For each test set, use SNPs from 1. as covariates and run OmicKriging using global GRM. Leave one out. h2vec should be h2 computed for the gene in 2
Compare performance using 
LASSO, 
LASSO+globalGRM, 
LASSO+globalGRM+localGRM, 
localGRM+globalGRM
Compute h2local and h2global using globalGRM and localGRM

References

Price, Alkes L, Agnar Helgason, Gudmar Thorleifsson, Steven A McCarroll, Augustine Kong, and Kari Stefansson. 2011. “Single-Tissue and Cross-Tissue Heritability of Gene Expression via Identity-by-Descent in Related or Unrelated Individuals.” Edited by Greg Gibson. PLoS Genetics 7 (2). Public Library of Science: e1001317. doi:10.1371/journal.pgen.1001317.

Wright, Fred A, Patrick F Sullivan, Andrew I Brooks, Fei Zou, Wei Sun, Kai Xia, Vered Madar, et al. 2014. “Heritability and Genomics of Gene Expression in Peripheral Blood.” Nature Genetics 46 (5). Nature Publishing Group: 430–37. doi:10.1038/ng.2951.

