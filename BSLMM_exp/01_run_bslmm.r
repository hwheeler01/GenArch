####by Heather E. Wheeler 20150108####
date <- Sys.Date()
#args <- commandArgs(trailingOnly=T)
args <- '22'
"%&%" = function(a,b) paste(a,b,sep="")
library(dplyr)

###############################################
### Directories & Variables
pre.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
my.dir <- pre.dir %&% "expArch_DGN-WB_imputedGTs/"
rna.dir <- my.dir %&% "dgn-exp/"
annot.dir <- pre.dir %&% "gtex-annot/"
gc12.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/"
gt.dir <- "/group/im-lab/nas40t2/Data/Transcriptome/WB1K/imputed/DGN-imputed-for-PrediXcan/"
out.dir <- pre.dir %&% "BSLMM_exp/"

tis <- "DGN-WB"  
chrom <- as.numeric(args[1]) 
chrname <- "chr" %&% chrom

################################################
rpkmid <- rna.dir %&% tis %&% "exp.ID.list"
expid <- scan(rpkmid,"character")
rpkmgene <- rna.dir %&% tis %&% "exp.GENE.list"
geneid <- scan(rpkmgene,"character")
rpkmfile <- rna.dir %&% tis %&% ".rntransform.exp.IDxGENE"
expdata <- scan(rpkmfile) 
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- expid
colnames(expdata) <- geneid

t.expdata <- expdata #don't need to transpose DGN

gencodefile <- gc12.dir %&% 'gencode.v12.V1.summary.protein.nodup.genenames'
gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,6]
gencode <- gencode[gencode[,1]==chrname,] ##pull genes on chr of interest
t.expdata <- t.expdata[,intersect(colnames(t.expdata),rownames(gencode))] ###pull gene expression data w/gene info
                
expsamplelist <- rownames(t.expdata) ###samples with exp data###

bimfile <- gt.dir %&% "DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr" %&% chrom %&% ".bim" ###get SNP position information###
bim <- read.table(bimfile)
rownames(bim) <- bim$V2
colnames(bim) <- c("chr","snp","cm","bp","ref","alt")
                
famfile <- gt.dir %&% "DGN.hapmap2.chr1-22.QC.fam" ###samples with gt data###
fam <- read.table(famfile)
samplelist <- intersect(fam$V1,expsamplelist)
                        
exp.w.geno <- t.expdata[samplelist,] ###get expression of samples with genotypes###
explist <- colnames(exp.w.geno)

gtfile <- gt.dir %&% 'DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr' %&% chrom %&% '.SNPxID'
gtX <- scan(gtfile)
gtX <- matrix(gtX, ncol = length(fam$V1), byrow=TRUE)
rownames(gtX) <- bim$V2
colnames(gtX) <- fam$V1
X <- gtX[,samplelist]

for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  gene <- explist[i]
  geneinfo <- gencode[gene,]
  chr <- geneinfo[1]
  c <- substr(chr$V1,4,5)
  start <- geneinfo$V3 - 1e6 ### 1Mb lower bound for cis-eQTLS
  end <- geneinfo$V4 + 1e6 ### 1Mb upper bound for cis-eQTLs
  chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
  cissnps <- subset(chrsnps,chrsnps[,4]>=start & chrsnps[,4]<=end) ### pull cis-SNP info
  cisgenos <- X[intersect(rownames(X),cissnps[,2]),,drop=FALSE] ### pull cis-SNP genotypes
  cisgenos <- mutate(data.frame(cisgenos),snp=rownames(cisgenos))
  cisbim <- inner_join(bim,cisgenos,by='snp')
  annotfile <- cbind(cisbim[,1],cisbim[,4],cisbim[,2])
  genofile <- cbind(cisbim[,1],cisbim[,5:dim(cisbim)[2]])
  phenofile <- data.frame(exp.w.geno[,gene])

  write.table(annotfile, file=out.dir %&% "tmp.annot." %&% chrom, quote=F, row.names=F, col.names=F, sep=",")
  write.table(genofile, file=out.dir %&% "tmp.geno." %&% chrom, quote=F, row.names=F, col.names=F, sep=",")
  write.table(phenofile, file=out.dir %&% "tmp.pheno." %&% chrom, quote=F, row.names=F, col.names=F, sep=",")

}
