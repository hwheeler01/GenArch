####by Heather E. Wheeler 20150108####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('22','C')
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
whichlist <- args[2]

getquant <- function(x) quantile(x,c(0.5,0.025,0.975)) ##pulls the median and the 95% credible sets

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

if(whichlist == "A"){
  explist <- explist[1:floor(length(explist)/4)]
}else if (whichlist == "B"){
  explist <- explist[(floor(length(explist)/4)+1):(floor(length(explist)/4)*2)]
}else if (whichlist == "C"){
  explist <- explist[(floor(length(explist)/4)*2+1):(floor(length(explist)/4)*3)]
}else{
  explist <- explist[(floor(length(explist)/4)*3+1):length(explist)]
}

gtfile <- gt.dir %&% 'DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr' %&% chrom %&% '.SNPxID'
gtX <- scan(gtfile)
gtX <- matrix(gtX, ncol = length(fam$V1), byrow=TRUE)
rownames(gtX) <- bim$snp
colnames(gtX) <- fam$V1
X <- gtX[,samplelist]

resultsarray <- array(0,c(length(explist),19))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025","pi025","n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975")
dimnames(resultsarray)[[2]] <- resultscol
working1M <- out.dir %&% "working_" %&% tis %&% "_exp_BSLMM_1M_iterations_chr" %&% chrom %&% whichlist %&% "_" %&% date %&% ".txt"
write(resultscol,file=working1M,ncolumns=19,sep="\t")

working100K <- out.dir %&% "working_" %&% tis %&% "_exp_BSLMM_100K_iterations_chr" %&% chrom %&% whichlist %&% "_" %&% date %&% ".txt"
write(resultscol,file=working100K,ncolumns=19,sep="\t")

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
  if(dim(cisgenos)[1] > 0){
    cisgenos <- mutate(data.frame(cisgenos),snp=rownames(cisgenos))
    cisbim <- inner_join(bim,cisgenos,by='snp')
    annotfile <- cbind(cisbim[,1],cisbim[,4],cisbim[,2])
    genofile <- cbind(cisbim[,1],cisbim[,5:dim(cisbim)[2]])
    phenofile <- data.frame(exp.w.geno[,gene])

    write.table(annotfile, file=out.dir %&% "tmp.annot." %&% chrom %&% "." %&% whichlist, quote=F, row.names=F, col.names=F, sep=",")
    write.table(genofile, file=out.dir %&% "tmp.geno." %&% chrom %&% "." %&% whichlist, quote=F, row.names=F, col.names=F, sep=",")
    write.table(phenofile, file=out.dir %&% "tmp.pheno." %&% chrom %&% "." %&% whichlist, quote=F, row.names=F, col.names=F, sep=",")

    runBSLMM <- "gemma -g " %&% out.dir %&% "tmp.geno." %&% chrom %&% "." %&% whichlist %&% " -p " %&% out.dir %&% "tmp.pheno." %&% chrom %&% "." %&% whichlist %&% " -a " %&% out.dir %&% "tmp.annot." %&% chrom %&% "." %&% whichlist %&% " -bslmm 1 -seed 12345 -o tmp." %&% chrom %&% "." %&% whichlist
    system(runBSLMM)

    hyp <- read.table(out.dir %&% "output/tmp." %&% chrom %&% "." %&% whichlist %&% ".hyp.txt",header=T)
    hyp50 <- hyp[(dim(hyp)[1]/2+1):dim(hyp)[1],] #take second half of sampling iterations
    quantres <- apply(hyp50,2,getquant)
    res <- c(gene,quantres[1,],quantres[2,],quantres[3,])

    hyp10 <- hyp[5001:10000,] #for comparison to 1M (equivalent to second half of 100K iterations)
    quantres10 <- apply(hyp10,2,getquant)
    res10 <- c(gene,quantres10[1,],quantres10[2,],quantres10[3,])

  }else{
    res <- c(gene,rep(NA,18))
    res10 <- c(gene,rep(NA,18))
  }
  names(res) <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025","pi025","n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975") 
  resultsarray[gene,] <- res
  write(res,file=working1M,ncolumns=19,append=T,sep="\t")

  names(res10) <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025","pi025","n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975")
  write(res10,file=working100K,ncolumns=19,append=T,sep="\t")
}

write.table(resultsarray,file=out.dir %&% tis %&% "_exp_BSLMM_1M_iterations_chr" %&% chrom %&% whichlist %&% "_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
