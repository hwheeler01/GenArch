####by Heather E. Wheeler 20150108####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('22','1')
"%&%" = function(a,b) paste(a,b,sep="")

###############################################
### Directories & Variables

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
ct.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/"
annot.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/"
gt.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/DGN-WB_genotypes/by.chr/"
en.dir <- my.dir %&% "DGN-weights/"

k <- 10 ### k-fold CV
n <- 1 #number of k-fold CV replicates, remove nrep loop for this implementation
tis <- "DGN-WB"  

##alpha = The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
#(1-α)/2||β||_2^2+α||β||_1.
#alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.

alpha <- as.numeric(args[2]) #alpha to test in CV

gencodefile <- annot.dir %&% 'gencode.v12.V1.summary.protein.nodup.genenames'
chrom <- args[1]
chrname <- "chr" %&% chrom

snpset <- "hapmapSnpsCEU"

################################################
### Functions & Libraries

library(glmnet)
#library(doMC)
#registerDoMC(4)
#getDoParWorkers()

################################################
rpkmid <- ct.dir %&% tis %&% ".exp.ID.list"
expid <- scan(rpkmid,"character")
rpkmgene <- ct.dir %&% tis %&% ".exp.GENE.list"
geneid <- scan(rpkmgene,"character")
rpkmfile <- ct.dir %&% tis %&% ".exp.IDxGENE"
expdata <- scan(rpkmfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- expid
colnames(expdata) <- geneid

t.expdata <- expdata #don't need to transpose DGN

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,6]
gencode <- gencode[gencode[,1]==chrname,] ##pull genes on chr of interest
t.expdata <- t.expdata[,intersect(colnames(t.expdata),rownames(gencode))] ###pull gene expression data w/gene info
                
expsamplelist <- rownames(t.expdata) ###samples with exp data###
              	

bimfile <- gt.dir %&% "DGN.hapmap2.QC.chr" %&% chrom %&% ".bim" ###get SNP position information###
bim <- read.table(bimfile)
rownames(bim) <- bim$V2

famfile <- gt.dir %&% "DGN.hapmap2.QC.chr" %&% chrom %&% ".fam" ###samples with gt data###
fam <- read.table(famfile)
samplelist <- intersect(fam$V1,expsamplelist)

exp.w.geno <- t.expdata[samplelist,] ###get expression of samples with genotypes###
explist <- colnames(exp.w.geno)

gtfile <- gt.dir %&% tis %&% '.gt.chr' %&% chrom %&% '.IDxSNP'
gtX <- scan(gtfile)
gtX <- matrix(gtX, ncol = length(bim$V2), byrow=TRUE)
colnames(gtX) <- bim$V2
rownames(gtX) <- fam$V1
X <- gtX[samplelist,]

set.seed(42)
groupid <- sample(1:10,length(samplelist),replace=TRUE) ##need to use same folds to compare alphas

resultsarray <- array(0,c(length(explist),8))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval")
dimnames(resultsarray)[[2]] <- resultscol
workingbest <- "working_" %&% tis %&% "_exp_" %&% k %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% "_" %&% date %&% ".txt"
write(resultscol,file=workingbest,ncolumns=8,sep="\t")

weightcol = c("gene","SNP","refAllele","effectAllele","beta")
workingweight <- en.dir %&% tis %&% "_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_weights_chr" %&% chrom %&% "_" %&% date %&% ".txt"
write(weightcol,file=workingweight,ncol=5,sep="\t")


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
  cisgenos <- X[,intersect(colnames(X),cissnps[,2])] ### pull cis-SNP genotypes
  if(is.null(dim(cisgenos))){
    bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
  }else{
    minorsnps <- subset(colMeans(cisgenos), colMeans(cisgenos,na.rm=TRUE)>0) ###pull snps with at least 1 minor allele###
    minorsnps <- names(minorsnps)
    cisgenos <- cisgenos[,minorsnps]
    if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs
      bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else{

      exppheno <- exp.w.geno[,gene] ### pull expression data for gene
      exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
      exppheno[is.na(exppheno)] <- 0
      rownames(exppheno) <- rownames(exp.w.geno)
  
      ##run Cross-Validation over alphalist
      fit <- cv.glmnet(cisgenos,exppheno,nfolds=k,alpha=alpha,keep=T,foldid=groupid,parallel=F) ##parallel=T is slower on tarbell, not sure why

      fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm)) ##pull info to find best lambda
      best.lam <- fit.df[which.min(fit.df[,1]),] # needs to be min or max depending on cv measure (MSE min, AUC max, ...)
      cvm.best = best.lam[,1]
      lambda.best = best.lam[,2]
      nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output
      
      ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best]) # get betas from best lambda
      ret[ret == 0.0] <- NA
      bestbetas = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
      names(bestbetas) = rownames(ret)[which(!is.na(ret))]

      pred.mat <- fit$fit.preval[,nrow.best] # pull out predictions at best lambda

    }
  }
  if(length(bestbetas) > 0){
    res <- summary(lm(exppheno~pred.mat))
    genename <- as.character(gencode[gene,6])
    rsq <- res$r.squared
    pval <- res$coef[2,4]

    resultsarray[gene,] <- c(genename, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval)

    
    ### output best shrunken betas for PrediXcan
    bestbetalist <- names(bestbetas)
    bestbetainfo <- bim[bestbetalist,]
    betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
    betafile<-cbind(genename,betatable[,2],betatable[,5],betatable[,6],betatable[,7]) ##output "gene","SNP","refAllele","effectAllele","beta"
    write(t(betafile),file=workingweight,ncolumns=5,append=T,sep="\t") # t() necessary for correct output from write() function

  }else{
    genename <- as.character(gencode[gene,6])
    resultsarray[gene,1] <- genename
    resultsarray[gene,2:8] <- c(NA,NA,NA,NA,0,NA,NA)

  }
  write(resultsarray[gene,],file=workingbest,ncolumns=8,append=T,sep="\t")
}


write.table(resultsarray,file=tis %&% "_exp_" %&% k %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% "_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
