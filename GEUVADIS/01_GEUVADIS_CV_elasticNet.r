####by Heather E. Wheeler 20160412####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('22','0.5','GBR')
"%&%" = function(a,b) paste(a,b,sep="")

###############################################
### Directories & Variables

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
exp.dir <- my.dir %&% "GEUVADIS/"
annot.dir <- my.dir %&% "gtex-annot/"
gt.dir <- exp.dir %&% "geu_dosage_1KGP_MAF0.01/"
en.dir <- exp.dir %&% "geu_weights/"

k <- 10 ### k-fold CV
n <- 1 #number of k-fold CV replicates
pop <- args[3] 

##alpha = The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
#(1-α)/2||β||_2^2+α||β||_1.
#alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.

alpha <- as.numeric(args[2]) #alpha to test in CV

chrom <- args[1]

snpset <- "1KGP_geuMAF0.01"

################################################
### Functions & Libraries
library(glmnet)
library(dplyr)
################################################

##make list of needed population samples
samples <- read.table(gt.dir %&% "samples.txt")
gtsamplelist <- as.character(samples$V1)
if(pop=='ALL'){
  popsamplelist = gtsamplelist
}else{
  popsamplelist <- as.character(dplyr::filter(samples,V2==pop)$V1)
}

##read exp and chr gt dosages
exp <- read.table(exp.dir %&% "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz", header=T)
gt <- read.table(gt.dir %&% "chr" %&% chrom %&% ".dosage.gz")
colnames(gt) <- c("chr", "snp", "pos", "a2freq", "a1", "a2", gtsamplelist)
#a2 is the dose allele

expsamplelist <- colnames(exp)[5:dim(exp)[2]]
gtsamplelist <- colnames(gt)[7:dim(gt)[2]]

#get list of pop samples with gt & exp
newpop <- popsamplelist[popsamplelist %in% gtsamplelist]
popsamplelist <- newpop[newpop %in% expsamplelist]

popgt <- dplyr::select(gt, chr, snp, pos, a1, a2, one_of(popsamplelist))

#pull gene info & expression from pop of interest
popexp <- dplyr::select(exp, TargetID, Chr, Coord, one_of(popsamplelist)) %>% filter(Chr==chrom)
explist <- as.character(popexp$TargetID)

set.seed(42)
groupid <- sample(1:10,length(popsamplelist),replace=TRUE) ##need to use same folds to compare alphas

resultsarray <- array(0,c(length(explist),8))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval")
dimnames(resultsarray)[[2]] <- resultscol
workingbest <- exp.dir %&% "GEUVADIS_" %&% pop %&% "_exp_" %&% k %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% "_" %&% date %&% ".txt"
write(resultscol,file=workingbest,ncolumns=8,sep="\t")

weightcol = c("gene","SNP","refAllele","effectAllele","beta")
workingweight <- en.dir %&% "GEUVADIS_" %&% pop %&% "_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_weights_chr" %&% chrom %&% "_" %&% date %&% ".txt"
write(weightcol,file=workingweight,ncol=5,sep="\t")

for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  gene <- explist[i]
  start <- popexp$Coord[i] - 5e5 ### 500kb TSS lower bound for cis-eQTLS
  end <- popexp$Coord[i] + 5e5 ### 500kb TSS upper bound for cis-eQTLs
  cisgenos <- subset(popgt,popgt[,3]>=start & popgt[,3]<=end) ### pull cis-SNP genotypes
  rownames(cisgenos) <- cisgenos$snp #carry rsids along
  cismat <- as.matrix(cisgenos[,6:dim(cisgenos)[2]]) #get dosages only in matrix format for glmnet
  cismat <- t(cismat) #transpose to match previous code
  expmat <- as.matrix(popexp[,4:dim(popexp)[2]]) #make exp only in matrix format for glmnet
  expmat <- t(expmat) #transpose to match previous code
  colnames(expmat) <- popexp$TargetID #carry gene IDs along
###START HERE
  if(is.null(dim(cismat))){
    bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
  }else{
    minorsnps <- subset(colMeans(cismat), colMeans(cismat,na.rm=TRUE)>0) ###pull snps with at least 1 minor allele###
    minorsnps <- names(minorsnps)
    cismat <- cismat[,minorsnps]
    if(is.null(dim(cismat)) | dim(cismat)[2] == 0){###effectively skips genes with <2 cis-SNPs
      bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else{
      exppheno <- expmat[,gene] ### pull expression data for gene
      exppheno <- scale(exppheno, center=T, scale=T)  ###scale to compare across genes
      exppheno[is.na(exppheno)] <- 0
  
      ##run Cross-Validation over alphalist
      fit <- cv.glmnet(cismat,exppheno,nfolds=k,alpha=alpha,keep=T,foldid=groupid,parallel=F) ##parallel=T is slower on tarbell, not sure why

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
    rsq <- res$r.squared
    pval <- res$coef[2,4]

    resultsarray[gene,] <- c(gene, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval)

    
    ### output best shrunken betas for PrediXcan
    bestbetalist <- names(bestbetas)
    bestbetainfo <- cisgenos[bestbetalist,1:5]
    betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
    betafile<-cbind(gene,betatable[,2],betatable[,4],betatable[,5],betatable[,6]) ##output "gene","SNP","refAllele","effectAllele","beta"
    write(t(betafile),file=workingweight,ncolumns=5,append=T,sep="\t") # t() necessary for correct output from write() function

  }else{
    resultsarray[gene,1] <- gene
    resultsarray[gene,2:8] <- c(NA,NA,NA,NA,0,NA,NA)

  }
  write(resultsarray[gene,],file=workingbest,ncolumns=8,append=T,sep="\t")
}

write.table(resultsarray,file=exp.dir %&% "GEUVADIS_" %&% pop %&% "_exp_" %&% k %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% "_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
