####by Heather E. Wheeler 20150108####
args <- commandArgs(trailingOnly=T)
#args <- '22'
"%&%" = function(a,b) paste(a,b,sep="")

###############################################
### Directories & Variables
#pre <- "/Users/heather/Dropbox/elasticNet_testing"
pre <- ""

my.dir <- pre %&% "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/"
ct.dir <- pre %&% "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/"
gt.dir <- pre %&% "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/DGN-WB_genotypes/by.chr/"
en.dir <- pre %&% "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/EN/hapmap2/transcriptome-DGN-WB/"

k <- 10 ### k-fold CV
n <- 3 #number of k-fold CV replicates
tis <- "DGN-WB"  
chrom <- as.numeric(args[1]) 
chrname <- "chr" %&% chrom

##alpha = The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as
#(1-α)/2||β||_2^2+α||β||_1.
#alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.

alphalist <- 0:20/20 #vector of alphas to test in CV

################################################
### Functions & Libraries

library(glmnet)

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(x))
lower <- function(x) quantile(x,0.025,na.rm=TRUE)
upper <- function(x) quantile(x,0.975,na.rm=TRUE)

## convenience function to select best lambda over n k-fold cv replicates for linear model by Keston edited by Heather to get predicted values over different alphas
glmnet.select <- function(response, covariates, nrep.set = 10, nfold.set = 10, alpha.set, foldid, ...) {
  require(glmnet)
  fullout <- list()
  for(h in 1:length(alphalist)){
    best.lam.sim = vector()
    best.cvm.sim = vector()
    pred.matrix = matrix(0,nrow=dim(covariates)[1],ncol=nrep.set)
    for (i in 1:nrep.set) {
      glmnet.fit = cv.glmnet(covariates, response, nfolds = nfold.set, alpha = alpha.set[h], foldid = foldid[,i], keep = TRUE) 
      new.df = data.frame(glmnet.fit$cvm, glmnet.fit$lambda, glmnet.fit$glmnet.fit$df, 1:length(glmnet.fit$lambda))
      best.lam = new.df[which.min(new.df[,1]),] # needs to be min or max depending on cv measure (MSE min, AUC max, ...)
      cvm.best = best.lam[,1]
      nrow.max = best.lam[,4]
      best.lam.sim[i] = nrow.max
      best.cvm.sim[i] = cvm.best
      pred.matrix[,i] = glmnet.fit$fit.preval[,nrow.max]
    }
    cvm.avg = mean(best.cvm.sim) # average cvm
    nrow.max = as.integer(mean(best.lam.sim)) # best lambda over cv bootstraps
    ret <- as.data.frame(glmnet.fit$glmnet.fit$beta[,nrow.max])
    ret[ret == 0.0] <- NA
    ret.vec = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
    names(ret.vec) = rownames(ret)[which(!is.na(ret))]
    min.lambda <- glmnet.fit$glmnet.fit$lambda[nrow.max]
    pred.avg <- rowMeans(pred.matrix) #when nrep.set>1
#   pred.avg <- pred.matrix #when nrep.set=1
    output = list(ret.vec, cvm.avg, nrow.max, min.lambda, pred.avg, alpha.set[h])
    fullout <- c(fullout,output)
  }
  return(fullout)
}


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

gencodefile <- my.dir %&% 'gencode.v12.V1.summary.protein.nodup.genenames'
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

grouplist <- read.table(ct.dir %&% tis %&% '.10reps.10fold.group.list',header=T)
rownames(grouplist) <- grouplist[,1]
groupid <- grouplist[,2:dim(grouplist)[2]]

resultsarray <- array(0,c(length(explist),8))
dimnames(resultsarray)[[1]] <- explist
dimnames(resultsarray)[[2]] <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval")

allR2array <- array(0,c(length(explist),1,length(alphalist)))
dimnames(allR2array)[[1]] <- explist
dimnames(allR2array)[[2]] <- c("R2")
dimnames(allR2array)[[3]] <- alphalist

allParray <- array(0,c(length(explist),length(alphalist)))
dimnames(allParray)[[1]] <- explist
dimnames(allParray)[[2]] <- alphalist

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
    cisgenos <- scale(cisgenos, center=T, scale=T)
    cisgenos[is.na(cisgenos)] <- 0
    if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs
      bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else{

      exppheno <- exp.w.geno[,gene] ### pull expression data for gene
      exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
      exppheno[is.na(exppheno)] <- 0
      rownames(exppheno) <- rownames(exp.w.geno)
      ##run Cross-Validation over alphalist
      cv <- glmnet.select(exppheno,cisgenos,nrep.set=n,nfold.set=k,alpha.set=alphalist,foldid=groupid) ###run lasso k-fold CV n times to determine best lambda & betas

      allbetas <- list() ##non-zero betas for each alpha 1:length(alpha.set)
      allcvm.avg <- vector() ##mean minimum cross-validated MSE for each alpha
      allnrow.max <- vector() ##mean (as integer) of best lambda's vector position for each alpha
      alllambdas <- vector() ##mean best lambda for each alpha
      pred.avg.mat <- matrix(NA,nrow=length(exppheno),ncol=length(alphalist)) ##predicted values at each alpha
        
      for(j in 1:length(alphalist)*6){
        allbetas <- c(allbetas,cv[j-5])
        allcvm.avg <- c(allcvm.avg,cv[[j-4]])
        allnrow.max <- c(allnrow.max,cv[[j-3]])
        alllambdas <- c(alllambdas,cv[[j-2]])
        pred.avg.mat[,j/6] <- cv[[j-1]]
      }
          
      indexbestbetas <- which.min(allcvm.avg)
      bestbetas <- allbetas[[indexbestbetas]] ###how many SNPs in best predictor?
    }
  }
  if(length(bestbetas) > 0){
    for(a in 1:length(alphalist)){        
      pred.en <- pred.avg.mat[,a] ###mean k-fold CV predictions from n reps
      cvm <- allcvm.avg[a]
      ### calculate correlation between predicted and observed expression
      res <- summary(lm(exppheno~pred.en))
      genename <- as.character(gencode[gene,6])
      allR2array[gene,1,a] <- res$r.squared
      allParray[gene,a] <- res$coef[2,4]
    }  
    idxR2 <- which.max(allR2array[gene,1,]) ##determine alpha that gives max R2
    bestbetas <- allbetas[[idxR2]] ##may differ from min cvm betas
    
    ##for the best alpha, find output
    resultsarray[gene,1] <- genename
    resultsarray[gene,2] <- alphalist[idxR2]
    resultsarray[gene,3] <- allcvm.avg[idxR2] ###add mean minimum cvm (cross-validated mean-squared error) to results
    resultsarray[gene,4] <- allnrow.max[idxR2] ###add mean of best lambda iteration to results
    resultsarray[gene,5] <- alllambdas[idxR2] ###add best lambda to results
    resultsarray[gene,6] <- length(bestbetas) ###add #snps in prediction to results
    resultsarray[gene,7] <- allR2array[gene,1,idxR2] ###lm R2
    resultsarray[gene,8] <- allParray[gene,idxR2] ###lm p-value

    ### output best shrunken betas for PrediXcan
    bestbetalist <- names(bestbetas)
    bestbetainfo <- bim[bestbetalist,]
    betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
    betafile<-cbind(betatable[,2],betatable[,5],betatable[,7]) ###middle column: [,6] for GEUVADIS, [,5] for GTEx/other plink bed/bim/bam files
    colnames(betafile) <- c("SNP","eff.allele","beta")
    rownames(betafile) <- bestbetalist
	  write.table(betafile, file=en.dir %&% gene %&% "-" %&% tis %&% ".txt",quote=F,row.names=F,sep="\t")

  }else{
	  genename <- as.character(gencode[gene,6])
    resultsarray[gene,1] <- genename
    resultsarray[gene,2:8] <- c(NA,NA,NA,NA,0,NA,NA)
  }
}


date <- Sys.Date()
write.table(resultsarray,file=tis %&% "_exp_" %&% k %&% "-foldCV_" %&% n %&% "-reps_elasticNet_bestAlpha_hapmap2snps_predictionInfo_chr" %&% chrom %&% "_" %&% date %&% ".txt",quote=F,row.names=F)
write.table(allR2array, file=tis %&% "_exp_" %&% k %&% "-foldCV_" %&% n %&% "-reps_elasticNet_eachAlphaR2_hapmap2snps_chr" %&% chrom %&% "_" %&% date %&% ".txt",quote=F)
