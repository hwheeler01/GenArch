####by Heather E. Wheeler 20160421####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('22','GBR')
"%&%" = function(a,b) paste(a,b,sep="")

###############################################
### Directories & Variables

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
exp.dir <- my.dir %&% "GEUVADIS/"
annot.dir <- my.dir %&% "gtex-annot/"
gt.dir <- exp.dir %&% "geu_dosage_1KGP_MAF0.01/"
grm.dir <- exp.dir %&% "geu_grms/"

chrom <- args[1]
pop <- args[2] 

################################################
### Functions & Libraries
library(glmnet)
library(dplyr)
library(data.table)
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

expfile <- "gunzip -c " %&% exp.dir %&% "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"
exp <- fread(expfile) #faster than read.table, requires string of above command for gz files

gtfile <- "gunzip -c " %&% gt.dir %&% "chr" %&% chrom %&% ".dosage.gz"
gt <- fread(gtfile) #faster than read.table, requires string of above command for gz files
setnames(gt,c("chr", "snp", "pos", "a2freq", "a1", "a2", gtsamplelist))
#a2 is the dose allele

expsamplelist <- colnames(exp)[5:dim(exp)[2]]
gtsamplelist <- colnames(gt)[7:dim(gt)[2]]

#get list of pop samples with gt & exp
newpop <- popsamplelist[popsamplelist %in% gtsamplelist]
popsamplelist <- newpop[newpop %in% expsamplelist]
popout <- cbind(popsamplelist,popsamplelist)

popgt <- data.frame(dplyr::select(gt, chr, snp, pos, a1, a2, one_of(popsamplelist)))

#pull gene info & expression from pop of interest
popexp <- data.frame(dplyr::select(exp, TargetID, Chr, Coord, one_of(popsamplelist)) %>% filter(Chr==chrom))
explist <- as.character(popexp$TargetID)
expmat <- as.matrix(popexp[,4:dim(popexp)[2]]) #make exp only in matrix format
expmat <- t(expmat) #transpose to match previous code
colnames(expmat) <- popexp$TargetID #carry gene IDs along

loc.mat <- matrix(0,nrow=length(explist),ncol=6)
colnames(loc.mat) <- c("pop","N","ensid","local.h2","local.se","local.p")

for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  ensid <- explist[i]
  exppheno <- expmat[,ensid] 
  geneexp <- cbind(popout, exppheno)
  nsubj <- dim(geneexp)[1]
  write.table(geneexp, file="tmp.pheno." %&% pop %&% chrom, row.names=F, col.names=F, quote=F) #output pheno for gcta
  
  runLOC <- "gcta64 --grm " %&% grm.dir %&% pop %&% "." %&% ensid %&% " --reml --pheno tmp.pheno." %&% pop %&% chrom %&% " --out tmp." %&% 
    pop %&% chrom
  system(runLOC)
  if(file.exists("tmp." %&% pop %&% chrom %&% ".hsq")==TRUE){
    hsq <- scan("tmp." %&% pop %&% chrom %&% ".hsq","character")
    res <- c(pop, nsubj, ensid, hsq[14], hsq[15], hsq[25])
  }else{ #gcta did not coverge or Error: the information matrix is not invertible.
    res <- c(pop, nsubj, ensid, NA, NA, NA)
  }
  loc.mat[i,] <- res
  system("rm tmp." %&% pop %&% chrom %&% ".hsq")
}

write.table(loc.mat,file=exp.dir %&% "GEUVADIS_" %&% pop %&% "_exp_gcta_reml_local_h2_1KGP_geuMAF0.01_chr" %&% chrom %&% "_" %&% date %&% 
              ".txt",quote=F,row.names=F,sep="\t")
