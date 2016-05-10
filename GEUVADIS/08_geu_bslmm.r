####by Heather E. Wheeler 20150108####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('22','A','GBR')
"%&%" = function(a,b) paste(a,b,sep="")
library(dplyr)
library(data.table)
###############################################
### Directories & Variables
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
exp.dir <- my.dir %&% "GEUVADIS/"
annot.dir <- my.dir %&% "gtex-annot/"
gt.dir <- exp.dir %&% "geu_dosage_1KGP_MAF0.01/"
out.dir <- exp.dir %&% "geu_bslmm/"

pop <- args[3]
chrom <- as.numeric(args[1]) 
chrname <- "chr" %&% chrom
whichlist <- args[2]

getquant <- function(x) quantile(x,c(0.5,0.025,0.975)) ##pulls the median and the 95% credible sets

################################################
#make list of needed population samples
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

expsamplelist <- colnames(exp)[5:dim(exp)[2]]
gtsamplelist <- colnames(gt)[7:dim(gt)[2]]

#get list of pop samples with gt & exp
newpop <- popsamplelist[popsamplelist %in% gtsamplelist]
popsamplelist <- newpop[newpop %in% expsamplelist]

popgt <- data.frame(dplyr::select(gt, chr, snp, pos, a1, a2, one_of(popsamplelist)))

#pull gene info & expression from pop of interest
popexp <- data.frame(dplyr::select(exp, TargetID, Chr, Coord, one_of(popsamplelist)) %>% filter(Chr==chrom))
explist <- as.character(popexp$TargetID)


ngroups = 8 ##split chr into ngroups, edit run*.sh files to match groups

if(whichlist == "A"){
  explist <- explist[1:floor(length(explist)/ngroups)]
}else{
  for(i in 2:(ngroups-1)){
    if(whichlist == LETTERS[i]){
      explist <- explist[(floor(length(explist)/ngroups)*(i-1)+1):(floor(length(explist)/ngroups)*i)]
    }
  }
}
if(whichlist == LETTERS[ngroups]){
  explist <- explist[(floor(length(explist)/ngroups)*(ngroups-1)+1):length(explist)]
}


resultsarray <- array(0,c(length(explist),19))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025","pi025",
	"n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975")
dimnames(resultsarray)[[2]] <- resultscol

working100K <- out.dir %&% "working_" %&% pop %&% "_BSLMM-s100K_iterations_" %&% pop %&% "_chr" %&% chrom %&% whichlist %&% 
	"_" %&% date %&% ".txt"
write(resultscol,file=working100K,ncolumns=19,sep="\t")

for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  gene <- explist[i]
  start <- popexp$Coord[i] - 5e5 ### 500kb TSS lower bound for cis-eQTLS
  end <- popexp$Coord[i] + 5e5 ### 500kb TSS upper bound for cis-eQTLs
  cisgenos <- subset(popgt,popgt[,3]>=start & popgt[,3]<=end) ### pull cis-SNP genotypes
  rownames(cisgenos) <- cisgenos$snp #carry rsids along
  expmat <- as.matrix(popexp[,4:dim(popexp)[2]]) #make exp only in matrix format for glmnet
  expmat <- t(expmat) #transpose to match previous code
  colnames(expmat) <- popexp$TargetID #carry gene IDs along

  if(dim(cisgenos)[1] > 0){
    annotfile <- cbind(cisgenos[,2],cisgenos[,3],chrom)
    genofile <- cbind(cisgenos[,2],cisgenos[,4:dim(cisgenos)[2]])
    phenofile <- data.frame(expmat[,gene])

    write.table(annotfile, file=out.dir %&% "tmp.annot." %&% pop %&% chrom %&% whichlist, quote=F, row.names=F, 
	col.names=F, sep=",")
    write.table(genofile, file=out.dir %&% "tmp.geno." %&% pop %&% chrom %&% whichlist, quote=F, row.names=F, 
	col.names=F, sep=",")
    write.table(phenofile, file=out.dir %&% "tmp.pheno." %&% pop %&% chrom %&% whichlist, quote=F, row.names=F,
	col.names=F, sep=",")

    runBSLMM <- "gemma -g " %&% out.dir %&% "tmp.geno." %&% pop %&% chrom %&% whichlist %&% " -p " %&% out.dir %&% 
	"tmp.pheno." %&% pop %&% chrom %&% whichlist %&% " -a " %&% out.dir %&% "tmp.annot." %&% pop %&% chrom %&% 
	whichlist %&% " -bslmm 1 -seed 12345 -s 100000 -o tmp." %&% pop %&% chrom %&%  whichlist
    system(runBSLMM)

    hyp <- read.table(exp.dir %&% "output/tmp." %&% pop %&% chrom %&% whichlist %&% ".hyp.txt",header=T)
    hyp50 <- hyp[(dim(hyp)[1]/2+1):dim(hyp)[1],] #take second half of sampling iterations
    quantres <- apply(hyp50,2,getquant)
    res <- c(gene,quantres[1,],quantres[2,],quantres[3,])

  }else{
    res <- c(gene,rep(NA,18))
  }
  names(res) <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025",
	"pi025","n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975") 
  resultsarray[gene,] <- res
  write(res,file=working100K,ncolumns=19,append=T,sep="\t")
}

write.table(resultsarray,file=out.dir %&% "GEUVADIS_" %&% pop %&% "_BSLMM-s100K_iterations_" %&% pop %&% "_chr" %&% chrom %&% whichlist %&% 
	"_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
