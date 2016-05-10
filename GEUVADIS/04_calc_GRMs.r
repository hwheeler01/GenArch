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
write.table(popout, file = exp.dir %&% "tmp." %&% pop %&% ".samplelist" %&% chrom, quote=F,col.names = F,row.names = F)

popgt <- data.frame(dplyr::select(gt, chr, snp, pos, a1, a2, one_of(popsamplelist)))

#pull gene info & expression from pop of interest
popexp <- data.frame(dplyr::select(exp, TargetID, Chr, Coord, one_of(popsamplelist)) %>% filter(Chr==chrom))
explist <- as.character(popexp$TargetID)

##make chrGRM

runGCTAchr <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% "chr" %&% chrom %&%  ".mldose.gz " %&% gt.dir %&% "chr" %&% chrom %&%
  ".mlinfo.gz --keep " %&% exp.dir %&% pop %&% ".samplelist --make-grm-bin --out " %&% grm.dir %&% pop %&% ".chr" %&% chrom
system(runGCTAchr)

###make localGRMs, for subset of genes in gencodeset

for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  gene <- explist[i]
  start <- popexp$Coord[i] - 5e5 ### 500kb TSS lower bound for cis-eQTLS
  end <- popexp$Coord[i] + 5e5 ### 500kb TSS upper bound for cis-eQTLs
  cisgenos <- subset(popgt,popgt[,3]>=start & popgt[,3]<=end) ### pull cis-SNP genotypes
  snplist <- cisgenos$snp #carry rsids along   
  write.table(snplist, file= exp.dir %&% "tmp.SNPlist." %&% pop %&% chrom,quote=F,col.names=F,row.names=F)
  runGCTAgrm <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% "chr" %&% chrom %&%  ".mldose.gz " %&% gt.dir %&% "chr" %&% chrom %&% 
    ".mlinfo.gz --keep " %&% exp.dir %&% "tmp." %&% pop %&% ".samplelist" %&% chrom %&% " --make-grm-bin --extract tmp.SNPlist." %&% pop %&% chrom %&% " --out " %&% grm.dir %&% 
    pop %&% "." %&% gene
  system(runGCTAgrm)
}
