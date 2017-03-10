####by Heather E. Wheeler 20170308####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)
library(preprocessCore)

thresh <- 'fdr0.05' ##threshold for inclusion of trans FHS SNPs
tis <- "DGN-WB"

pre.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
my.dir <- pre.dir %&% "expArch_DGN-WB_imputedGTs/"
annot.dir <- pre.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "dgn-grms/FHS_eQTLs_" %&% thresh %&% "/"
rna.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/"
tpca.dir <- my.dir %&% "transPCA/"

source(pre.dir %&% 'GenABEL/R/ztransform.R')
source(pre.dir %&% 'GenABEL/R/rntransform.R')

gencodeset <- args[1]
gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% gencodeset
Nk <- args[2]

###############################################
### Scan expression data
rpkmid <- rna.dir %&% tis %&% ".exp.ID.list"
subjid <- scan(rpkmid,"character")
rpkmgene <- rna.dir %&% tis %&% ".exp.GENE.list"
geneid <- scan(rpkmgene,"character")
rpkmfile <- rna.dir %&% tis %&% ".exp.IDxGENE"
expdata <- scan(rpkmfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)

###quantile normalize and transform to standard normal expdata matrix, as in GTEx paper###
gex <- t(expdata) #get genes in cols
qn.expdata <- normalize.quantiles(gex) ##quantile normalize
expdata <- apply(qn.expdata,1,"rntransform") ##rank transform to normality & transposes##
rownames(expdata) <- subjid
colnames(expdata) <- geneid

### Get gene subset to analyze
localfile <- my.dir %&% "localGRM.list"
locallist <- as.data.frame(scan(localfile,"character")) ##list of local GRMs
colnames(locallist)<-'gene'

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5]
colnames(gencode) <- c('chr','str','start','end','ensid','gene','func','known')

grmlist <- inner_join(locallist,gencode,by='gene') ##genes in genecodeset with grm
geneidlist <- as.data.frame(geneid)
colnames(geneidlist) <- 'gene'
ensidlist <- inner_join(grmlist,geneidlist,by='gene') ##genes with grm and exp
rownames(ensidlist) <- ensidlist$ensid
ensgene <- as.character(ensidlist$gene)

### Get individual subset to analyze
grm.id <- read.table(grm.dir %&% "DGN.global_Chr" %&% gencodeset %&% ".grm.id")
indidlist <- intersect(subjid,grm.id[,1]) ##subjects with exp and grm
nsubj <- length(indidlist)

### Get expression data from intersected genes and subjects
localexp <- expdata[indidlist,ensgene]
localensid <- ensidlist$ensid

### Output matrices
loc.mat <- matrix(0,nrow=length(localensid),ncol=7)
colnames(loc.mat) <- c("tissue","N","ensid","gene","local.h2","local.se","local.p")

glo.mat <- matrix(0,nrow=length(localensid),ncol=4)
colnames(glo.mat) <- c("gene","global.h2","global.se","global.p")

loc.adj.mat <- matrix(0,nrow=length(localensid),ncol=4)
colnames(loc.adj.mat) <- c("gene","tpca.local.h2","tpca.locall.se","tpca.local.p")

glo.adj.mat <- matrix(0,nrow=length(localensid),ncol=4)
colnames(glo.adj.mat) <- c("gene","tpca.global.h2","tpca.global.se","tpca.global.p")

for(i in 1:length(localensid)){
	cat(i,"of",length(localensid),"\n")
	ensid <- as.character(localensid[i])
	gene <- as.character(gencode[ensid,6])
	chr <- as.character(gencode[ensid,1])
	c <- as.numeric(substr(chr,4,5))

	#output expression pheno for gcta
	geneexp <- cbind(rownames(localexp),localexp[,gene])
	write.table(geneexp, file="tmp.pheno." %&% Nk %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

	## Y ~ localGRM + PFs
	runLOC <- "gcta64 --grm " %&% grm.dir %&% "local-" %&% gene %&% " --reml-no-constrain --pheno tmp.pheno." %&% 
	  Nk %&% gencodeset %&% " --out tmp." %&% Nk %&% gencodeset %&% " --qcovar DGN-WB." %&% Nk %&% 
    ".PEER.factors.2017-03-06.txt"
	system(runLOC)
	if(file.exists("tmp." %&% Nk %&% gencodeset %&% ".hsq")==TRUE){
		hsq <- scan("tmp." %&% Nk %&% gencodeset %&% ".hsq","character")
        	res <- c(tis, nsubj, ensid, gene, hsq[14], hsq[15], hsq[25])
	}else{ #gcta did not coverge or Error: the information matrix is not invertible.
		res <- c(tis, nsubj, ensid, gene, NA, NA, NA)
	}
	loc.mat[i,] <- res
	system("rm tmp." %&% Nk %&% gencodeset %&% ".hsq")

	## Y ~ globalGRM (all SNPs) + PFs
	runGLO <- "gcta64 --grm " %&% my.dir %&% "dgn-grms/DGN.global_Chr1-22 --reml-no-constrain --pheno tmp.pheno." %&% 
	  Nk %&% gencodeset %&% " --out tmp." %&% Nk %&% gencodeset %&% " --qcovar DGN-WB." %&% Nk %&% 
	  ".PEER.factors.2017-03-06.txt"
	system(runGLO)
	if(file.exists("tmp." %&% Nk %&% gencodeset %&% ".hsq")==TRUE){
		hsq <- scan("tmp." %&% Nk %&% gencodeset %&% ".hsq","character")
		res <- c(gene, hsq[14], hsq[15], hsq[25])
	}else{
                res <- c(gene, NA, NA, NA) ##gcta did not converge
        }
	glo.mat[i,] <- res 
	system("rm tmp." %&% Nk %&% gencodeset %&% ".hsq")

	## Y ~ localGRM + adjPFs
	runLOC <- "gcta64 --grm " %&% grm.dir %&% "local-" %&% gene %&% " --reml-no-constrain --pheno tmp.pheno." %&% 
	  Nk %&% gencodeset %&% " --out tmp." %&% Nk %&% gencodeset %&% " --qcovar adj_unconstrained_DGN-WB." %&% Nk %&% 
	  ".PEER.factors.2017-03-06.txt"
	system(runLOC)
	if(file.exists("tmp." %&% Nk %&% gencodeset %&% ".hsq")==TRUE){
	  hsq <- scan("tmp." %&% Nk %&% gencodeset %&% ".hsq","character")
	  res <- c(gene, hsq[14], hsq[15], hsq[25])
	}else{ #gcta did not coverge or Error: the information matrix is not invertible.
	  res <- c(gene, NA, NA, NA)
	}
	loc.adj.mat[i,] <- res
	system("rm tmp." %&% Nk %&% gencodeset %&% ".hsq")
	
	## Y ~ globalGRM (all SNPs) + adjPFs
	runGLO <- "gcta64 --grm " %&% my.dir %&% "dgn-grms/DGN.global_Chr1-22 --reml-no-constrain --pheno tmp.pheno." %&% 
	  Nk %&% gencodeset %&% " --out tmp." %&% Nk %&% gencodeset %&% " --qcovar adj_unconstrained_DGN-WB." %&% Nk %&% 
	  ".PEER.factors.2017-03-06.txt"
	system(runGLO)
	if(file.exists("tmp." %&% Nk %&% gencodeset %&% ".hsq")==TRUE){
	  hsq <- scan("tmp." %&% Nk %&% gencodeset %&% ".hsq","character")
	  res <- c(gene, hsq[14], hsq[15], hsq[25])
	}else{
	  res <- c(gene, NA, NA, NA) ##gcta did not converge
	}
	glo.adj.mat[i,] <- res 
	system("rm tmp." %&% Nk %&% gencodeset %&% ".hsq")
}

full.mat <- cbind(loc.mat,loc.adj.mat,glo.mat,glo.adj.mat)
#first run: constrained adj PFs, output:
#write.table(full.mat,file=tis %&% ".h2_transPCA-constrained_" %&% Nk %&% "_PFs_Chr" %&% gencodeset %&% 
#              "_globalAll_reml-no-constrain." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
#second run: unconstrained adj PFs, output:
write.table(full.mat,file=tis %&% ".h2_transPCA-unconstrained_" %&% Nk %&% "_PFs_Chr" %&% gencodeset %&% 
              "_globalAll_reml-no-constrain." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")