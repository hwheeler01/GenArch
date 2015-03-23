####by Heather E. Wheeler 20141120####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

thresh <- 'fdr0.05' ##threshold for inclusion of trans FHS SNPs

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
lmer.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "gtex-grms/"
grm.trans.dir <- grm.dir %&% "FHS_eQTLs_" %&% thresh %&% "/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein." %&% args[1] ##genes split into 20 files
gencodeset <- args[1]

###############################################
### Scan expression data
idfile <- lmer.dir %&% "ranef.mean0.1.SUBJID"
genefile <- lmer.dir %&% "ranef.mean0.1.GENE"
expfile <- lmer.dir %&% "ranef.mean0.1.SUBJIDxGENE"

subjid <- scan(idfile,"character")
geneid <- scan(genefile, "character")
expdata <- scan(expfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- subjid
colnames(expdata) <- geneid

### Get gene subset to analyze
localfile <- grm.dir %&% "localGRM.list"
locallist <- scan(localfile,"character")

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5]

grmlist <- intersect(locallist,rownames(gencode)) ##genes in genecodeset with grm
ensidlist <- intersect(grmlist,colnames(expdata)) ##genes with grm and exp

### Get individual subset to analyze
grm.id <- read.table(grm.dir %&% "GTEx.global.grm.id")
indidlist <- intersect(subjid,grm.id[,1]) ##subjects with exp and grm
nsubj <- length(indidlist)
tis <- "cross-tissue"

### Get expression data from intersected genes and subjects
localexp <- expdata[indidlist,ensidlist]
localensid <- colnames(localexp)

### Local h2 output matrix

local.mat <- matrix(0,nrow=length(localensid),ncol=7)
colnames(local.mat) <- c("tissue","N","ensid","gene","h2","se","p")

for(i in 1:length(localensid)){
	cat(i,"of",length(localensid),"\n")
	ensid <- localensid[i]
	gene <- as.character(gencode[ensid,6])
	chr <- as.character(gencode[ensid,1])

	#output expression pheno for gcta
	geneexp <- cbind(rownames(localexp),localexp[,i])
	write.table(geneexp, file="tmp.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

	## Y ~ localGRM
	runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
	system(runLOC)
	hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
        res <- c(tis, nsubj, ensid, gene, hsq[14], hsq[15], hsq[25])
	local.mat[i,] <- res
}

write.table(local.mat,file="GTEx.cross-tissue.h2.marginal.local_subset" %&% gencodeset %&% "." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")

########################################################
### Estimate trans h2 for genes with FHS trans eQTLs ###
########################################################

### Get gene subset to analyze
transfile <- grm.dir %&% "transGRM.list"
translist <- scan(transfile,"character")

### Get expression data from intersected transGRM genes and subjects
transoverlap<-intersect(translist,colnames(localexp))
transexp <- localexp[,transoverlap]
transensid <- colnames(transexp)

### Output matrices
loc.mat <- matrix(0,nrow=length(transensid),ncol=7)
colnames(loc.mat) <- c("tissue","N","ensid","gene","local.h2","local.se","local.p")

glo.mat <- matrix(0,nrow=length(transensid),ncol=4)
colnames(glo.mat) <- c("gene","trans.h2","trans.se","trans.p")

locglo.mat <- matrix(0,nrow=length(transensid),ncol=5)
colnames(locglo.mat) <- c("gene","loc.jt.h2","loc.jt.se","trans.jt.h2","trans.jt.se")

for(i in 1:length(transensid)){
        cat(i,"of",length(transensid),"\n")
        ensid <- as.character(transensid[i])
        gene <- as.character(gencode[ensid,6])
        chr <- as.character(gencode[ensid,1])
        c <- as.numeric(substr(chr,4,5))

        #output expression pheno for gcta
        geneexp <- cbind(rownames(transexp),transexp[,ensid])
        write.table(geneexp, file="tmp.pheno." %&% thresh %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

        ## Y ~ localGRM
        runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml --pheno tmp.pheno." %&% thresh %&% gencodeset %&% " --out tmp." %&% thresh %&% gencodeset
        system(runLOC)
        hsq <- scan("tmp." %&% thresh %&% gencodeset %&% ".hsq","character")
        if(hsq[23]=="of"){
                res <- c(tis, nsubj, ensid, gene, NA, NA, NA) ##gcta did not converge, script is reading in previous Y~localGRM+transGRM result
        }else{
              	res <- c(tis, nsubj, ensid, gene, hsq[14], hsq[15], hsq[25])
        }
	loc.mat[i,] <- res

        ## Y ~ transGRM (FHS trans-eQTLs for gene)
        runGLO <- "gcta64 --grm " %&% grm.trans.dir %&% "trans-" %&% ensid %&% " --reml --pheno tmp.pheno." %&% thresh %&% gencodeset %&% " --out tmp." %&% thresh %&% gencodeset
        system(runGLO)
        hsq <- scan("tmp." %&% thresh %&% gencodeset %&% ".hsq","character")
        res <- c(gene, hsq[14], hsq[15], hsq[25])
        glo.mat[i,] <- res

        ## Y ~ localGRM + transGRM (joint model)
        runMULT <- "echo " %&% grm.dir %&% ensid %&% " > tmp.multiGRM." %&% thresh %&% gencodeset
        runMULT2 <- "echo " %&% grm.trans.dir %&% "trans-" %&% ensid %&% " >> tmp.multiGRM." %&% thresh %&% gencodeset
        system(runMULT)
        system(runMULT2)
        runLOCGLO <- "gcta64 --mgrm-bin tmp.multiGRM." %&% thresh %&% gencodeset %&% " --reml --pheno tmp.pheno." %&% thresh %&% gencodeset %&% " --out tmp." %&% thresh %&% gencodeset
        system(runLOCGLO)
        hsq <- scan("tmp." %&% thresh %&% gencodeset %&% ".hsq","character")
        if(hsq[18]=="logL0"){
                res <- c(gene, NA, NA, NA, NA) ##gcta did not converge, script is reading in Y~transGRM result
        }else{
              	res <- c(gene, hsq[17], hsq[18], hsq[20], hsq[21])
        }
	locglo.mat[i,] <- res
}


full.mat <- cbind(loc.mat,glo.mat,locglo.mat)
write.table(full.mat,file=tis %&% ".h2.all.models_FHS" %&% thresh %&% ".subset" %&% gencodeset %&% "_transForGene." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")

