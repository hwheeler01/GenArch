####by Heather E. Wheeler 20141219####
args <- commandArgs(trailingOnly=T)
#args <- c('00','00')
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

thresh <- 'fdr0.05' ##threshold for dir to get local GRMs

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
exp.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/"
annot.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/"
grm.dir <- my.dir %&% "expArch_DGN-WB_imputedGTs/dgn-grms/FHS_eQTLs_" %&% thresh %&% "/"
grm3.dir <- my.dir %&% "expArch_DGN-WB_imputedGTs/dgn-grms/otherChrallSNPs/"

#gencodefile <- annot.dir %&% "gencode.v12.V1.summary.protein.nodup.genenames." %&% args[1] ##genes split into 20 files
#gencodeset <- args[1]
## for subsets that didn't finish in 100 hours, split gencode files into 4 (~250 genes/file)
gencodefile <- annot.dir %&% "gencode.v12.V1.summary.protein.nodup.genenames." %&% args[1] %&% "." %&% args[2]
gencodeset <- args[1] %&% "." %&% args[2]
nperms <- 100 ##number of perms per run, runtime >100 hours for DGN

###############################################
### Scan expression data
idfile <- exp.dir %&% "DGN-WB.exp.ID.list"
genefile <- exp.dir %&% "DGN-WB.exp.GENE.list"
expfile <- exp.dir %&% "DGN-WB.exp.IDxGENE"

subjid <- scan(idfile,"character")
geneid <- scan(genefile, "character")
expdata <- scan(expfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- subjid
colnames(expdata) <- geneid

### Get gene subset to analyze
localfile <- exp.dir %&% "DGN-WB.localGRM.h2.exp.2014-08-30.genelist"
locallist <- scan(localfile,"character")

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,6]

grmlist <- intersect(locallist,rownames(gencode)) ##genes in genecodeset with grm
ensidlist <- intersect(grmlist,colnames(expdata)) ##genes with grm and exp

### Get individual subset to analyze
grm.id <- read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/localGRMs/DGN.global.grm.id")
indidlist <- intersect(subjid,grm.id[,1]) ##subjects with exp and grm
nsubj <- length(indidlist)
tis <- "DGN-WB"

### Get expression data from intersected genes and subjects
localexp <- expdata[indidlist,ensidlist]
localensid <- colnames(localexp)

loc.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms
glo.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms

for(i in 1:length(localensid)){
	cat(i,"of",length(localensid),"\n")
	gene <- localensid[i]
	ensid <- as.character(gencode[gene,5])
	chr <- as.character(gencode[gene,1])

	##observed hsq
	geneexp <- cbind(rownames(localexp),localexp[,i])
	write.table(geneexp, file="tmp.pheno." %&% thresh %&% gencodeset, col.names=F, quote=F) #output pheno for gcta
	
	## Y ~ localGRM + globalGRM (joint model)
        runMULT <- "echo " %&% grm.dir %&% "local-" %&% gene %&% " > tmp.multiGRM." %&% thresh %&% gencodeset
        runMULT2 <- "echo " %&% grm3.dir %&% "global-" %&% gene %&% "-otherChrallSNPs >> tmp.multiGRM." %&% thresh %&% gencodeset
        system(runMULT)
        system(runMULT2)
        runLOCGLO <- "gcta64 --mgrm-bin tmp.multiGRM." %&% thresh %&% gencodeset %&% " --reml-no-constrain --pheno tmp.pheno." %&% thresh %&% gencodeset %&% " --out tmp." %&% thresh %&% gencodeset
        system(runLOCGLO)
        if(file.exists("tmp." %&% thresh %&% gencodeset %&% ".hsq")==TRUE){
                hsq <- scan("tmp." %&% thresh %&% gencodeset %&% ".hsq", "character")
                loch2 <- c(tis, nsubj, ensid, gene, hsq[17])
                gloh2 <- c(tis, nsubj, ensid, gene, hsq[20])
        }else{
                 loch2 <- c(tis, nsubj, ensid, gene, NA)
                 gloh2 <- c(tis, nsubj, ensid, gene, NA)
        }
        loc.mat[i,1:5] <- loch2
        glo.mat[i,1:5] <- gloh2
        system("rm tmp." %&% thresh %&% gencodeset %&% ".hsq")

	

	##permutations
	set.seed(gencodeset)
	for(j in 1:nperms){ ##choose number of permutations
		#output expression pheno for gcta
		geneexp <- cbind(rownames(localexp),localexp[,i])
		perm <- sample(geneexp[,2],replace=F)
		permexp <- cbind(geneexp[,1],perm)
		write.table(permexp, file="tmp.pheno." %&% thresh %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

	      ## Y ~ localGRM + globalGRM (joint model)
	        runMULT <- "echo " %&% grm.dir %&% "local-" %&% gene %&% " > tmp.multiGRM." %&% thresh %&% gencodeset
        	runMULT2 <- "echo " %&% grm3.dir %&% "global-" %&% gene %&% "-otherChrallSNPs >> tmp.multiGRM." %&% thresh %&% gencodeset
	        system(runMULT)
        	system(runMULT2)
        	runLOCGLO <- "gcta64 --mgrm-bin tmp.multiGRM." %&% thresh %&% gencodeset %&% " --reml-no-constrain --pheno tmp.pheno." %&% thresh %&% gencodeset %&% " --out tmp." %&% thresh %&% gencodeset
        	system(runLOCGLO)
        	if(file.exists("tmp." %&% thresh %&% gencodeset %&% ".hsq")==TRUE){
			hsq <- scan("tmp." %&% thresh %&% gencodeset %&% ".hsq", "character")
			loch2 <- hsq[17]
			gloh2 <- hsq[20]
		}else{
			loch2 <- NA
			gloh2 <- NA
		}
		loc.mat[i,j+5] <- loch2
		glo.mat[i,j+5] <- gloh2
		system("rm tmp." %&% thresh %&% gencodeset %&% ".hsq")

	}
}

write.table(loc.mat,file="DGN-WB_h2_joint_local_reml-no-constrain_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(glo.mat,file="DGN-WB_h2_joint_global_reml-no-constrain_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
