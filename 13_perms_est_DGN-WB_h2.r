####by Heather E. Wheeler 20141219####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

my.dir <- "/group/im-lab/hwheeler/cross-tissue/"
exp.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/"
annot.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/"
grm.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/localGRMs/"

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
grm.id <- read.table(grm.dir %&% "DGN.global.grm.id")
indidlist <- intersect(subjid,grm.id[,1]) ##subjects with exp and grm
nsubj <- length(indidlist)
tis <- "DGN-WB"

### Get expression data from intersected genes and subjects
localexp <- expdata[indidlist,ensidlist]
localensid <- colnames(localexp)

glo.only.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms
loc.only.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms
loc.joint.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms
glo.joint.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms


for(i in 1:length(localensid)){
	cat(i,"of",length(localensid),"\n")
	gene <- localensid[i]
	ensid <- as.character(gencode[gene,5])
	chr <- as.character(gencode[gene,1])

	##observed hsq
	geneexp <- cbind(1,localexp[,i])
	write.table(geneexp, file="tmp.dgn.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta
	
	## Y ~ globalGRM
	runGLO <- "gcta64 --grm " %&% grm.dir %&% "DGN.global" %&% " --reml --pheno tmp.dgn.pheno." %&% gencodeset %&% " --out tmp.dgn." %&% gencodeset
	system(runGLO)
	hsq <- scan("tmp.dgn." %&% gencodeset %&% ".hsq","character")
        if(hsq[23]=="of"){
                res <- c(tis, nsubj, ensid, gene, NA) ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
        }else{
              	res <- c(tis, nsubj, ensid, gene, hsq[14])
        }
	glo.only.mat[i,1:5] <- res


        ## Y ~ localGRM
        runLOC <- "gcta64 --grm " %&% grm.dir %&% "DGN." %&% gene %&% " --reml --pheno tmp.dgn.pheno." %&% gencodeset %&% " --out tmp.dgn." %&% gencodeset
	system(runLOC)
        hsq <- scan("tmp.dgn." %&% gencodeset %&% ".hsq","character")
        if(hsq[23]=="of"){
                res <- c(tis, nsubj, ensid, gene, NA) ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
        }else{
                res <- c(tis, nsubj, ensid, gene, hsq[14])
        }
        loc.only.mat[i,1:5] <- res

        ## Y ~ localGRM + globalGRM
        runMULT <- "echo " %&% grm.dir %&% "DGN." %&% gene %&% "> tmp.dgn.multiGRM." %&% gencodeset
        runMULT2 <- "echo " %&% grm.dir %&% "DGN.global" %&% ">> tmp.dgn.multiGRM." %&% gencodeset
        system(runMULT)
        system(runMULT2)
	runLOCGLO <- "gcta64 --mgrm-bin tmp.dgn.multiGRM." %&% gencodeset %&% " --reml --pheno tmp.dgn.pheno." %&% gencodeset %&% " --out tmp.dgn." %&% gencodeset
        system(runLOCGLO)
        hsq <- scan("tmp.dgn." %&% gencodeset %&% ".hsq","character")
        if(hsq[18]=="logL0"){
                locres <- c(tis, nsubj, ensid, gene, NA) ##gcta did not converge, script is reading in Y~localGRM result
		glores <- c(tis, nsubj, ensid, gene, NA) ##gcta did not converge, script is reading in Y~localGRM result
        }else{
		locres <- c(tis, nsubj, ensid, gene, hsq[17])
		glores <- c(tis, nsubj, ensid, gene, hsq[20])
        }
        loc.joint.mat[i,1:5] <- locres
	glo.joint.mat[i,1:5] <- glores

	##permutations
	set.seed(gencodeset)
	for(j in 1:nperms){ ##choose number of permutations
		#output expression pheno for gcta
		geneexp <- cbind(1,localexp[,i])
		perm <- sample(geneexp[,2],replace=F)
		permexp <- cbind(geneexp[,1],perm)
		write.table(permexp, file="tmp.dgn.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

		## Y ~ globalGRM
                runGLO <- "gcta64 --grm " %&% grm.dir %&% "DGN.global" %&% " --reml --pheno tmp.dgn.pheno." %&% gencodeset %&% " --out tmp.dgn." %&% gencodeset
                system(runGLO)
                hsq <- scan("tmp.dgn." %&% gencodeset %&% ".hsq","character")
                if(hsq[23]=="of"){
                        res <- NA ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
                }else{
                        res <- hsq[14]
                }
                glo.only.mat[i,j+5] <- res


		## Y ~ localGRM
		runLOC <- "gcta64 --grm " %&% grm.dir %&% "DGN." %&% gene %&% " --reml --pheno tmp.dgn.pheno." %&% gencodeset %&% " --out tmp.dgn." %&% gencodeset
		system(runLOC)
		hsq <- scan("tmp.dgn." %&% gencodeset %&% ".hsq","character")
		if(hsq[23]=="of"){
			res <- NA ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
		}else{
	        	res <- hsq[14]
        	}
		loc.only.mat[i,j+5] <- res

		## Y ~ localGRM + globalGRM
		runMULT <- "echo " %&% grm.dir %&% "DGN." %&% gene %&% "> tmp.dgn.multiGRM." %&% gencodeset
		runMULT2 <- "echo " %&% grm.dir %&% "DGN.global" %&% ">> tmp.dgn.multiGRM." %&% gencodeset
		system(runMULT)
        	system(runMULT2)
		runLOCGLO <- "gcta64 --mgrm-bin tmp.dgn.multiGRM." %&% gencodeset %&% " --reml --pheno tmp.dgn.pheno." %&% gencodeset %&% " --out tmp.dgn." %&% gencodeset
		system(runLOCGLO)
		hsq <- scan("tmp.dgn." %&% gencodeset %&% ".hsq","character")
		if(hsq[18]=="logL0"){
			locres <- NA ##gcta did not converge, script is reading in Y~localGRM result
			glores <- NA ##gcta did not converge, script is reading in Y~localGRM result
		}else{
			locres <- hsq[17]
			glores <- hsq[20]
		}
		loc.joint.mat[i,j+5] <- locres
		glo.joint.mat[i,j+5] <- glores
	}
}

write.table(glo.only.mat,file="DGN-WB_h2_globalonly_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(loc.only.mat,file="DGN-WB_h2_localonly_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(loc.joint.mat,file="DGN-WB_h2_localjoint_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(glo.joint.mat,file="DGN-WB_h2_globaljoint_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
