####by Heather E. Wheeler 20141219####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
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

loc.only.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms

for(i in 1:length(localensid)){
	cat(i,"of",length(localensid),"\n")
	gene <- localensid[i]
	ensid <- as.character(gencode[gene,5])
	chr <- as.character(gencode[gene,1])

	##observed hsq
	geneexp <- cbind(1,localexp[,i])
	write.table(geneexp, file="tmp.dgn.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta
	
        ## Y ~ localGRM
        runLOC <- "gcta64 --grm " %&% grm.dir %&% "DGN." %&% gene %&% " --reml-no-constrain --pheno tmp.dgn.pheno." %&% gencodeset %&% " --out tmp.dgn." %&% gencodeset
	system(runLOC)
	if(file.exists("tmp.dgn." %&% gencodeset %&% ".hsq")==TRUE){
	        hsq <- scan("tmp.dgn." %&% gencodeset %&% ".hsq","character")
		res <- c(tis, nsubj, ensid, gene, hsq[14])
	}else{        
                res <- c(tis, nsubj, ensid, gene, NA) ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
        }
        loc.only.mat[i,1:5] <- res
	system("rm tmp.dgn." %&% gencodeset %&% ".hsq")

	##permutations
	set.seed(gencodeset)
	for(j in 1:nperms){ ##choose number of permutations
		#output expression pheno for gcta
		geneexp <- cbind(1,localexp[,i])
		perm <- sample(geneexp[,2],replace=F)
		permexp <- cbind(geneexp[,1],perm)
		write.table(permexp, file="tmp.dgn.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

		## Y ~ localGRM
		runLOC <- "gcta64 --grm " %&% grm.dir %&% "DGN." %&% gene %&% " --reml-no-constrain --pheno tmp.dgn.pheno." %&% gencodeset %&% " --out tmp.dgn." %&% gencodeset
		system(runLOC)
		if(file.exists("tmp.dgn." %&% gencodeset %&% ".hsq")==TRUE){
			hsq <- scan("tmp.dgn." %&% gencodeset %&% ".hsq","character")
			res <- hsq[14]
		}else{
			res <- NA ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
        	}
		loc.only.mat[i,j+5] <- res
		system("rm tmp.dgn." %&% gencodeset %&% ".hsq")
	}
}

write.table(loc.only.mat,file="DGN-WB_h2_localonly_reml-no-constrain_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
