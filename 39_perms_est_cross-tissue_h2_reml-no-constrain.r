####by Heather E. Wheeler 20141120####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
lmer.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "gtex-grms/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein." %&% args[1] ##genes split into 20 files
gencodeset <- args[1]
nperms <- 100 ##number of perms per run, runtime ~26-36 hours

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

loc.only.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms

for(i in 1:length(localensid)){
	cat(i,"of",length(localensid),"\n")
	ensid <- localensid[i]
	gene <- as.character(gencode[ensid,6])
	chr <- as.character(gencode[ensid,1])

	##observed hsq
	geneexp <- cbind(rownames(localexp),localexp[,i])
	write.table(geneexp, file="tmp.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta
	

        ## Y ~ localGRM
        runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml-no-constrain --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
	system(runLOC)
	if(file.exists("tmp." %&% gencodeset %&% ".hsq")==TRUE){
	        hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
		res <- c(tis, nsubj, ensid, gene, hsq[14])
	}else{
                res <- c(tis, nsubj, ensid, gene, NA) ##Error: the information matrix is not invertible.
        }
        loc.only.mat[i,1:5] <- res
	system("rm tmp." %&% gencodeset %&% ".hsq")

	##permutations
	set.seed(gencodeset)
	for(j in 1:nperms){ ##choose number of permutations
		#output expression pheno for gcta
		geneexp <- cbind(rownames(localexp),localexp[,i])
		perm <- sample(geneexp[,2],replace=F)
		permexp <- cbind(geneexp[,1],perm)
		write.table(permexp, file="tmp.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

		## Y ~ localGRM
		runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml-no-constrain --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
		system(runLOC)
		if(file.exists("tmp." %&% gencodeset %&% ".hsq")==TRUE){
                	hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
                	res <- hsq[14]
        	}else{
                	res <- NA ##gcta did not converge               
        	}
		loc.only.mat[i,j+5] <- res
		system("rm tmp." %&% gencodeset %&% ".hsq")
	}
}

write.table(loc.only.mat,file="GTEx_cross-tissue_h2_localonly_reml-no-constrain_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
