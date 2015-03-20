####by Heather E. Wheeler 20141120####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

thresh <- 'fdr0.05' ##threshold for inclusion of trans FHS SNPs
tis <- "DGN-WB"

pre.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
my.dir <- pre.dir %&% "expArch_DGN-WB_imputedGTs/"
annot.dir <- pre.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "dgn-grms/FHS_eQTLs_" %&% thresh %&% "/"
exp.dir <- my.dir %&% "dgn-exp/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% args[1]
gencodeset <- args[1]

###############################################
### Scan expression data
idfile <- exp.dir %&% "DGN-WBexp.ID.list"
genefile <- exp.dir %&% "DGN-WBexp.GENE.list"
expfile <- exp.dir %&% "DGN-WB.rntransform.exp.IDxGENE"

subjid <- scan(idfile,"character")
geneid <- scan(genefile, "character")
expdata <- scan(expfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- subjid
colnames(expdata) <- geneid

### Get gene subset to analyze
system("ls " %&% grm.dir %&% "trans*id > tmp." %&% thresh %&% gencodeset)
system("awk -F \"trans-\" \'{print $2}\' < tmp." %&% thresh %&% gencodeset %&% " | awk -F \".grm\" \'{print $1}\' >transGRM.list." %&% thresh %&% gencodeset)
transfile <- "transGRM.list."	%&% thresh %&% gencodeset
translist <- as.data.frame(scan(transfile,"character")) ##list of trans GRMs
colnames(translist)<-'gene'

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5]
colnames(gencode) <- c('chr','str','start','end','ensid','gene','func','known')

grmlist <- inner_join(translist,gencode,by='gene') ##genes in genecodeset with grm
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
transexp <- expdata[indidlist,ensgene]
transensid <- ensidlist$ensid

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
	geneexp <- cbind(rownames(transexp),transexp[,gene])
	write.table(geneexp, file="tmp.pheno." %&% thresh %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

	## Y ~ localGRM
	runLOC <- "gcta64 --grm " %&% grm.dir %&% "local-" %&% gene %&% " --reml --pheno tmp.pheno." %&% thresh %&% gencodeset %&% " --out tmp." %&% thresh %&% gencodeset
	system(runLOC)
	hsq <- scan("tmp." %&% thresh %&% gencodeset %&% ".hsq","character")
	if(hsq[23]=="of"){
		res <- c(tis, nsubj, ensid, gene, NA, NA, NA) ##gcta did not converge, script is reading in previous Y~localGRM+transGRM result
	}else{
	        res <- c(tis, nsubj, ensid, gene, hsq[14], hsq[15], hsq[25])
        }
	loc.mat[i,] <- res


	## Y ~ transGRM (FHS trans-eQTLs for gene) 
	runGLO <- "gcta64 --grm " %&% grm.dir %&% "trans-" %&% gene %&% " --reml --pheno tmp.pheno." %&% thresh %&% gencodeset %&% " --out tmp." %&% thresh %&% gencodeset
	system(runGLO)
	hsq <- scan("tmp." %&% thresh %&% gencodeset %&% ".hsq","character")
	res <- c(gene, hsq[14], hsq[15], hsq[25])
	glo.mat[i,] <- res 


	## Y ~ localGRM + transGRM (joint model)
	runMULT <- "echo " %&% grm.dir %&% "local-" %&% gene %&% " > tmp.multiGRM." %&% thresh %&% gencodeset
	runMULT2 <- "echo " %&% grm.dir %&% "trans-" %&% gene %&% " >> tmp.multiGRM." %&% thresh %&% gencodeset
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
write.table(full.mat,file=tis %&% ".h2.all.models_FHS" %&% thresh %&% ".Chr" %&% gencodeset %&% "_transForGene." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
