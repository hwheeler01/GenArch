####by Heather E. Wheeler 20141120####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

my.dir <- "/group/im-lab/hwheeler/cross-tissue/"
lmer.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "gtex-grms/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein." %&% args[1] ##genes split into 20 files
gencodeset <- args[1]

###############################################
### Scan expression data
idfile <- lmer.dir %&% "resid.SAMPID"
genefile <- lmer.dir %&% "resid.GENE"
expfile <- lmer.dir %&% "resid.SAMPIDxGENE"

sampid <- scan(idfile,"character")
geneid <- scan(genefile, "character")
expdata <- scan(expfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- sampid
colnames(expdata) <- geneid

### Get gene subset to analyze
localfile <- grm.dir %&% "localGRM.list"
locallist <- scan(localfile,"character")

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5]

grmlist <- intersect(locallist,rownames(gencode)) ##genes in genecodeset with grm
ensidlist <- intersect(grmlist,colnames(expdata)) ##genes with grm and exp

gensetexp <- expdata[,ensidlist] ##get exp data from intersected genes

### Get tissue sets to analyze
tissues <- read.table (annot.dir %&% "GTEx_Analysis_2014-06-13.SampleTissue.annot",header=T,sep="\t")
tissuesets <- table(tissues$SMTSD)
largesets <- tissuesets[tissuesets >= 400] ###start with tissues with n>=400
tislist <- names(largesets)

for(i in 1:length(tislist)){
    tis <- tislist[i]
    samplelist <- subset(tissues,SMTSD == tis)
    tissue.exp <- gensetexp[intersect(rownames(gensetexp),samplelist$SAMPID),] ###pull expression data for chosen tissue###
    
    ###take mean of exp for samples with >1 RNA Seq dataset
    subjidtable <- read.table(annot.dir %&% "GTEx_Data_2014-06-13_Annotations_SubjectSampleMappingDS.txt",header=T,sep="\t")
    rownames(subjidtable) <- subjidtable$SAMPID
    subjidlist <- subjidtable[rownames(tissue.exp),1]
    tissue.exp.substr <- tissue.exp
    rownames(tissue.exp.substr) <- subjidlist
    uniqsubjid <- unique(subjidlist)
    for(id in uniqsubjid){  ####take mean of exp for samples with >1 RNA Seq dataset
        matchexp <- tissue.exp.substr[rownames(tissue.exp.substr)==id,]
        if(is.array(matchexp)=='TRUE'){
                expmean <- colMeans(matchexp)
		tissue.exp.substr[as.character(id),] <- expmean
	}
    }

    ### Get subject subset to analyze
    grm.id <- read.table(grm.dir %&% "GTEx.global.grm.id")
    indidlist <- intersect(uniqsubjid,grm.id[,1]) ##subjects with exp and grm
    nsubj <- length(indidlist)

    ### Get expression data from intersected genes and subjects
    localexp <- tissue.exp.substr[indidlist,]
    localensid <- colnames(localexp) #genelist to analyze

    loc.mat <- matrix(0,nrow=length(localensid),ncol=7)
    colnames(loc.mat) <- c("tissue","N","ensid","gene","h2","se","p")

    locglo.mat <- matrix(0,nrow=length(localensid),ncol=6)
    colnames(locglo.mat) <- c("ensid","gene","local.h2","local.se","global.h2","global.se")

    locchrglo.mat <- matrix(0,nrow=length(localensid),ncol=8)
    colnames(locchrglo.mat) <- c("ensid","gene","local.h2","local.se","chr.h2","chr.se","global.h2","global.se")

    for(j in 1:length(localensid)){
	cat(i,"of",length(localensid),"\n")
	ensid <- localensid[j]
	gene <- as.character(gencode[ensid,6])
	chr <- as.character(gencode[ensid,1])

	#output expression pheno for gcta
	geneexp <- cbind(rownames(localexp),localexp[,j])
	write.table(geneexp, file="tmp.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

	## Y ~ localGRM
	runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
	system(runLOC)
	hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
	if(hsq[25]=="V(G3)/Vp"){
                res <- c(tis, nsubj, ensid, gene, NA, NA, NA) ##gcta did not converge, script is reading in previous Y~localGRM+chrGRM+globalGRM result
        }else{
	        res <- c(tis, nsubj, ensid, gene, hsq[14], hsq[15], hsq[25])
        }
	loc.mat[j,] <- res

	## Y ~ localGRM + globalGRM
	runMULT <- "echo " %&% grm.dir %&% ensid %&% "> tmp.multiGRM." %&% gencodeset
	runMULT2 <- "echo " %&% grm.dir %&% "GTEx.global" %&% ">> tmp.multiGRM." %&% gencodeset
	system(runMULT)
        system(runMULT2)
	runLOCGLO <- "gcta64 --mgrm-bin tmp.multiGRM." %&% gencodeset %&% " --reml --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
	system(runLOCGLO)
	hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
	if(hsq[18]=="logL0"){
		res <- c(ensid, gene, NA, NA, NA, NA) ##gcta did not converge, script is reading in Y~localGRM result
	}else{
		res <- c(ensid, gene, hsq[17], hsq[18], hsq[20], hsq[21])
	}
	locglo.mat[j,] <- res

	## Y ~ localGRM + chrGRM + globalGRM
	runMULT <- "echo " %&% grm.dir %&% ensid %&% "> tmp.multiGRM." %&% gencodeset
        runMULT2 <- "echo " %&% grm.dir %&% "GTEx." %&% chr %&% ">> tmp.multiGRM." %&% gencodeset
	runMULT3 <- "echo " %&% grm.dir %&% "GTEx.global" %&% ">> tmp.multiGRM." %&% gencodeset
        system(runMULT)
        system(runMULT2)
	system(runMULT3)     
	runLOCCHRGLO <- "gcta64 --mgrm-bin tmp.multiGRM." %&% gencodeset %&% " --reml --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
        system(runLOCCHRGLO)
        hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
	if(hsq[20]=="LRT"){
                res <- c(ensid, gene, NA, NA, NA, NA, NA, NA) ##gcta did not converge, script is reading in Y~localGRM result
        }else if(hsq[23]=="of"){
		res <- c(ensid, gene, NA, NA, NA, NA, NA, NA) ##gcta did not converge, script is reading in Y~localGRM+globalGRM result
	}else{
	        res <- c(ensid, gene, hsq[20], hsq[21], hsq[23], hsq[24], hsq[26], hsq[27])
        }
	locchrglo.mat[j,] <- res
    }

    full.mat <- cbind(loc.mat,locglo.mat,locchrglo.mat) 
    write.table(full.mat,file="GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_subset" %&% gencodeset %&% "." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
}
