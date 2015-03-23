####by Heather E. Wheeler 20141120####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

pasteid = function(list,x) paste(list[[x]][1:2],collapse='-') 

thresh <- 'fdr0.05' ##threshold for inclusion of trans FHS SNPs

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
lmer.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "gtex-grms/"
grm.trans.dir <- grm.dir %&% "FHS_eQTLs_" %&% thresh %&% "/"
exp.dir <- my.dir %&% "gtex-rnaseq/ind-tissues-from-nick/" ##tissue exp adj by 15 PEER, 3PCs, gender

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein." %&% args[1] ##genes split into 20 files
gencodeset <- args[1]

###############################################

### Get gene subset to analyze
localfile <- grm.dir %&% "localGRM.list"
locallist <- scan(localfile,"character")

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5]

grmlist <- intersect(locallist,rownames(gencode)) ##genes in genecodeset with grm

### Get tissue sets to analyze
tissues <- read.table (annot.dir %&% "GTEx_Analysis_2014-06-13.SampleTissue.annot",header=T,sep="\t")
tissuesets <- table(tissues$SMTSD)
largesets <- tissuesets[tissuesets >= 100] ###start with tissues with n>=100
tislist <- names(largesets)

for(i in 1:length(tislist)){
    tis <- tislist[i]
    samplelist <- subset(tissues,SMTSD == tis)

    squishtis <- gsub(' ','',tis) ##removes all whitespace to match .RDS files

    ##shorten ID to "GTEX-\w+" to match *ID.list files
    samlis <- as.character(samplelist[,1])
    splitsamlis <- strsplit(samlis,'-')
    shortsample <- character()
    for(h in 1:length(samlis)){
	shortsamlis<-as.character(pasteid(splitsamlis,h))
	shortsample <- c(shortsample,shortsamlis)
    }

    ##read exp data
    obs<-readRDS(exp.dir %&% "GTEx_PrediXmod." %&% squishtis %&% ".exp.adj.15PEERfactors.3PCs.gender.IDxGENE.RDS")
    genefile <- exp.dir %&% "GTEx_PrediXmod." %&% squishtis %&% ".exp.adj.15PEERfactors.3PCs.gender.GENE.list"
    idfile <- exp.dir %&% "GTEx_PrediXmod." %&% squishtis %&% ".exp.adj.15PEERfactors.3PCs.gender.ID.list"
    genelist <- scan(genefile,"character")
    idlist <- scan(idfile,"character")
    colnames(obs) <- genelist
    rownames(obs) <- idlist    

    tissue.exp <- obs[intersect(rownames(obs),shortsample),intersect(colnames(obs),grmlist)] ###pull expression data for chosen tissue & grmlist###
    
    nsubj <- dim(tissue.exp)[1]

    ### Get expression data from intersected genes and subjects
    localexp <- tissue.exp
    localensid <- colnames(localexp) #genelist to analyze

    local.mat <- matrix(0,nrow=length(localensid),ncol=7)
    colnames(local.mat) <- c("tissue","N","ensid","gene","h2","se","p")

    for(j in 1:length(localensid)){
	cat(i,"of",length(localensid),"\n")
	ensid <- localensid[j]
	gene <- as.character(gencode[ensid,6])
	chr <- as.character(gencode[ensid,1])

	#output expression pheno for gcta
	geneexp <- cbind(rownames(localexp),localexp[,j])
	write.table(geneexp, file="tmp.tis-wide.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

	## Y ~ localGRM
	runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml --pheno tmp.tis-wide.pheno." %&% gencodeset %&% " --out tmp.tis-wide." %&% gencodeset
	system(runLOC)
	hsq <- scan("tmp.tis-wide." %&% gencodeset %&% ".hsq","character")
        res <- c(tis, nsubj, ensid, gene, hsq[14], hsq[15], hsq[25])
	local.mat[j,] <- res
    }
    write.table(local.mat,file="GTEx.resid.tissue-wide.h2_" %&% tis %&% "_marginal.local_subset" %&% gencodeset %&% "." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")

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

    for(j in 1:length(transensid)){
        cat(j,"of",length(transensid),"\n")
        ensid <- as.character(transensid[j])
        gene <- as.character(gencode[ensid,6])
        chr <- as.character(gencode[ensid,1])
        c <- as.numeric(substr(chr,4,5))

        #output expression pheno for gcta
        geneexp <- cbind(rownames(transexp),transexp[,ensid])
        write.table(geneexp, file="tmp.tis-wide.pheno." %&% thresh %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

        ## Y ~ localGRM
        runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml --pheno tmp.tis-wide.pheno." %&% thresh %&% gencodeset %&% " --out tmp.tis-wide." %&% thresh %&% gencodeset
        system(runLOC)
        hsq <- scan("tmp.tis-wide." %&% thresh %&% gencodeset %&% ".hsq","character")
        if(hsq[23]=="of"){
                res <- c(tis, nsubj, ensid, gene, NA, NA, NA) ##gcta did not converge, script is reading in previous Y~localGRM+transGRM result
        }else{
              	res <- c(tis, nsubj, ensid, gene, hsq[14], hsq[15], hsq[25])
        }
	loc.mat[j,] <- res

        ## Y ~ transGRM (FHS trans-eQTLs for gene)
        runGLO <- "gcta64 --grm " %&% grm.trans.dir %&% "trans-" %&% ensid %&% " --reml --pheno tmp.tis-wide.pheno." %&% thresh %&% gencodeset %&% " --out tmp.tis-wide." %&% thresh %&% gencodeset
        system(runGLO)
        hsq <- scan("tmp.tis-wide." %&% thresh %&% gencodeset %&% ".hsq","character")
        res <- c(gene, hsq[14], hsq[15], hsq[25])
        glo.mat[j,] <- res

        ## Y ~ localGRM + transGRM (joint model)
        runMULT <- "echo " %&% grm.dir %&% ensid %&% " > tmp.tis-wide.multiGRM." %&% thresh %&% gencodeset
        runMULT2 <- "echo " %&% grm.trans.dir %&% "trans-" %&% ensid %&% " >> tmp.tis-wide.multiGRM." %&% thresh %&% gencodeset
        system(runMULT)
        system(runMULT2)
        runLOCGLO <- "gcta64 --mgrm-bin tmp.tis-wide.multiGRM." %&% thresh %&% gencodeset %&% " --reml --pheno tmp.tis-wide.pheno." %&% thresh %&% gencodeset %&% " --out tmp.tis-wide." %&% thresh %&% gencodeset
        system(runLOCGLO)
        hsq <- scan("tmp.tis-wide." %&% thresh %&% gencodeset %&% ".hsq","character")
        if(hsq[18]=="logL0"){
                res <- c(gene, NA, NA, NA, NA) ##gcta did not converge, script is reading in Y~transGRM result
        }else{
              	res <- c(gene, hsq[17], hsq[18], hsq[20], hsq[21])
        }
	locglo.mat[j,] <- res
}


full.mat <- cbind(loc.mat,glo.mat,locglo.mat)
write.table(full.mat,file="GTEx.resid.tissue-wide.h2_" %&% tis %&% "_FHS" %&% thresh %&% ".subset" %&% gencodeset %&% "_transForGene." %&% date %&% ".txt",quote=F,row.names=F,sep="\t")


}


