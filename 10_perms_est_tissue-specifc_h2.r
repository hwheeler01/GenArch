####by Heather E. Wheeler 20141120####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

my.dir <- "/group/im-lab/hwheeler/cross-tissue/"
lmer.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "gtex-grms/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein." %&% args[1] ##genes split into 20 files
gencodeset <- args[1] ##gene set to run
i <- as.numeric(args[2]) ##place in tissue list to run
nperms <- 100 ##number of perms per run

###############################################
### Scan expression data
idfile <- lmer.dir %&% "head.resid.SAMPID"
genefile <- lmer.dir %&% "resid.GENE"
expfile <- lmer.dir %&% "head.resid.SAMPIDxGENE"

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

#for(i in 1:length(tislist)){
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

    glo.only.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms
    loc.only.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms
    loc.joint.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms
    glo.joint.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms

    for(j in 1:length(localensid)){
        cat(j,"of",length(localensid),"\n")
	ensid <- localensid[j]
	gene <- as.character(gencode[ensid,6])
	chr <- as.character(gencode[ensid,1])

	#output expression pheno for gcta
	geneexp <- cbind(rownames(localexp),localexp[,j])
	write.table(geneexp, file="tmp.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

	## Y ~ globalGRM
        runGLO <- "gcta64 --grm " %&% grm.dir %&% "GTEx.global" %&% " --reml --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
        system(runGLO)
        hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
        if(hsq[23]=="of"){
                res <- c(tis, nsubj, ensid, gene, NA) ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
        }else{
                res <- c(tis, nsubj, ensid, gene, hsq[14])
        }
        glo.only.mat[j,1:5] <- res


        ## Y ~ localGRM
        runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
        system(runLOC)
        hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
        if(hsq[23]=="of"){
                res <- c(tis, nsubj, ensid, gene, NA) ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
        }else{
                res <- c(tis, nsubj, ensid, gene, hsq[14])
        }
        loc.only.mat[j,1:5] <- res

        ## Y ~ localGRM + globalGRM
        runMULT <- "echo " %&% grm.dir %&% ensid %&% "> tmp.multiGRM." %&% gencodeset
        runMULT2 <- "echo " %&% grm.dir %&% "GTEx.global" %&% ">> tmp.multiGRM." %&% gencodeset
        system(runMULT)
        system(runMULT2)
        runLOCGLO <- "gcta64 --mgrm-bin tmp.multiGRM." %&% gencodeset %&% " --reml --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
        system(runLOCGLO)
        hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
        if(hsq[18]=="logL0"){
                locres <- c(tis, nsubj, ensid, gene, NA) ##gcta did not converge, script is reading in Y~localGRM result
                glores <- c(tis, nsubj, ensid, gene, NA) ##gcta did not converge, script is reading in Y~localGRM result
        }else{
                locres <- c(tis, nsubj, ensid, gene, hsq[17])
                glores <- c(tis, nsubj, ensid, gene, hsq[20])
        }
        loc.joint.mat[j,1:5] <- locres
        glo.joint.mat[j,1:5] <- glores
	
	##permutations
        set.seed(gencodeset)
        for(k in 1:nperms){ ##choose number of permutations
	    #output expression pheno for gcta
            perm <- sample(geneexp[,2],replace=F)
            permexp <- cbind(geneexp[,1],perm)
            write.table(permexp, file="tmp.pheno." %&% gencodeset, col.names=F, quote=F) #output pheno for gcta

            ## Y ~ globalGRM
            runGLO <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
            system(runGLO)
            hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
            if(hsq[23]=="of"){
                    res <- NA ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
            }else{
                    res <- hsq[14]
            }
            glo.only.mat[j,k+5] <- res


            ## Y ~ localGRM
            runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
            system(runLOC)
            hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
            if(hsq[23]=="of"){
                    res <- NA ##gcta did not converge, script is reading in previous Y~localGRM+globalGRM result
            }else{
                    res <- hsq[14]
            }
            loc.only.mat[j,k+5] <- res

            ## Y ~ localGRM + globalGRM
            runMULT <- "echo " %&% grm.dir %&% ensid %&% "> tmp.multiGRM." %&% gencodeset
            runMULT2 <- "echo " %&% grm.dir %&% "GTEx.global" %&% ">> tmp.multiGRM." %&% gencodeset
            system(runMULT)
            system(runMULT2)
            runLOCGLO <- "gcta64 --mgrm-bin tmp.multiGRM." %&% gencodeset %&% " --reml --pheno tmp.pheno." %&% gencodeset %&% " --out tmp." %&% gencodeset
            system(runLOCGLO)
            hsq <- scan("tmp." %&% gencodeset %&% ".hsq","character")
            if(hsq[18]=="logL0"){
                    locres <- NA ##gcta did not converge, script is reading in Y~localGRM result
                    glores <- NA ##gcta did not converge, script is reading in Y~localGRM result
            }else{
                    locres <- hsq[17]
                    glores <- hsq[20]
            }
            loc.joint.mat[j,k+5] <- locres
            glo.joint.mat[j,k+5] <- glores
        }
    }

    write.table(glo.only.mat,file="GTEx_" %&% tis %&% "_h2_globalonly_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
    write.table(loc.only.mat,file="GTEx_" %&% tis %&% "_h2_localonly_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
    write.table(loc.joint.mat,file="GTEx_" %&% tis %&% "_h2_localjoint_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
    write.table(glo.joint.mat,file="GTEx_" %&% tis %&% "_h2_globaljoint_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")

#}
