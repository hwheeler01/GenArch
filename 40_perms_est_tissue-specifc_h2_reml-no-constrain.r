####by Heather E. Wheeler 20141120####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
lmer.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "gtex-grms/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein." %&% args[1] ##genes split into 20 files
gencodeset <- args[1] ##gene set to run
i <- as.numeric(args[2]) ##place in tissue list to run
nperms <- 100 ##number of perms per run

###############################################
### Scan expression data
idfile <- lmer.dir %&% "resid.mean0.1.SAMPID"
genefile <- lmer.dir %&% "resid.mean0.1.GENE"
expfile <- lmer.dir %&% "resid.mean0.1.SAMPIDxGENE"

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
tislist <- scan(my.dir %&% "nine.spaces.tissue.list","c",sep="\n")

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

    loc.only.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms

    for(j in 1:length(localensid)){
        cat(j,"of",length(localensid),"\n")
	ensid <- localensid[j]
	gene <- as.character(gencode[ensid,6])
	chr <- as.character(gencode[ensid,1])

	#output expression pheno for gcta
	geneexp <- cbind(rownames(localexp),localexp[,j])
	write.table(geneexp, file="tmp.pheno." %&% gencodeset %&% "." %&% i, col.names=F, quote=F) #output pheno for gcta

        ## Y ~ localGRM
        runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml-no-constrain --pheno tmp.pheno." %&% gencodeset %&% "." %&% i %&% " --out tmp." %&% gencodeset %&% "." %&% i
        system(runLOC)
	if(file.exists("tmp." %&% gencodeset %&% "." %&% i %&% ".hsq")==TRUE){
	        hsq <- scan("tmp." %&% gencodeset %&% "." %&% i %&% ".hsq","character")
		res <- c(tis, nsubj, ensid, gene, hsq[14])
	}else{
                res <- c(tis, nsubj, ensid, gene, NA) ##Error: the information matrix is not invertible.
        }
        loc.only.mat[j,1:5] <- res
	system("rm tmp." %&% gencodeset %&% "." %&% i %&% ".hsq")
	
	##permutations
        set.seed(gencodeset)
        for(k in 1:nperms){ ##choose number of permutations
	    #output expression pheno for gcta
            perm <- sample(geneexp[,2],replace=F)
            permexp <- cbind(geneexp[,1],perm)
            write.table(permexp, file="tmp.pheno." %&% gencodeset %&% "." %&% i, col.names=F, quote=F) #output pheno for gcta

            ## Y ~ localGRM
            runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml-no-constrain --pheno tmp.pheno." %&% gencodeset %&% "." %&% i %&% " --out tmp." %&% gencodeset %&% "." %&% i
            system(runLOC)
	    if(file.exists("tmp." %&% gencodeset %&% "." %&% i %&% ".hsq")==TRUE){
	            hsq <- scan("tmp." %&% gencodeset %&% "." %&% i %&% ".hsq","character")
		    res <- hsq[14]
	    }else{        
                    res <- NA ##gcta did not converge
            }
            loc.only.mat[j,k+5] <- res
	    system("rm tmp." %&% gencodeset %&% "." %&% i %&% ".hsq")
        }
    }

    write.table(loc.only.mat,file="GTEx_TS_" %&% tis %&% "_h2_localonly_reml-no-constrain_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")

#}
