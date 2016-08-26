####by Heather E. Wheeler 20141120####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
lmer.dir <- my.dir %&% "lmer.fits/"
annot.dir <- my.dir %&% "gtex-annot/"
grm.dir <- my.dir %&% "gtex-grms/"
exp.dir <- my.dir %&% "gtex-rnaseq/ind-tissues-from-nick/" ##tissue exp adj by 15 PEER, 3PCs, gender

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein." %&% args[1] ##genes split into 20 files
gencodeset <- args[1] ##gene set to run
i <- as.numeric(args[2]) ##place in tissue list to run
nperms <- 100 ##number of perms per run

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

tiswidefile <- exp.dir %&% "GTEx_PrediXmod.tissue.list"
peerlist<- scan(tiswidefile,"character")
squishlist<- gsub(' ','',tislist) ##removes all whitespace to match peerlist & .RDS files

tislist<-intersect(peerlist,squishlist)

tis <- tislist[i]

##read exp data
obs<-readRDS(exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.IDxGENE.RDS")
genefile <- exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.GENE.list"
idfile <- exp.dir %&% "GTEx_PrediXmod." %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.ID.list"
genelist <- scan(genefile,"character")
idlist <- scan(idfile,"character")
colnames(obs) <- genelist
rownames(obs) <- idlist
tissue.exp <- obs[,intersect(colnames(obs),grmlist)] ###pull expression data for chosen tissue & grmlist###

nsubj <- dim(tissue.exp)[1]
### Get expression data from intersected genes and subjects
localexp <- tissue.exp
localensid <- colnames(localexp) #genelist to analyze

##output matrix
loc.only.mat <- matrix(0,nrow=length(localensid),ncol=5+nperms) ##ncol=5 + nperms

for(j in 1:length(localensid)){
	cat(j,"of",length(localensid),"\n")
	ensid <- localensid[j]
	gene <- as.character(gencode[ensid,6])
	chr <- as.character(gencode[ensid,1])

	#output expression pheno for gcta
	geneexp <- cbind(rownames(localexp),localexp[,j])
	write.table(geneexp, file="tmp.TW.pheno." %&% gencodeset %&% "." %&% i, col.names=F, quote=F) #output pheno for gcta

        ## Y ~ localGRM
        runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml-no-constrain --pheno tmp.TW.pheno." %&% gencodeset %&% "." %&% i %&% " --out tmp.TW." %&% gencodeset %&% "." %&% i
        system(runLOC)
	if(file.exists("tmp.TW." %&% gencodeset %&% "." %&% i %&% ".hsq")==TRUE){
	        hsq <- scan("tmp.TW." %&% gencodeset %&% "." %&% i %&% ".hsq","character")
		res <- c(tis, nsubj, ensid, gene, hsq[14])
	}else{
                res <- c(tis, nsubj, ensid, gene, NA) ##Error: the information matrix is not invertible.
        }
        loc.only.mat[j,1:5] <- res
	system("rm tmp.TW." %&% gencodeset %&% "." %&% i %&% ".hsq")
	
	##permutations
        set.seed(gencodeset)
        for(k in 1:nperms){ ##choose number of permutations
		#output expression pheno for gcta
		perm <- sample(geneexp[,2],replace=F)
            	permexp <- cbind(geneexp[,1],perm)
            	write.table(permexp, file="tmp.TW.pheno." %&% gencodeset %&% "." %&% i, col.names=F, quote=F) #output pheno for gcta

            	## Y ~ localGRM
            	runLOC <- "gcta64 --grm " %&% grm.dir %&% ensid %&% " --reml-no-constrain --pheno tmp.TW.pheno." %&% gencodeset %&% "." %&% i %&% " --out tmp.TW." %&% gencodeset %&% "." %&% i
            	system(runLOC)
	    	if(file.exists("tmp.TW." %&% gencodeset %&% "." %&% i %&% ".hsq")==TRUE){
	        	hsq <- scan("tmp.TW." %&% gencodeset %&% "." %&% i %&% ".hsq","character")
		    	res <- hsq[14]
	    	}else{        
                    	res <- NA ##gcta did not converge
            	}
            	loc.only.mat[j,k+5] <- res
	    	system("rm tmp.TW." %&% gencodeset %&% "." %&% i %&% ".hsq")
        }
}

write.table(loc.only.mat,file="GTEx_TW_" %&% tis %&% "_h2_localonly_reml-no-constrain_subset" %&% gencodeset %&% "_" %&% nperms %&% "perms_" %&% date %&% ".txt",quote=F,row.names=F,col.names=F,sep="\t")
