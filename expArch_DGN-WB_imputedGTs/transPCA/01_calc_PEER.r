#before R call run: R CMD INSTALL ~/R_peer_source_1.3.tgz, see my.dir %&% "run_scripts/run_02_calc_PEER.sh"
####by Heather E. Wheeler 2017-03-06####
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
args <- commandArgs(trailingOnly=T)

###############################################
### Directories & Variables
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
rna.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/"
trans.dir <- my.dir %&% "expArch_DGN-WB_imputedGTs/transPCA/"
tis <- "DGN-WB"

Nk <- args[1] ##number of peer factors to calculate, recommend 25% of sample size, but no more than 100, GTEx included 15 in pilot analyses

################################################
### Functions & Libraries

library(peer)
library(preprocessCore)
source(my.dir %&% 'GenABEL/R/ztransform.R')
source(my.dir %&% 'GenABEL/R/rntransform.R')

################################################
###Scan data
rpkmid <- rna.dir %&% tis %&% ".exp.ID.list"
expid <- scan(rpkmid,"character")
rpkmgene <- rna.dir %&% tis %&% ".exp.GENE.list"
geneid <- scan(rpkmgene,"character")
rpkmfile <- rna.dir %&% tis %&% ".exp.IDxGENE"
expdata <- scan(rpkmfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- expid
colnames(expdata) <- geneid

###quantile normalize and transform to standard normal expdata matrix, as in GTEx paper###
gex <- t(expdata) #get genes in cols
qn.expdata <- normalize.quantiles(gex) ##quantile normalize
rn.qn.expdata <- apply(qn.expdata,1,"rntransform") ##rank transform to normality & transposes##

###Now we can create the model object, ### from https://github.com/PMBio/peer/wiki/Tutorial

model = PEER()

###set the observed data,

PEER_setPhenoMean(model,as.matrix(rn.qn.expdata))

dim(PEER_getPhenoMean(model))

###(NULL response means no error here), say we want to infer K=20 hidden confounders,

PEER_setNk(model,Nk)

PEER_getNk(model)

####and perform the inference. ###for Nk=20 and DGN-WB, it took 29 iterations

PEER_update(model)

factors = PEER_getX(model)

###make outputfile
rownames(factors) <- colnames(gex)
outdf <- data.frame(colnames(gex), colnames(gex)) #FID IID
outdf <- cbind(outdf, factors)
####make header in preparation for gcta
headerstring <- "FID IID "
for(i in 1:Nk){
    headerstring <- headerstring %&% "PF" %&% i %&% " "
}
colnames(outdf) <- strsplit(headerstring, " ")[[1]]
write.table(outdf,file=trans.dir %&% tis %&% "." %&% Nk %&% ".PEER.factors." %&% date %&% ".txt", quote=F, row.names=F, col.names=F) #no header for gcta

weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

pdf(file=trans.dir %&% tis %&% "." %&% Nk %&% ".PEER.factors.plotmodel." %&% date %&% ".pdf")
PEER_plotModel(model)
dev.off()
