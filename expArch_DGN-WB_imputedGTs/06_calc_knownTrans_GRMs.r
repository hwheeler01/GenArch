####by Heather E. Wheeler 20141114####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

###make DGN transGRMs using known trans-eQTLs from FHS at a given threshold (unique to gene: trans-GENE.grm)

pre.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
my.dir <- pre.dir %&% "expArch_DGN-WB_imputedGTs/"
rna.dir <- my.dir %&% "dgn-exp/"
annot.dir <- pre.dir %&% "gtex-annot/"
gt.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-imputation/DGN-imputed-for-PrediXcan/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein." %&% args[1]
gencodeset <- args[1]

thresh <- 'fdr0.05'
eqtlfile <- my.dir %&% "Framingham_trans-eqtl-gene_" %&% thresh %&% "_hapmapSnpsCEU.list"
eqtllist <- read.table(eqtlfile)

grm.dir <- my.dir %&% "dgn-grms/FHS_eQTLs_" %&% thresh %&% "/"

machpre <- "DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.FHSfdr0.05_trans-eQTLs."

###make transGRMs using known trans-eQTLs from FHS 

gencode <- read.table(gencodefile)
colnames(gencode) <- c('chr','str','start','end','ensid','gene','func','known')
rownames(gencode) <- gencode[,5]
ensidlist <- gencode[,5]

#finished.grms <- scan("done.grms." %&% thresh,"character") ###already calculated a bunch of grms in first run (quit at 24hr), don't run them again
#geneidlist <- gencode[,6]
#todolist <- setdiff(geneidlist,finished.grms)

#todo<-data.frame(todolist)
#colnames(todo) <- 'gene'
#todoinfo <- inner_join(gencode,todo,by='gene')
#ensidlist <- todoinfo[,6]

for(i in 1:length(ensidlist)){
    cat(i,"/",length(ensidlist),"\n")
    ensid <- ensidlist[i]
    geneinfo <- gencode[as.character(ensid),]
    gene <- as.character(geneinfo[1,6])
    transsnps <- subset(eqtllist,eqtllist[,2]==gene) ### pull trans-eQTL list
    snplist <- as.vector(transsnps[,1])
    if(length(snplist) > 0){
        write.table(snplist, file= my.dir %&% "tmp." %&% thresh %&% ".SNPlist." %&% gencodeset,quote=F,col.names=F,row.names=F)
        runGCTAgrm <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "chr1-22.mldose.gz " %&% gt.dir %&% machpre %&% "chr1-22.mlinfo.gz --make-grm-bin --extract tmp." %&% thresh %&% ".SNPlist." %&% gencodeset %&% " --out " %&% grm.dir %&% "trans-" %&% gene %&% " --thread-num 8"
        system(runGCTAgrm)
    }
}
