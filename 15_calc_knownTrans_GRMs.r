####by Heather E. Wheeler 20150323####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)

### FHS trans-eQTLs (fdr < 0.05),  non-ambiguous SNPs were extracted from the GTEx genotype data with:
### /group/im-lab/nas40t2/hwheeler/cross-tissue/gtex-genotypes/01_vcf2dosage.mach_gtex_FHStrans.pl

### make GTEx transGRMs using known trans-eQTLs from FHS at a given threshold (unique to gene: trans-GENE.grm):

thresh <- 'fdr0.05'

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
rna.dir <- my.dir %&% "gtex-rnaseq/"
annot.dir <- my.dir %&% "gtex-annot/"
gt.dir <- my.dir %&% "gtex-genotypes/"
grm.dir <- my.dir %&% "gtex-grms/FHS_eQTLs_" %&% thresh %&% "/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein." %&% args[1]
gencodeset <- args[1]

eqtlfile <- my.dir %&% "expArch_DGN-WB_imputedGTs/Framingham_trans-eqtl-gene_" %&% thresh %&% "_hapmapSnpsCEU.list"
eqtllist <- read.table(eqtlfile)

machpre <- "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.FHSfdr0.05_trans-eQTLs.chr1-22." 

###make transGRMs using known trans-eQTLs from FHS

gencode <- read.table(gencodefile)
colnames(gencode) <- c('chr','str','start','end','ensid','gene','func','known')
rownames(gencode) <- gencode[,5]
ensidlist <- gencode[,5]

for(i in 1:length(ensidlist)){
    cat(i,"/",length(ensidlist),"\n")
    ensid <- ensidlist[i]
    geneinfo <- gencode[as.character(ensid),]
    gene <- as.character(geneinfo[1,6])
    transsnps <- subset(eqtllist,eqtllist[,2]==gene) ### pull trans-eQTL list
    snplist <- as.vector(transsnps[,1])
    if(length(snplist) > 0){
        write.table(snplist, file= my.dir %&% "tmp." %&% thresh %&% ".SNPlist." %&% gencodeset,quote=F,col.names=F,row.names=F)
        runGCTAgrm <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "mldose.gz " %&% gt.dir %&% machpre %&% "mlinfo.gz --make-grm-bin --extract tmp." %&% thresh %&% ".SNPlist." %&% gencodeset %&% " --out " %&% grm.dir %&% "trans-" %&% ensid %&% " --thread-num 8"
        system(runGCTAgrm)
    }
}
