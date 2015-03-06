####by Heather E. Wheeler 20141114####
args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()

###make DGN GRMs: local (unique to gene: localGENE.grm), global (all eQTLs from Framingham (FHS): globalChr.grm), global-local (unique to gene:global-GENE_Chr.grm)

pre.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
my.dir <- pre.dir %&% "expArch_DGN-WB_imputedGTs/"
rna.dir <- my.dir %&% "dgn-exp/"
annot.dir <- pre.dir %&% "gtex-annot/"
gt.dir <- "/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-imputation/DGN-imputed-for-PrediXcan/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% args[1]
gencodeset <- args[1]

thresh <- 'fdr0.05'
eqtlfile <- my.dir %&% "Framingham_eqtl-gene_" %&% thresh %&% "_hapmapSnpsCEU.rsIDlist"
eqtllist <- scan(eqtlfile,"character")

grm.dir <- my.dir %&% "dgn-grms/FHS_eQTLs_" %&% thresh %&% "/"

machpre <- "DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU."

bimfile <- gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".bim" ###get SNP position information###
bim <- read.table(bimfile)
rownames(bim) <- bim$V2

###make globalGRM for chr, only need to do once

runGCTAglo <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mldose.gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mlinfo.gz --extract " %&% eqtlfile %&% " --make-grm-bin --out " %&% grm.dir %&% "DGN.global_Chr" %&% gencodeset
system(runGCTAglo)

###make localGRMs, for subset of genes in gencodeset (chr) and global minus local SNPs for the chr, these are both gene specific

gencode <- read.table(gencodefile)
rownames(gencode) <- gencode[,5]

ensidlist <- gencode[,5]

for(i in 1:length(ensidlist)){
    cat(i,"/",length(ensidlist),"\n")
    ensid <- ensidlist[i]
    geneinfo <- gencode[ensid,]
    gene <- geneinfo[1,6]
    chr <- geneinfo[1,1]
    c <- substr(chr,4,5)
    start <- geneinfo$V3 - 1e6 ### 1Mb lower bound for cis-eQTLS
    end <- geneinfo$V4 + 1e6 ### 1Mb upper bound for cis-eQTLs
    chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
    cissnps <- subset(chrsnps,chrsnps[,4]>=start & chrsnps[,4]<=end) ### pull cis-SNP info
    snplist <- cissnps[,2]    
    write.table(snplist, file= my.dir %&% "tmp." %&% thresh %&% ".SNPlist." %&% gencodeset,quote=F,col.names=F,row.names=F)
    runGCTAgrm <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&%  ".mldose.gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mlinfo.gz --make-grm-bin --extract tmp." %&% thresh %&% ".SNPlist." %&% gencodeset %&% " --out " %&% grm.dir %&% "local-" %&% gene
    system(runGCTAgrm)
    nonlocal <- setdiff(snplist,eqtllist)
    write.table(nonlocal, file= my.dir %&% "tmp." %&% thresh %&% ".SNPlist." %&% gencodeset,quote=F,col.names=F,row.names=F)
    runGCTAgrm <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&%  ".mldose.gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mlinfo.gz --make-grm-bin --extract tmp." %&% thresh %&% ".SNPlist." %&% gencodeset %&% " --out " %&% grm.dir %&% "global-" %&% gene %&% "-Chr" %&% gencodeset
    system(runGCTAgrm)
}
