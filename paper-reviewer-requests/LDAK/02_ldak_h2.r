####by Heather E. Wheeler 20141114####
args <- commandArgs(trailingOnly=T)
args <- "22"
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(dplyr)
library(data.table)

pre.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
my.dir <- pre.dir %&% "expArch_DGN-WB_imputedGTs/"
rna.dir <- my.dir %&% "dgn-exp/"
annot.dir <- pre.dir %&% "gtex-annot/"
gt.dir <- "/group/im-lab/nas40t2/Data/Transcriptome/WB1K/imputed/DGN-imputed-for-PrediXcan/"

gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein.chr" %&% args[1]
gencodeset <- args[1]

grm.dir <- my.dir %&% "dgn-grms/"

machpre <- "DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU."

bimfile <- gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".bim" ###get SNP position information###
bim <- read.table(bimfile)
rownames(bim) <- bim$V2

###get dosage data
#gt <- readRDS(gt.dir %&% "DGN.imputed_maf0.05_R20.8_1000G.chr" %&% args[1] %&% ".SNPxID.rds")
gt <- read.table(gt.dir %&% "DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr" %&% args[1] %&% ".SNPxID")

### Scan expression data
idfile <- rna.dir %&% "DGN-WBexp.ID.list"
genefile <- rna.dir %&% "DGN-WBexp.GENE.list"
expfile <- rna.dir %&% "DGN-WB.rntransform.exp.IDxGENE"

subjid <- scan(idfile,"character")
geneid <- scan(genefile, "character")
expdata <- scan(expfile)
expdata <- matrix(expdata, ncol = length(geneid), byrow=TRUE)
rownames(expdata) <- subjid
colnames(expdata) <- geneid

gencode <- read.table(gencodefile)
colnames(gencode) <- c('chr','str','start','end','ensid','gene','func','known')

### get chr gene exp data
genebool <- colnames(expdata) %in% gencode$gene
chrexp <- expdata[,genebool==TRUE]
expgencode <- gencode[gencode$gene %in% colnames(expdata),]
ensidlist <- expgencode[,5]

for(i in 1:length(ensidlist)){
    cat(i,"/",length(ensidlist),"\n")
    ensid <- ensidlist[i]
    geneinfo <- gencode[as.character(ensid),]
    gene <- geneinfo[1,6]
    chr <- geneinfo[1,1]
    c <- substr(chr,4,5)
    start <- geneinfo$start - 1e6 ### 1Mb lower bound for cis-eQTLS
    end <- geneinfo$end + 1e6 ### 1Mb upper bound for cis-eQTLs
    genefile <- data.frame(as.character(ensid),c,start,end)
    write.table(genefile,file="genefile",quote=F,row.names=F,col.names=F)

    exppheno <- data.frame(chrexp[,as.character(gene)])
    colnames(exppheno) <- c('exp')
    phenofile <- mutate(exppheno, FID=rownames(exppheno), IID=1) %>% dplyr::select(FID, IID, exp)
    write.table(phenofile,file="phenofile",quote=F,row.names=F)

    chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
    cissnps <- subset(chrsnps,chrsnps[,4]>=start & chrsnps[,4]<=end) ### pull cis-SNP info
    snplist <- cissnps[,2]    
    write.table(snplist, file= my.dir %&% "tmp." %&% thresh %&% ".SNPlist." %&% gencodeset,quote=F,col.names=F,row.names=F)
    runGCTAgrm <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&%  ".mldose.gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mlinfo.gz --make-grm-bin --extract tmp." %&% thresh %&% ".SNPlist." %&% gencodeset %&% " --out " %&% grm.dir %&% "local-" %&% gene
    system(runGCTAgrm)
    nonlocal <- setdiff(snplist,eqtllist)
    if(length(nonlocal)==0){ ##condition when no local SNPs are eQTLs in FHS, generate global grm using all FHS eQTLs on Chr to ensure global-GENE_Chr.grm* files will be made
        runGCTAgrm <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&%  ".mldose.gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mlinfo.gz --make-grm-bin --extract " %&% eqtlfile %&% " --out " %&% grm.dir %&% "global-" %&% gene %&% "-Chr" %&% gencodeset
    }else{
        write.table(nonlocal, file= my.dir %&% "tmp." %&% thresh %&% ".SNPlist." %&% gencodeset,quote=F,col.names=F,row.names=F)
    	runGCTAgrm <- "gcta64 --dosage-mach-gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&%  ".mldose.gz " %&% gt.dir %&% machpre %&% "chr" %&% gencodeset %&% ".mlinfo.gz --make-grm-bin --extract tmp." %&% thresh %&% ".SNPlist." %&% gencodeset %&% " --out " %&% grm.dir %&% "global-" %&% gene %&% "-Chr" %&% gencodeset
    }
    system(runGCTAgrm)
}


./ldak.4.9 --cut-genes testgene --sp dgn-sp-format/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr22 --genefile genefile --weights chr22/weightsALL 
./ldak.4.9 --calc-genes-reml testgene --pheno phenofile --sp dgn-sp-format/DGN.imputed_maf0.05_R20.8.hapmapSnpsCEU.chr22 --partition 1 --weights chr22/weightsALL
