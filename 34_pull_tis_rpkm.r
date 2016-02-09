####by Heather E. Wheeler 20140602####
"%&%" = function(a,b) paste(a,b,sep="")
args <- commandArgs(trailingOnly=T)
date <- Sys.Date() 
###############################################
### Directories & Variables
 
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
rna.dir <- my.dir %&% "gtex-rnaseq/"
annot.dir <- my.dir %&% "gtex-annot/"
out.dir <- rna.dir %&% "ind-tissues-RPKM/"
 
tissue <- "Nerve - Tibial" ###check GTEx_Analysis_2014-06-13.SampleTissue.annot for available tissues###

tislist <- scan(my.dir %&% "nine.spaces.tissue.list","c",sep="\n")

for(tissue in tislist){ 
	sam <- read.table(annot.dir %&% "GTEx_Analysis_2014-06-13.SampleTissue.annot",header=T,sep="\t") 
	### above file includes SAMPID and SMTSD from: 
	### /group/im-lab/nas40t2/haky/Data/dbGaP/GTEx/41400/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2014-06-13/sample_annotations/GTEx_Data_2014-06-13_Annotations_SampleAttributesDS.txt
	sample <- subset(sam,SMTSD == tissue) ### pull sample list of chosen tissue###
 
	expidlist <- scan(rna.dir %&% "GTEx_Analysis_2014-06-13.RNA-seq.ID.list","character")
	expgenelist <- scan(rna.dir %&% "GTEx_Analysis_2014-06-13.RNA-seq.GENE.list","character")
	exp <- scan(rna.dir %&% "GTEx_Analysis_2014-06-13.RNA-seq.GENExID")
	expdata <- matrix(exp, ncol=length(expidlist), byrow=T)
	t.expdata <- t(expdata)	
	rownames(t.expdata) <- expidlist
	colnames(t.expdata) <- expgenelist

	gencodefile <- annot.dir %&% "gencode.v18.genes.patched_contigs.summary.protein"
	gencode <- read.table(gencodefile) 
	rownames(gencode) <- gencode[,5]
	t.expdata <- t.expdata[,intersect(colnames(t.expdata),rownames(gencode))] ###pull protein coding gene expression data
	
	tissue.exp <- t.expdata[intersect(rownames(t.expdata),sample$SAMPID),] ###pull expression data for chosen tissue###

	tissue.exp <- t(tissue.exp) #for merging in R

	write.table(tissue.exp, file= out.dir %&% tissue %&% ".RPKM.GENExID", quote=F, row.names=F, col.names=F)
	write(colnames(tissue.exp), file = out.dir %&% tissue %&% ".RPKM.GENE.list", ncolumns=1)
	write(rownames(tissue.exp), file = out.dir %&% tissue %&% ".RPKM.ID.list", ncolumns=1)

	##calc mean for each gene
	meanRPKM <- apply(tissue.exp,1,mean)
	meanRPKM <- data.frame(ensid=names(meanRPKM),meanRPKM)
	write.table(meanRPKM, file= out.dir %&% tissue %&% ".meanRPKM.txt", quote=F, row.names=F)
}
