From: David DeLuca [mailto:ddeluca@broadinstitute.org] 
Sent: Monday, November 10, 2014 10:35 AM
To: Ayellet Segre
Cc: Gamazon, Eric [BSD] - MED
Subject: Re: quick question
 
Most relevant normalization steps in bold
 
# Reads in gct expression values and produces the  matrix used to load into eQTL analysis
  
args <- commandArgs(trailingOnly = TRUE)

# $1 = gct
# $2 = output file
# $3 = variation threshold (try 0.3)
# $4 = genotype files
# $5 = expression theshold (use -5) or 1 for super strict
# $6 = percent of expressing samples (use 0.1)
# $7 = trim headers to idividual level (true or false)
# $8 = quantile normalize data (true or false)
# $9 = include description field? true or false
# $10 = standard normalize within genes
# $11 = median preservation normalization within genes
 
 
gct<-read.delim(args[1],as.is=TRUE,header=T,skip=2,check.names=FALSE)       
outFile <- args[2]
variationThresh <- as.numeric(args[3])
genos<-read.delim(args[4],as.is=TRUE,header=F,check.names=FALSE)     
exprThresh <- as.numeric(args[5])
a<-args[6]
percentExpressed<-NULL
if(length(grep("%",a))>0){
        percentExpressed <- as.numeric(strtrim(a,nchar(a)-1))/100
}
cat("Parameters: " ,args[1], outFile, variationThresh,"\n")
 
getIndivid<-function(s){
        spl<-strsplit(s,"-")
        pasteIt<-function(toks) {paste(toks[1],toks[2],sep="-")      }
        unlist(lapply(spl,pasteIt))
}
 
 
# remove the descirption column:
#gct[[2]]<-NULL
names(gct)[1]<-"Id"
 
idFieldCount <- 1
if(args[9]){
      idFieldCount<-2
}
 
 
DATA_COLS<-(idFieldCount+1):ncol(gct)
gct<- gct[,c(rep(TRUE,idFieldCount),getIndivid(names(gct[DATA_COLS])) %in% genos$V1)]
cat("Number of samples with genos: ",ncol(gct)-idFieldCount,"\n")
 
DATA_COLS<-(idFieldCount+1):ncol(gct)
m<-data.matrix(gct[,DATA_COLS])
 
if(is.null(percentExpressed)){
        numExpressed<-as.integer(a)
}else{
        numExpressed<-as.integer(percentExpressed*length(DATA_COLS))
}
 
cat("Number of expressed genes required: ",numExpressed,"\n")
 
# Create a filter based on standard deviation
devs<-apply(m,1,sd)
f1<- (devs>variationThresh)
cat(sum(f1)," pass sd filter\n")
 
# Create a filter based on the number of expressed transcripts
 
exprFilter<-function(x, exprThresh, numExpressedThresh){
        x<-x[!is.na(x)]
        l<-length(x)
sum(x>exprThresh) > numExpressedThresh
}
 
f2<-apply(m,1,exprFilter,exprThresh=exprThresh,numExpressedThresh=numExpressed)
cat(sum(f2)," pass min expr filter\n")
 
f<- f1 & f2
cat("Writing output file with ",sum(f),"transcripts\n")
 
# filter rows and order columns alphabetically (so they  match genotypes):
gct<-gct[f,c(seq(1,idFieldCount),idFieldCount+order(names(gct)[DATA_COLS]))]
 
 
# trim headers)
if(args[7]){
        names(gct)[DATA_COLS]<-getIndivid(names(gct[DATA_COLS]))
}
 
# log and quantile normalize
if(args[8]){
        m<-data.matrix(gct[,DATA_COLS])
        shift<-0.01
        m<-log2(m+shift)
        suppressMessages(library(preprocessCore))
        mNorm<-normalize.quantiles(m)
        gct[,DATA_COLS]<-m
}
 
# quantile normalize within genes
if (args[10]){
        m<-data.matrix(gct[,DATA_COLS])
        m = t(apply( m, 1, rank, ties.method = "average"));
        m = qnorm(m / (ncol(m)+1));
        gct[,DATA_COLS]<-m
 
}
 
 
#Andrey's normalization procedure that preserves the median and the median absolute deviation
the.norm = function(x) {
               rks = rank(x, ties.method = "average");
               qnm = qnorm(rks / (length(rks)+1))
               med = median(x);
               mad = median(abs(x - med))
               return(qnm*mad/0.674489750196082 + med); # qnorm(3/4)==0.674489750196082
        }
 
# quantile normalize within genes
if (args[11]){
        m<-data.matrix(gct[,DATA_COLS])
 
        m = t(apply( m, 1, the.norm));
 
        gct[,DATA_COLS]<-m
 
}
 
# output table
 
if(args[9]){
        outFile<-sub("txt","gct",outFile)
        cat("#1.2\n",file=outFile)
        cat(c(nrow(gct),(ncol(gct)-2)),"\n",file=outFile,sep="\t",append=TRUE)
}
 
suppressWarnings(write.table(gct,outFile,append=args[9],quote=F,row.names=F,sep="\t"))
 
cat("Done.\n")

 
	

