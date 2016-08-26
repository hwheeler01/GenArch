args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(data.table)
library(dplyr)
my.dir <- "/Volumes/im-lab/nas40t2/hwheeler/cross-tissue/"
tislist <- scan(my.dir %&% "tissue.list.prefix","c",sep='\t')

sevec <- vector()
for(tis in tislist){
  print(tis)
  permdata <- data.frame(fread(my.dir %&% tis %&% "_h2_localonly_reml-no-constrain_allgenes_100perms_2016-08-03.txt"))
  perms <- as.matrix(permdata[,-1:-5])
  genemeanh2 <- data.frame(meanh2=rowMeans(perms,na.rm=TRUE))
  print(dim(genemeanh2))
  print(summary(genemeanh2))
  genemeanh2 <- dplyr::filter(genemeanh2,meanh2<10 & meanh2>-10) #rm outliers
  print(dim(genemeanh2))
  se <- signif(sd(genemeanh2$meanh2,na.rm=TRUE),2)
  sevec <- c(sevec, se)
}

res<-cbind(tislist,sevec)
colnames(res) <- c('tissue','se')
write.table(res,file=my.dir %&% "SE_estimate_from_h2_localonly_reml-no-constrain_allgenes_100perms.txt", quote=F,row.names = F,sep='\t')
