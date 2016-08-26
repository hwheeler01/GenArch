args <- commandArgs(trailingOnly=T)
"%&%" = function(a,b) paste(a,b,sep="")
date = Sys.Date()
library(data.table)
library(dplyr)
my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"

tislist <- c("DGN-WB_h2_joint_global","DGN-WB_h2_joint_local")
sevec <- vector()

for(tis in tislist){
  print(tis)
  permdata <- data.frame(fread(my.dir %&% tis %&% "_reml-no-constrain_allgenes_100perms_2016-08-12.txt"))
  print(summary(permdata$V5))
  perms <- as.matrix(permdata[,-1:-5])
  genemeanh2 <- data.frame(meanh2=rowMeans(perms,na.rm=TRUE))
  print(dim(genemeanh2))
  print(summary(genemeanh2))
  genemeanh2 <- dplyr::filter(genemeanh2,meanh2<1 & meanh2>-1) #rm outliers
  print(dim(genemeanh2))
  se <- signif(sd(genemeanh2$meanh2,na.rm=TRUE),2)
  sevec <- c(sevec, se)
}

res<-cbind(tislist,sevec)
colnames(res) <- c('tissue','se')
write.table(res,file=my.dir %&% "SE_estimate_from_h2_joint_reml-no-constrain_allgenes_100perms.txt", quote=F,row.names = F,sep='\t')
