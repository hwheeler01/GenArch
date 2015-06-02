library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")
date<-Sys.Date()

###convert Nick's GTEx 10-fold CV elastic net results from largest 9 tissues for plotting 

my.dir <- '/group/im-lab/nas40t2/hwheeler/cross-tissue/'
#nick.dir <- '/group/im-lab/nas40t2/nknoblauch/newresults/' ##Adipose-Subcutaneous was missing header row, so I copied everything over and added header
nick.dir <- my.dir %&% 'fromNick_GTEx_TW_elasticNet/'

tislist <- scan(my.dir %&% 'GTEx_nine_tissues','c')

outdata <- data.frame()

for(i in 1:length(tislist)){
  tis <- tislist[i]
  tisdata <- read.table(nick.dir %&% tis %&% '.allResults.txt',header=T)
  interest <- select(tisdata,ensid,R2,alpha) %>% mutate(tissue=tis)
  outdata <- rbind(outdata,interest)
}

write.table(outdata, file="GTEx_tissue-wide_elasticNet_for_ggplot2_" %&% date %&% ".txt",quote=F,row.names=F)
