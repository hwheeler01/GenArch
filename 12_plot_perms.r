date <- Sys.Date()
library(dplyr)
library(ggplot2)
library(gridExtra)
se <- function(x) sqrt(var(x,na.rm=T)/sum(!is.na(x)))
med <- function(x) median(x,na.rm=T)
lci<-function(x) quantile(x, c(.025, 0.975),na.rm=T)[1]
uci<-function(x) quantile(x, c(.025, 0.975),na.rm=T)[2]
"%&%" = function(a,b) paste(a,b,sep="")
mydir <- "/Users/heather/Dropbox/cross-tissue/"
tissues <- scan('tissue.list.prefix','character',sep='\t')


###plot all genes
for(i in 1:length(tissues)){
  tis <- tissues[i]

  ##local only
  loconly <- read.table(mydir %&% tis %&% '_h2_localonly_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- loconly[,5]
  perms <- loconly[,6:105]
  meanperms <- rowMeans(perms,na.rm=T)
  seperms <- apply(perms,1,se)
    
  data <- data.frame(obs,meanperms,seperms)
  data <- data %>% mutate(ymin = pmax(0, meanperms - 2 * seperms), ymax = pmin(1, meanperms + 2 * seperms) )
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]

  p1 <- ggplot(data,aes(x=1:nrow(data),y=meanperms,ymin=ymin, ymax=ymax) ) + geom_pointrange(col='gray')+geom_point(col='gray')+ggtitle(tis %&% ' Local (Only) Heritability')
  loplot <- p1 + geom_point(data=data,aes(x=1:nrow(data),y=obs),color='red') + xlab("Genes sorted by observed h2") + ylab("h2") + ylim(0,1)
  print(loplot)
  
  ##global only
  gloonly <- read.table(tis %&% '_h2_globalonly_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- gloonly[,5]
  perms <- gloonly[,6:105]
  meanperms <- rowMeans(perms,na.rm=T)
  seperms <- apply(perms,1,se)
  
  data <- data.frame(obs,meanperms,seperms)
  data <- data %>% mutate(ymin = pmax(0, meanperms - 2 * seperms), ymax = pmin(1, meanperms + 2 * seperms) )
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  
  p1 <- ggplot(data,aes(x=1:nrow(data),y=meanperms,ymin=ymin, ymax=ymax) ) + geom_pointrange(col='gray')+geom_point(col='gray')+ggtitle(tis %&% ' Global (Only) Heritability')
  goplot <- p1 + geom_point(data=data,aes(x=1:nrow(data),y=obs),color='red') + xlab("Genes sorted by observed h2") + ylab("h2") + ylim(0,1)
  print(goplot)

  ##local joint
  locjoint <- read.table(tis %&% '_h2_localjoint_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- locjoint[,5]
  perms <- locjoint[,6:105]
  meanperms <- rowMeans(perms,na.rm=T)
  seperms <- apply(perms,1,se)
  
  data <- data.frame(obs,meanperms,seperms)
  data <- data %>% mutate(ymin = pmax(0, meanperms - 2 * seperms), ymax = pmin(1, meanperms + 2 * seperms) )
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  
  p1 <- ggplot(data,aes(x=1:nrow(data),y=meanperms,ymin=ymin, ymax=ymax) ) + geom_pointrange(col='gray')+geom_point(col='gray')+ggtitle(tis %&% ' Local (Joint) Heritability')
  ljplot <- p1 + geom_point(data=data,aes(x=1:nrow(data),y=obs),color='red') + xlab("Genes sorted by observed h2") + ylab("h2") + ylim(0,1)
  print(ljplot)
  
  ##global joint
  glojoint <- read.table(tis %&% '_h2_globaljoint_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- glojoint[,5]
  perms <- glojoint[,6:105]
  meanperms <- rowMeans(perms,na.rm=T)
  seperms <- apply(perms,1,se)
    
  data <- data.frame(obs,meanperms,seperms)
  data <- data %>% mutate(ymin = pmax(0, meanperms - 2 * seperms), ymax = pmin(1, meanperms + 2 * seperms) )
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  
  p1 <- ggplot(data,aes(x=1:nrow(data),y=meanperms,ymin=ymin, ymax=ymax) ) + geom_pointrange(col='gray')+geom_point(col='gray')+ggtitle(tis %&% ' Global (Joint) Heritability')
  gjplot <- p1 + geom_point(data=data,aes(x=1:nrow(data),y=obs),color='red') + xlab("Genes sorted by observed h2") + ylab("h2") + ylim(0,1)
  print(gjplot)
}

###plot genes w/h2 > 0.05
for(i in 1:length(tissues)){
  tis <- tissues[i]
  
  ##local only
  loconly <- read.table(mydir %&% tis %&% '_h2_localonly_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- loconly[,5]
  perms <- loconly[,6:105]
  meanperms <- rowMeans(perms,na.rm=T)
  seperms <- apply(perms,1,se)
  
  data <- data.frame(obs,meanperms,seperms)
  data <- data %>% mutate(ymin = pmax(0, meanperms - 2 * seperms), ymax = pmin(1, meanperms + 2 * seperms) )
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  gt05 <- data[data[,1]>0.05,] ##only genes with obs h2 > 0.05

  p1 <- ggplot(gt05,aes(x=1:nrow(gt05),y=meanperms,ymin=ymin, ymax=ymax) ) + geom_pointrange(col='gray')+geom_point(col='gray')+ggtitle(tis %&% ' Local (Only) Heritability')
  losub <- p1 + geom_point(data=gt05,aes(x=1:nrow(gt05),y=obs),color='red') + xlab("Genes sorted by observed h2 > 0.05") + ylab("h2") + ylim(0,1)
  print(losub)
  
  ##global only
  gloonly <- read.table(tis %&% '_h2_globalonly_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- gloonly[,5]
  perms <- gloonly[,6:105]
  meanperms <- rowMeans(perms,na.rm=T)
  seperms <- apply(perms,1,se)
  
  data <- data.frame(obs,meanperms,seperms)
  data <- data %>% mutate(ymin = pmax(0, meanperms - 2 * seperms), ymax = pmin(1, meanperms + 2 * seperms) )
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  gt05 <- data[data[,1]>0.05,] ##only genes with obs h2 > 0.05
  
  p1 <- ggplot(gt05,aes(x=1:nrow(gt05),y=meanperms,ymin=ymin, ymax=ymax) ) + geom_pointrange(col='gray')+geom_point(col='gray')+ggtitle(tis %&% ' Global (Only) Heritability')
  gosub <- p1 + geom_point(data=gt05,aes(x=1:nrow(gt05),y=obs),color='red') + xlab("Genes sorted by observed h2 > 0.05") + ylab("h2") + ylim(0,1)
  print(gosub)
  
  ##local joint
  locjoint <- read.table(tis %&% '_h2_localjoint_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- locjoint[,5]
  perms <- locjoint[,6:105]
  meanperms <- rowMeans(perms,na.rm=T)
  seperms <- apply(perms,1,se)
  
  data <- data.frame(obs,meanperms,seperms)
  data <- data %>% mutate(ymin = pmax(0, meanperms - 2 * seperms), ymax = pmin(1, meanperms + 2 * seperms) )
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  gt05 <- data[data[,1]>0.05,] ##only genes with obs h2 > 0.05
  
  p1 <- ggplot(gt05,aes(x=1:nrow(gt05),y=meanperms,ymin=ymin, ymax=ymax) ) + geom_pointrange(col='gray')+geom_point(col='gray')+ggtitle(tis %&% ' Local (Joint) Heritability')
  ljsub <- p1 + geom_point(data=gt05,aes(x=1:nrow(gt05),y=obs),color='red') + xlab("Genes sorted by observed h2 > 0.05") + ylab("h2") + ylim(0,1)
  print(ljsub)
  
  ##global joint
  glojoint <- read.table(tis %&% '_h2_globaljoint_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- glojoint[,5]
  perms <- glojoint[,6:105]
  meanperms <- rowMeans(perms,na.rm=T)
  seperms <- apply(perms,1,se)
  
  data <- data.frame(obs,meanperms,seperms)
  data <- data %>% mutate(ymin = pmax(0, meanperms - 2 * seperms), ymax = pmin(1, meanperms + 2 * seperms) )
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  gt05 <- data[data[,1]>0.05,] ##only genes with obs h2 > 0.05
  
  p1 <- ggplot(gt05,aes(x=1:nrow(gt05),y=meanperms,ymin=ymin, ymax=ymax) ) + geom_pointrange(col='gray')+geom_point(col='gray')+ggtitle(tis %&% ' Global (Joint) Heritability')
  gjsub <- p1 + geom_point(data=gt05,aes(x=1:nrow(gt05),y=obs),color='red') + xlab("Genes sorted by observed h2 > 0.05") + ylab("h2") + ylim(0,1)
  print(gjsub)

}

#############################################################################################################

###plot genes w/h2 > 0.05, plot median permuted h2 and 95% CI of median
for(i in 1:length(tissues)){
  tis <- tissues[i]
  
  ##local only
  loconly <- read.table(mydir %&% tis %&% '_h2_localonly_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- loconly[,5]
  perms <- loconly[,6:105]
  medperms <- apply(perms,1,med)
  seperms <- apply(perms,1,se)
  minperms <- apply(perms,1,lci)
  maxperms <- apply(perms,1,uci)
  
  data <- data.frame(obs,medperms,seperms,minperms,maxperms)
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  gt05 <- data[data[,1]>0.05,] ##only genes with obs h2 > 0.05
  
  p1 <- ggplot(gt05,aes(x=1:nrow(gt05),y=medperms,ymin=minperms, ymax=maxperms) ) + geom_pointrange(col='gray')+geom_point(col='black')+ggtitle(tis %&% ' Local (Only) Heritability')
  losub <- p1 + geom_point(data=gt05,aes(x=1:nrow(gt05),y=obs),color='red') + xlab("Genes sorted by observed h2 > 0.05") + ylab("h2") + ylim(0,1)
  print(losub)
  
  ##global only
  gloonly <- read.table(tis %&% '_h2_globalonly_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- gloonly[,5]
  perms <- gloonly[,6:105]
  medperms <- apply(perms,1,med)
  seperms <- apply(perms,1,se)
  minperms <- apply(perms,1,lci)
  maxperms <- apply(perms,1,uci)
  
  data <- data.frame(obs,medperms,seperms,minperms,maxperms)
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  gt05 <- data[data[,1]>0.05,] ##only genes with obs h2 > 0.05
  
  p1 <- ggplot(gt05,aes(x=1:nrow(gt05),y=medperms,ymin=minperms, ymax=maxperms) ) + geom_pointrange(col='gray')+geom_point(col='black')+ggtitle(tis %&% ' Global (Only) Heritability')
  gosub <- p1 + geom_point(data=gt05,aes(x=1:nrow(gt05),y=obs),color='red') + xlab("Genes sorted by observed h2 > 0.05") + ylab("h2") + ylim(0,1)
  print(gosub)
  
  ##local joint
  locjoint <- read.table(tis %&% '_h2_localjoint_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- locjoint[,5]
  perms <- locjoint[,6:105]
  medperms <- apply(perms,1,med)
  seperms <- apply(perms,1,se)
  minperms <- apply(perms,1,lci)
  maxperms <- apply(perms,1,uci)
  
  data <- data.frame(obs,medperms,seperms,minperms,maxperms)
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  gt05 <- data[data[,1]>0.05,] ##only genes with obs h2 > 0.05
  
  p1 <- ggplot(gt05,aes(x=1:nrow(gt05),y=medperms,ymin=minperms, ymax=maxperms) ) + geom_pointrange(col='gray')+geom_point(col='black')+ggtitle(tis %&% ' Local (Joint) Heritability')
  ljsub <- p1 + geom_point(data=gt05,aes(x=1:nrow(gt05),y=obs),color='red') + xlab("Genes sorted by observed h2 > 0.05") + ylab("h2") + ylim(0,1)
  print(ljsub)
  
  ##global joint
  glojoint <- read.table(tis %&% '_h2_globaljoint_allgenes_100perms_2014-12-19.txt',sep='\t')
  obs <- glojoint[,5]
  perms <- glojoint[,6:105]
  medperms <- apply(perms,1,med)
  seperms <- apply(perms,1,se)
  minperms <- apply(perms,1,lci)
  maxperms <- apply(perms,1,uci)
  
  data <- data.frame(obs,medperms,seperms,minperms,maxperms)
  data <- data %>% arrange(obs) 
  data <- data[complete.cases(data),]
  gt05 <- data[data[,1]>0.05,] ##only genes with obs h2 > 0.05
  
  p1 <- ggplot(gt05,aes(x=1:nrow(gt05),y=medperms,ymin=minperms, ymax=maxperms) ) + geom_pointrange(col='gray')+geom_point(col='black')+ggtitle(tis %&% ' Global (Joint) Heritability')
  gjsub <- p1 + geom_point(data=gt05,aes(x=1:nrow(gt05),y=obs),color='red') + xlab("Genes sorted by observed h2 > 0.05") + ylab("h2") + ylim(0,1)
  print(gjsub)
  
}

#pdf(file="100perms_h2plot_sorted_" %&% tis %&% ".pdf")
#grid.arrange(ljplot,gjplot)
#dev.off()
