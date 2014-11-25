"%&%" = function(a,b) paste(a,b,sep="")
library(scatterplot3d)
date <- Sys.Date()
grms<-read.table('GTEx.cross-tissue.h2.all.models_2014-11-25.txt',header=T,sep="\t")

pdf(file="GTEx_cross-tissue_global_v_local_h2_" %&% date %&% ".pdf",height=5,width=5)
plot(grms$local.h2,grms$global.h2,xlim=c(0,1),ylim=c(0,1),cex=0.5,xlab="localGRM h2",ylab="globalGRM h2",
     main="GTEx Cross-tissue expression GCTA",col=rgb(100,100,100,100,maxColorValue=255))
abline(1,-1,col="red")
dev.off()

grm3<-cbind(grms$local.h2.1,grms$chr.h2,grms$global.h2.1)
pdf(file="GTEx_cross-tissue_global_v_local_v_chr_h2_" %&% date %&% ".pdf",height=5,width=5)
scatterplot3d(grm3[,1],grm3[,2],grm3[,3],pch=16, highlight.3d=TRUE,cex.symbols=0.5,
              type="p",xlab="localGRM h2",ylab="chrGRM h2",zlab="globalGRM h2",main="GTEx Cross-tissue expression GCTA")
dev.off()

tissues<-read.table("tissue.list",sep="\n")

for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  
  pdf(file="GTEx_resid_" %&% tis %&% "_global_v_local_h2_" %&% date %&% ".pdf",height=5,width=5)
  plot(tisgrms$local.h2,tisgrms$global.h2,xlim=c(0,1),ylim=c(0,1),cex=0.5,xlab="localGRM h2",ylab="globalGRM h2",
       main="GTEx " %&% tis %&% " expression GCTA",col=rgb(100,100,100,100,maxColorValue=255))
  abline(1,-1,col="red")
  dev.off()
  
  tisgrm3<-cbind(tisgrms$local.h2.1,tisgrms$chr.h2,tisgrms$global.h2.1)
  pdf(file="GTEx_resid_" %&% tis %&% "_global_v_local_v_chr_h2_" %&% date %&% ".pdf",height=5,width=5)
  scatterplot3d(tisgrm3[,1],tisgrm3[,2],tisgrm3[,3],pch=16, highlight.3d=TRUE,cex.symbols=0.5,
                type="p",xlab="localGRM h2",ylab="chrGRM h2",zlab="globalGRM h2",main="GTEx " %&% tis %&% " expression GCTA")
  dev.off()
}

## Y ~ localGRM
###plot histograms of h2

png(file="hist_GTEx_localGRM_h2_" %&% date %&% ".png",height=6,width=10,units="in",res=100)
par(mfcol=c(2,4))
###no ylim
hist(grms$h2,ylim=c(0,15000),xlim=c(0,1),main="Cross-tissue",xlab="h2")
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$h2,xlim=c(0,1),main=tis,xlab="h2")
}
dev.off()
###ylim = 3000
png(file="hist_GTEx_localGRM_h2_ylim3000_" %&% date %&% ".png",height=6,width=10,units="in",res=100)
par(mfcol=c(2,4))
hist(grms$h2,ylim=c(0,3000),xlim=c(0,0.5),main="Cross-tissue",xlab="h2")
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$h2,ylim=c(0,3000),xlim=c(0,0.5),main=tis,xlab="h2")
}
dev.off()

###plot histograms of p
png(file="hist_GTEx_localGRM_p_" %&% date %&% ".png",height=6,width=10,units="in",res=100)
par(mfcol=c(2,4))
###no ylim
hist(grms$p,main="Cross-tissue",xlab="p-value",ylim=c(0,10000))
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$p,main=tis,xlab="p-value",ylim=c(0,10000))
}
dev.off()

## Y ~ localGRM + globalGRM
###plot histograms of h2

png(file="hist_GTEx_localGRM_globalGRM_h2_" %&% date %&% ".png",height=10,width=10,units="in",res=100)
par(mfcol=c(4,4))
###ylim = 2000
hist(grms$local.h2,ylim=c(0,2000),xlim=c(0,1),main="Cross-tissue",xlab="localGRM h2")
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$local.h2,ylim=c(0,2000),xlim=c(0,1),main=tis,xlab="localGRM h2")
}
hist(grms$global.h2,ylim=c(0,2000),xlim=c(0,1),main="Cross-tissue",xlab="globalGRM h2")
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$global.h2,ylim=c(0,2000),xlim=c(0,1),main=tis,xlab="globalGRM h2")
}
dev.off()

###plot histograms of se
png(file="hist_GTEx_localGRM_globalGRM_se_" %&% date %&% ".png",height=10,width=10,units="in",res=100)
par(mfcol=c(4,4))
###ylim = 5000
hist(grms$local.se,xlim=c(0,1),main="Cross-tissue",xlab="localGRM SE",ylim=c(0,5000))
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$local.se,xlim=c(0,1),main=tis,xlab="localGRM SE",ylim=c(0,5000))
}
hist(grms$global.se,xlim=c(0,1),main="Cross-tissue",xlab="globalGRM SE",ylim=c(0,5000))
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$global.se,xlim=c(0,1),main=tis,xlab="globalGRM SE",ylim=c(0,5000))
}
dev.off()

## Y ~ localGRM + chrGRM + globalGRM
###plot histograms of h2

png(file="hist_GTEx_localGRM_chrGRM_globalGRM_h2_" %&% date %&% ".png",height=10,width=15,units="in",res=100)
par(mfcol=c(4,6))
###ylim = 2000
hist(grms$local.h2.1,ylim=c(0,2000),xlim=c(0,1),main="Cross-tissue",xlab="localGRM h2")
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$local.h2.1,ylim=c(0,2000),xlim=c(0,1),main=tis,xlab="localGRM h2")
}
hist(grms$chr.h2,ylim=c(0,2000),xlim=c(0,1),main="Cross-tissue",xlab="chrGRM h2")
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$chr.h2,ylim=c(0,2000),xlim=c(0,1),main=tis,xlab="chrGRM h2")
}
hist(grms$global.h2.1,ylim=c(0,2000),xlim=c(0,1),main="Cross-tissue",xlab="globalGRM h2")
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$global.h2.1,ylim=c(0,2000),xlim=c(0,1),main=tis,xlab="globalGRM h2")
}
dev.off()

###plot histograms of se
png(file="hist_GTEx_localGRM_chrGRM_globalGRM_se_" %&% date %&% ".png",height=10,width=15,units="in",res=100)
par(mfcol=c(4,6))
###ylim = 4000
hist(grms$local.se.1,xlim=c(0,1),main="Cross-tissue",xlab="localGRM SE",ylim=c(0,4000))
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$local.se.1,xlim=c(0,1),main=tis,xlab="localGRM SE",ylim=c(0,4000))
}
hist(grms$chr.se,xlim=c(0,1),main="Cross-tissue",xlab="chrGRM SE",ylim=c(0,4000))
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$chr.se,xlim=c(0,1),main=tis,xlab="chrGRM SE",ylim=c(0,4000))
}
hist(grms$global.se.1,xlim=c(0,1),main="Cross-tissue",xlab="globalGRM SE",ylim=c(0,4000))
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  hist(tisgrms$global.se.1,xlim=c(0,1),main=tis,xlab="globalGRM SE",ylim=c(0,4000))
}
dev.off()


###plot scatterplots global v. local h2 on one png
png(file="scatter_GTEx_localGRM_globalGRM_h2_" %&% date %&% ".png",height=6,width=10,units="in",res=100)
par(mfcol=c(2,4))
###no ylim
plot(grms$local.h2,grms$global.h2,xlim=c(0,1),ylim=c(0,1),cex=0.5,xlab="localGRM h2",ylab="globalGRM h2",
     main="Cross-tissue",col=rgb(100,100,100,100,maxColorValue=255))
abline(1,-1,col="red")
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  plot(tisgrms$local.h2,tisgrms$global.h2,xlim=c(0,1),ylim=c(0,1),cex=0.5,xlab="localGRM h2",ylab="globalGRM h2",
       main=tis,col=rgb(100,100,100,100,maxColorValue=255))
  abline(1,-1,col="red")
}
dev.off()

###plot 3D scatterplots global v. local v. chr h2 on one png
png(file="scatter_GTEx_localGRM_chrGRM_globalGRM_h2_" %&% date %&% ".png",height=6,width=10,units="in",res=100)
par(mfcol=c(2,4))
###no ylim
scatterplot3d(grm3[,1],grm3[,2],grm3[,3],pch=16, highlight.3d=TRUE,cex.symbols=0.5,cex.lab=0.7,cex.axis=0.5,xlim=c(0,1),ylim=c(0,1),zlim=c(0,1),
              type="p",xlab="localGRM h2",ylab="chrGRM h2",zlab="globalGRM h2",main="Cross-tissue")
for(tis in tissues[,1]){
  file<- "GTEx.resid.tissue-specific.h2_" %&% tis %&% "_all.models_2014-11-24.txt"
  tisgrms<-read.table(file,header=T,sep="\t")
  tisgrm3<-cbind(tisgrms$local.h2.1,tisgrms$chr.h2,tisgrms$global.h2.1)
  scatterplot3d(tisgrm3[,1],tisgrm3[,2],tisgrm3[,3],pch=16, highlight.3d=TRUE,cex.symbols=0.5,cex.lab=0.7,cex.axis=0.5,
                xlim=c(0,1),ylim=c(0,1),zlim=c(0,1),type="p",xlab="localGRM h2",ylab="chrGRM h2",zlab="globalGRM h2",main=tis)
}
dev.off()