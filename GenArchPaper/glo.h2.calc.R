#calculate mean global h2 from the joint models in DGN
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(GGally)
library(grid)
library(corrplot)
source('/Volumes/im-lab/nas40t2/hwheeler/cross-tissue/GenArchPaper/multiplot.R')
"%&%" = function(a,b) paste(a,b,sep="")
my.dir <- '/Volumes/im-lab/nas40t2/hwheeler/cross-tissue/'
fig.dir <- '~/GitHub/GenArch/GenArchPaper/Figures/'
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(is.na(x)==FALSE))

## unconstrained model
dgn <- read.table(my.dir %&% 'expArch_DGN-WB_imputedGTs/DGN-WB.h2.all.models_FHSfdr0.05.Chr1-22_globalAll_reml-no-constrain.2015-12-15.txt',header=T)
mean(dgn$loc.jt.h2,na.rm=TRUE)
se(dgn$loc.jt.h2)
mean(dgn$glo.jt.h2,na.rm=TRUE)
mean(dgn$global.se,na.rm=TRUE)
data <- dgn[complete.cases(dgn),]
dim(dgn)
dim(data)

dgn <- read.table(my.dir %&% 'expArch_DGN-WB_imputedGTs/DGN-WB.h2.all.models_FHSfdr0.05.Chr1-22_globaleQTLOtherChr_reml-no-constrain.2015-12-13.txt',header=T)
mean(dgn$loc.jt.h2,na.rm=TRUE)
se(dgn$loc.jt.h2)
mean(dgn$glo.jt.h2,na.rm=TRUE)
mean(dgn$global.se,na.rm=TRUE)
data <- dgn[complete.cases(dgn),]
dim(dgn)
dim(data)

##Fig1 UNCONSTRAINED
#DGN-WB joint heritability. Local h^2^ is estimated with SNPs within 1 Mb of each gene. 
#distal h^2^ is estimated with either all non-chr SNPs or SNPs that are eQTLs in the Framingham Heart Study on other chromosomes (FDR < 0.05).
otherfile<-my.dir %&% 'expArch_DGN-WB_imputedGTs/DGN-WB.h2.all.models_FHSfdr0.05.Chr1-22_globaleQTLOtherChr_reml-no-constrain.2015-12-13.txt'

fdrother<-read.table(otherfile,header=T) ##FHS eQTLs w/fdr<0.05 on non-gene chromosomes used to define global GRM
d  <- fdrother %>% mutate(ymin = glo.jt.h2 - 2 * glo.jt.se, ymax = glo.jt.h2 + 2 * glo.jt.se)
fdrother <- mutate(d, loc.jt.P = ifelse(loc.jt.h2 > 0, 2*pnorm(-abs(loc.jt.h2/loc.jt.se)),1) , glo.jt.P = ifelse(glo.jt.h2 > 0, 2*pnorm(-abs(glo.jt.h2/glo.jt.se)),1)) 
fdrother <- mutate(fdrother, loc.jt.P=ifelse(is.na(loc.jt.P),1,loc.jt.P),glo.jt.P=ifelse(is.na(glo.jt.P),1,glo.jt.P)) %>% 
  mutate(locPlt05=loc.jt.P < 0.05,gloPlt05=glo.jt.P < 0.05) %>% mutate(`distal P`=factor(gloPlt05,labels=c('\u2265 0.05','< 0.05')))
table(fdrother$gloPlt05)

##Plot FDR based results
a<-ggplot(fdrother,aes(x=loc.jt.h2,y=glo.jt.h2,color=`distal P`)) + geom_point(cex=1,alpha=2/3)  + 
  xlab(expression("local h"^2)) + ylab(expression("distal h"^2)) + coord_cartesian(xlim=c(-0.05,1.05),ylim=c(-0.7,1.3)) + theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),legend.justification=c(1,1),legend.position=c(1,1))

##plot joint h2 estimates
local <- fdrother %>% select(loc.jt.h2,loc.jt.se,loc.jt.P,locPlt05)%>% arrange(loc.jt.h2) 
local <- local %>% mutate(loc.jt.h2=ifelse(is.na(loc.jt.h2),0,loc.jt.h2), loc.jt.se=ifelse(is.na(loc.jt.se),base::sample(loc.jt.se[is.na(loc.jt.se)==FALSE][1:100],size=length(loc.jt.se[is.na(loc.jt.se)==TRUE]),replace=TRUE),loc.jt.se))%>% arrange(loc.jt.h2) 
names(local) = c('h2','se','jt.P','Plt05')
data <- local %>% mutate(ymin =  h2 - 2 * se, ymax =  h2 + 2 * se )
cigt0 <- data$ymin>0
table(data$jt.P < 0.05,useNA='i')
sum(table(data$jt.P < 0.05,useNA='i'))
ptrue <-round(table(data$jt.P < 0.05,useNA='i')/sum(table(data$jt.P < 0.05,useNA='i')),3)*100
ptrue
meanh2<-round(mean(data$h2),3)
meanh2
meanse <- round(mean(data$se),3)
meanse
data <- mutate(data,P=factor(Plt05,labels=c("\u2265 0.05","< 0.05")),position=1:nrow(data))
my_grob2 = grobTree(textGrob(substitute(paste("% P < 0.05: ", m),list(m=ptrue[2])), x=0.05,  y=0.70, hjust=0,gp=gpar(fontsize=14)))

b<-ggplot(data,aes(x=position,y=h2,ymin=ymin, ymax=ymax, color=P) ) + geom_pointrange(col='gray')+geom_point()+ylab(expression("local h"^2)) + 
  xlab(expression("genes ordered by local h"^2))+theme_bw()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),legend.justification=c(0,1),legend.position=c(0,1))+
  annotation_custom(my_grob2)

global <- fdrother %>% select(gene,glo.jt.h2,glo.jt.se,glo.jt.P,gloPlt05) %>% arrange(glo.jt.h2) %>% mutate(glo.jt.h2=ifelse(is.na(glo.jt.h2),0,glo.jt.h2), glo.jt.se=ifelse(is.na(glo.jt.se),base::sample(glo.jt.se[is.na(glo.jt.se)==FALSE][1:1000],size=length(glo.jt.se[is.na(glo.jt.se)==TRUE]),replace=TRUE),glo.jt.se))%>% arrange(glo.jt.h2)  
names(global) = c('gene','h2','se','jt.P','Plt05')
data <- global %>% mutate(ymin =  h2 - 2 * se, ymax =  h2 + 2 * se)
table(data$jt.P < 0.05,useNA='i')
sum(table(data$jt.P < 0.05,useNA='i'))
ptrue <-signif(table(data$jt.P < 0.05,useNA='i')/sum(table(data$jt.P < 0.05,useNA='i')),3)*100
ptrue
meanh2<-round(mean(data$h2),3)
meanh2
meanse <- round(mean(data$se),3)
meanse
data <- mutate(data,P=factor(Plt05,labels=c("\u2265 0.05","< 0.05")),position=1:nrow(data))

glopriorlist <- dplyr::filter(data,Plt05==TRUE) %>% dplyr::select(gene)

my_grob2 = grobTree(textGrob(substitute(paste("% P < 0.05: ", m),list(m=ptrue[2])), x=0.05,  y=0.70, hjust=0,gp=gpar(fontsize=14)))

c<-ggplot(data,aes(x=position,y=h2,ymin=ymin, ymax=ymax, color=P) ) + geom_pointrange(col='gray')+geom_point()+ylab(expression("distal h"^2)) + 
  xlab(expression("genes ordered by distal h"^2))+theme_bw()+coord_cartesian(ylim=c(-1.05,1.55))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),legend.justification=c(0,1),legend.position=c(0,1))+
  annotation_custom(my_grob2)

##make same plots using all SNPs on other chrs for global
otherfile<-my.dir %&% 'expArch_DGN-WB_imputedGTs/DGN-WB.h2.all.models_FHSfdr0.05.Chr1-22_globalAll_reml-no-constrain.2015-12-15.txt'

other<-read.table(otherfile,header=T) ##all SNPs on non-gene chromosomes used to define global GRM
d  <- other %>% mutate(ymin =  glo.jt.h2 - 2 * glo.jt.se, ymax =  glo.jt.h2 + 2 * glo.jt.se )

other <- mutate(d, loc.jt.P = ifelse(loc.jt.h2 > 0, 2*pnorm(-abs(loc.jt.h2/loc.jt.se)),1), glo.jt.P = ifelse(glo.jt.h2 > 0, 2*pnorm(-abs(glo.jt.h2/glo.jt.se)),1)) %>% mutate(loc.jt.P=ifelse(is.na(loc.jt.P),1,loc.jt.P),glo.jt.P=ifelse(is.na(glo.jt.P),1,glo.jt.P)) %>%   mutate(locPlt05=loc.jt.P < 0.05,gloPlt05=glo.jt.P < 0.05) %>% mutate(`distal P`=factor(gloPlt05,labels=c('\u2265 0.05','< 0.05')))
table(other$gloPlt05)

aother<-ggplot(other,aes(x=loc.jt.h2,y=glo.jt.h2,color=`distal P`)) + geom_point(cex=1,alpha=2/3) + 
  xlab(expression("local h"^2)) + ylab(expression("distal h"^2)) + coord_cartesian(xlim=c(-0.05,1.05),ylim=c(-0.7,1.3)) + theme_bw() + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),legend.justification=c(1,1),legend.position=c(1,1))

local <- other %>% select(loc.jt.h2,loc.jt.se,loc.jt.P,locPlt05) %>% arrange(loc.jt.h2) %>% 
  mutate(loc.jt.h2=ifelse(is.na(loc.jt.h2),0,loc.jt.h2), loc.jt.se=ifelse(is.na(loc.jt.se),
  base::sample(loc.jt.se[is.na(loc.jt.se)==FALSE][1:100],size=length(loc.jt.se[is.na(loc.jt.se)==TRUE]),replace=TRUE),loc.jt.se))%>% arrange(loc.jt.h2)  
names(local) = c('h2','se','jt.P','Plt05')
data <- local %>% mutate(ymin =  h2 - 2 * se, ymax =  h2 + 2 * se )
table(data$jt.P < 0.05,useNA='i')
sum(table(data$jt.P < 0.05,useNA='i'))
ptrue <-round(table(data$jt.P < 0.05,useNA='i')/sum(table(data$jt.P < 0.05,useNA='i')),3)*100
ptrue
meanh2<-round(mean(data$h2),3)
meanh2
meanse <- round(mean(data$se),3)
meanse
data <- mutate(data,P=factor(Plt05,labels=c("\u2265 0.05","< 0.05")),position=1:nrow(data))
my_grob2 = grobTree(textGrob(substitute(paste("% P < 0.05: ", m),list(m=ptrue[2])), x=0.05,  y=0.70, hjust=0,gp=gpar(fontsize=14)))

bother<-ggplot(data,aes(x=position,y=h2,ymin=ymin, ymax=ymax, color=P) ) + geom_pointrange(col='gray')+geom_point()+ylab(expression("local h"^2)) + 
  xlab(expression("genes ordered by local h"^2))+theme_bw()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),legend.justification=c(0,1),legend.position=c(0,1))+
  annotation_custom(my_grob2)

global <- other %>% select(gene,glo.jt.h2,glo.jt.se,glo.jt.P,gloPlt05) %>% arrange(glo.jt.h2) %>% 
  mutate(glo.jt.h2=ifelse(is.na(glo.jt.h2),0,glo.jt.h2), glo.jt.se=ifelse(is.na(glo.jt.se),
  base::sample(glo.jt.se[is.na(glo.jt.se)==FALSE][1:1000],size=length(glo.jt.se[is.na(glo.jt.se)==TRUE]),replace=TRUE),glo.jt.se))%>% arrange(glo.jt.h2) 
names(global) = c('gene','h2','se','jt.P','Plt05')
data <- global %>% mutate(ymin =  h2 - 2 * se, ymax =  h2 + 2 * se )
table(data$jt.P < 0.05,useNA='i')
sum(table(data$jt.P < 0.05,useNA='i'))
ptrue <-signif(table(data$jt.P < 0.05,useNA='i')/sum(table(data$jt.P < 0.05,useNA='i')),3)*100
ptrue
meanh2<-round(mean(data$h2,useNA="T"),3)
meanh2
meanse <- round(mean(data$se),3)
meanse
data <- mutate(data,P=factor(Plt05,labels=c("\u2265 0.05","< 0.05")),position=1:nrow(data))

glolist <- dplyr::filter(data,Plt05==TRUE) %>% dplyr::select(gene)
table(glolist$gene %in% glopriorlist$gene)

my_grob2 = grobTree(textGrob(substitute(paste("% P < 0.05: ", m),list(m=ptrue[2])), x=0.05,  y=0.70, hjust=0,gp=gpar(fontsize=14)))

cother<-ggplot(data,aes(x=position,y=h2,ymin=ymin, ymax=ymax, color=P) ) + geom_pointrange(col='gray')+geom_point()+ylab(expression("distal h"^2)) + 
  xlab(expression("genes ordered by distal h"^2))+theme_bw()+coord_cartesian(ylim=c(-1.05,1.55))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),legend.justification=c(0,1),legend.position=c(0,1))+
  annotation_custom(my_grob2)

png(filename=fig.dir %&% "Fig-DGN-jt-h2-UNCONSTRAINED.png",width=720,height=960)
multiplot(aother+ggtitle('local = SNPs within 1Mb of gene\ndistal = SNPs on non-gene chrs\n') + theme(plot.title=element_text(face="bold")), bother, cother, a+ggtitle('local = SNPs within 1Mb of gene\ndistal = known eQTLs on non-gene chrs\n') + theme(plot.title=element_text(face="bold")),b, c,cols=2)
dev.off()

tiff(filename=fig.dir %&% "Fig-DGN-jt-h2-UNCONSTRAINED.tiff",width=720,height=960)
multiplot(aother+ggtitle('local = SNPs within 1Mb of gene\ndistal = SNPs on non-gene chrs\n') + theme(plot.title=element_text(face="bold")), bother, cother, a+ggtitle('local = SNPs within 1Mb of gene\ndistal = known eQTLs on non-gene chrs\n') + theme(plot.title=element_text(face="bold")),b, c,cols=2)
dev.off()