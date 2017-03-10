####################################
# Calc adj PEER factors using GCTA #
# adj PF = residuals from BLUP     #
# by Heather E. Wheeler            #
# 2017-03-06                       #
####################################

"%&%" = function(a,b) paste(a,b,sep="")
library(dplyr)

my.dir <- "/group/im-lab/nas40t2/hwheeler/cross-tissue/"
dgn.dir <- my.dir %&% "expArch_DGN-WB_imputedGTs/"
trans.dir <- dgn.dir %&% "transPCA/"
tis <- "DGN-WB"
data <- Sys.Date()
args <- commandArgs(trailingOnly=T)

Nk <- args[1] ##number of peer factors
pf.file <- tis %&% "." %&% Nk %&% ".PEER.factors.2017-03-06.txt"

### Output matrix and data.frame
hsq.mat <- matrix(0,nrow=as.numeric(Nk),ncol=4)
colnames(hsq.mat) <- c("PF","h2","se","p")

pf <- read.table(trans.dir %&% pf.file)
adj.pf <- pf[,1:2]

for(i in 1:Nk){
##try --reml and --reml-no-constrain
  rungcta <- "gcta64 --reml-no-constrain --grm " %&% dgn.dir %&% 
  "dgn-grms/DGN.global_Chr1-22 --pheno " %&% trans.dir %&% pf.file %&% 
  " --mpheno " %&% i %&% " --reml-pred-rand --out tmp." %&% tis %&% Nk %&%
  " --reml-maxit 1000"
  system(rungcta)
  #check for .hsq file, collect results
  if(file.exists("tmp." %&% tis %&% Nk %&% ".hsq")==TRUE){
    hsq <- scan("tmp." %&% tis %&% Nk %&% ".hsq","character")
    res <- c(i, hsq[14], hsq[15], hsq[25])
  }else{ #gcta did not coverge or Error: the information matrix is not invertible.
    res <- c(i, NA, NA, NA)
  }
  hsq.mat[i,] <- res
  system("rm tmp." %&% tis %&% Nk %&% ".hsq")
  #check for .indi.blp file, collect residuals (col 6)
  if(file.exists("tmp." %&% tis %&% Nk %&% ".indi.blp")==TRUE){
    blup <- read.table("tmp." %&% tis %&% Nk %&% ".indi.blp")
    adj.pf <- cbind(adj.pf,blup[,6])
  }else{ #gcta did not converge, thus keep original PF
    adj.pf <- cbind(adj.pf,pf[i+2])
  }
  system("rm tmp." %&% tis %&% Nk %&% ".indi.blp")
}

write.table(hsq.mat, file = trans.dir %&% "hsq_unconstrained_" %&% pf.file, quote=F, row.names=F)
write.table(adj.pf, file = trans.dir %&% "adj_unconstrained_" %&% pf.file, quote=F, row.names=F, col.names=F)
