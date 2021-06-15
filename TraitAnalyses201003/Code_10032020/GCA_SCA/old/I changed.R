


setwd('/Users/maohuang/Desktop/Kelp/GCA_SCA/GCAPM.TS/')  # I changed
load("HcI.rda")
source('/Users/maohuang/Desktop/Kelp/GCA_SCA/GCAPM.TS/splitingSNPs.cI/input/parameters.R')

load(mm.file)  ### not sure where

SNPs <- HcI
NaNs <- matrix(NA,nrow=1, ncol=dim(SNPs)[2])
for(i in 1:dim(SNPs)[2])
{
  NaNs[i] <- length(which(is.na(SNPs[,i])))
}
   # created a NA matrix of 1 row, ncol(SNPs)

index.1 <- order(NaNs) 
NaNs2 <- NaNs[index.1]

NaN.freq <- unique(NaNs2)   
percentage.NaNs <- NaN.freq * 100 / dim(SNPs)[1]


write.table(NaN.freq, file='NaNs.freq.csv', sep=',',row.names=F, col.names=F)
write.table(percentage.NaNs, file='percentage.NaNs.csv', sep=',',row.names=F, col.names=F)



index.2 <- which(NaNs <= NaN.freq.i)
X <- SNPs[,index.2]

p <- colMeans(X,na.rm=T)/2    
p <- ifelse(p <= 0.5,p,1-p)


print( length( which( p >= prop.MAF.j ) ) )
index.3 <- which( p >= prop.MAF.j ) 
write.csv(X[,index.3], file='X.csv')
print( paste('NaNs_',NaN.freq.i,'_MAF_',prop.MAF.j,sep='') )




quit(save='no')








