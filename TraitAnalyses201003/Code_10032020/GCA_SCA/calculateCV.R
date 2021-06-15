## CV of checks
rm(list=ls())
WD<-"/Users/maohuang/Desktop/Kelp/GCA_SCA/" # local
datafdr<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/data/"

# load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot-- OLD
load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23_UpdateAsh_0309_2021.rdata"))  ## Plot -- Updated Ash
load(paste0(datafdr,"dataNHim_withChk_3_sets_PhotoScore0123.rdata"))  ## Individual

dataNHpi<-dataNHpiBoth_C  ##!!!!!

exptlSP <- as.character(dataNHpi$popChk) == "ES"
dataNHpi$entry <- as.character(dataNHpi$popChk)
dataNHpi$entry[exptlSP] <- as.character(dataNHpi$plotNo[exptlSP])
dataNHpi$group <- as.factor(ifelse(exptlSP, 1, 0))

for (col in c( "line", "block","popChk","group","entry","Year")){ 
  dataNHpi[,col] <- factor(dataNHpi[,col])}
dataNHpi<-droplevels(dataNHpi)

traits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades") 
dataNH<-dataNHpi   ### !!!!!


dataCHK<-droplevels(dataNH[dataNH$popChk%in%c("Z1","Z2","Z3","Z4"),])
  head(dataCHK)
  dim(dataCHK)  
  dim(dataNHpi)  
  
Y1CHK<-dataCHK[dataCHK$Year==2019,]  
  dim(Y1CHK)
Y2CHK<-dataCHK[dataCHK$Year==2020,]  
  dim(Y2CHK)
  
cv<-function(x){
  (sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE))*100
}  

cv1<-cv(Y1CHK$dryWgtPerM)
cv2<-cv(Y1CHK$AshFDwPM)

cv3<-cv(Y2CHK$dryWgtPerM)
cv4<-cv(Y2CHK$AshFDwPM)

c(cv1,cv2,cv3,cv4)
