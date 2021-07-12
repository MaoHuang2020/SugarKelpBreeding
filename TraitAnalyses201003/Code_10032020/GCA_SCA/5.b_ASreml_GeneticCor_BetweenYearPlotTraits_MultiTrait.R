#ASreml_Genetic_Cor_BetweenYear_Plot Level only


rm(list=ls())
#WD<-"/Users/maohuang/Desktop/Kelp/GCA_SCA/" # local
#datafdr<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/data/"

WD<-"/local/workdir/mh865/GCA_SCA/"
datafdr<-paste0(WD,"OneTime1920/data/")

###1. Input outCovComb4_dipOrder
load(paste0("/local/workdir/mh865/outCovComb/outCovComb4_Mix_Conden_0527_2021.Rdata"))

All_grm<-outCovComb4_MixOrder  # 866

Gencor_Yr<-function (data=dataNH,Trait_grm=Trait_grm3){
  mod4<-asreml(Trait~Year+line+block+popChk,
               random= ~ us(Year):vm(Crosses,Trait_grm3),
               data = data, maxiter=100, trace=TRUE)
  return(list(mod=mod4))
}
###2. Input dataNH

#### Plot level
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

### Plot Level----Add blades!!!!
traits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades") 
dataNH<-dataNHpi   ### !!!!! Only the plot traits

### !!!!! Change the corresponding traits for each data set too!!!!
#R
library(asreml)
asreml.license.activate()
#enter this code CDEA-HECC-CDAH-FIED

modSum<-NULL
cor<-NULL
for (j in 1:length(traits)){
  
  dataNH$Trait<-dataNH[,traits[j]]
  dataNH<-droplevels(dataNH[!is.na(dataNH$Trait),])
  
  if (traits[j]%in%c("dryWgtPerM","AshFDwPM")){
    dataNH$Trait2<-ifelse(dataNH$Year==2019,dataNH$Trait*sqrt(10),dataNH$Trait)  # Yr1 phenotype * sqrt(10)
    dataNH$Trait<-dataNH$Trait2
  }
  
  Trait_Crosses<-as.character(dataNH$Crosses)
  Trait_grm<-All_grm[rownames(All_grm)%in%Trait_Crosses,colnames(All_grm)%in%Trait_Crosses]
  print(dim(dataNH))
  print(sum(is.na(dataNH$Trait)))
  print(dim(Trait_grm))
  
  ### Adding the checks into the Trait_grm, all 1s in diagonal,all 0s for others
  data<-dataNH    
  data[data$popChk=="ES",]$Crosses
  droplevels(data[data$popChk=="ES",])$Crosses
  ChkCross<-unique(droplevels(data[!data$popChk=="ES",])$Crosses) # 33 plots of checks, 4 unique ones
  Col0<-matrix(0,nrow=nrow(Trait_grm),ncol=length(ChkCross))
  colnames(Col0)<-ChkCross
  
  Trait_grm2<-cbind(Trait_grm,Col0)
  
  Row0<-matrix(0,nrow=length(ChkCross),ncol=ncol(Trait_grm))
  Chk1<-diag(x=1,nrow=length(ChkCross),ncol=length(ChkCross))
  Row0_Chk1<-cbind(Row0,Chk1)
  rownames(Row0_Chk1)<-ChkCross
  
  Trait_grm3<-rbind(Trait_grm2,Row0_Chk1)
  
  print(nrow(Trait_grm3)==length(unique(droplevels(data$Crosses))))
  modSum[[j]]<-summary(Gencor_Yr(data=data,Trait_grm=Trait_grm3)$mod)$varcomp
  #covarianceA_B/sqrt(varianceA*varianceB)
  
  cor<-c(cor,modSum[[j]][,"component"][2]/sqrt(modSum[[j]][,"component"][1]*modSum[[j]][,"component"][3]))
  
}

names(cor)<-traits

#### !!! Between Yr, Same Trait Plot, Gen_cor
Plotcor<-cor
write.csv(cor,"Genetic_Correlations_BetweenYears_PlotLevel_05272021.csv")
save(modSum,file="Genetic_Correlations_BetweenYears_PlotLevel_moddelSummary_05272021.Rdata")

