###############################################  
#### Making BLUEs within Year and between Years
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

#### Individual level to get their experimental factors.
#### The dataNHpi is already filtered for their phenotypic PHOTO SCORE (>2)
dataNHim<-dataNHimboth_C
dataNHim$line<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "line")
dataNHim$block<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "block")
dataNHim$Year<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Year")
dataNHim$popChk<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "popChk")
dataNHim$Crosses<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Crosses")
dataNHim$PhotoScore<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "PhotoScore")
dataNHim<-dataNHim[which(dataNHim$PhotoScore >1),]  # 3969 rows with PhotoScore >1
  str(dataNHim)
  dim(dataNHim)

traits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")

for (col in traits){
  dataNHim[,col]<-as.numeric(dataNHim[,col])}

dataNH<-dataNHim   ### !!!!!
#######
#######

traits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades") 
dataNH<-dataNHpi   ### !!!!!

#######
library(sommer)
  dim(dataNH)
###################
### Both Years
dataNH<-dataNH
  
###################
### Within Year
  # 2020_273 is the same as 2020_219
  comX<-c("SL18-OI-15-FG1xSA18-CB-2-MG3","SL18-LD-13-FG2xSL18-LD-13-MG2","SL18-LD-2-FG2xSL18-SF-19-MG1","SL18-NL-2-FG3xSA18-CB-4-MG1","SL18-NL-3-FG1xSA18-CB-7-MG2")
Yr<-2020  ###!!!!
#
Yr<-2019  ###!!!!
WithinYear=TRUE
dataNH<-dataNH[dataNH$Year==Yr,]
dataNH<-dataNH[order(dataNH$plotNo),]


  dim(dataNH)
if (WithinYear==FALSE){
  ### Both Years
  NumCrosses<-length(unique(dataNH$Crosses))
  dBLUPs<-BLUPs<-matrix(nrow=NumCrosses,ncol=length(traits))
  H2<-NULL
  convInfo<-vector()
  
for (j in 1:length(traits)){
  Coltrait<-traits[j]
  dataNH$Trait<-dataNH[,Coltrait]  ### !!!!! TraitBLUE
  data<-dataNH
  # #for (Year in c(2019,2020)){
  #   i=Year-2018
  #   data<-dataNH[dataNH$Year==Year,]
  data<-droplevels(data)
  
  tmp.mod <- mmer(fixed = Trait ~ popChk+Year, 
                  random = ~block + line + block:line+Crosses, 
                  data = data)
  
  PEV <- diag(tmp.mod$PevU$`Crosses`$Trait)
  varG <- as.numeric(tmp.mod$sigma$`Crosses`)
  tmp.BLUPs <- tmp.mod$U$`Crosses`$Trait  ##blups
  names(tmp.BLUPs)<-stringr::str_replace(names(tmp.BLUPs),"Crosses","")
  tmp.dBLUPs<-tmp.BLUPs / (1 - PEV/varG)  ##deregressed blups checks will have "inf" values
  print(identical(names(tmp.BLUPs),names(tmp.dBLUPs)))
  CheckCrosses<-unique(droplevels(data[!data$popChk=="ES",]$Crosses))  ## 4 Crosses
  CheckCrosses
  tmp.BLUPs[names(tmp.BLUPs)%in%CheckCrosses]<-NA
  tmp.dBLUPs[names(tmp.dBLUPs)%in%CheckCrosses]<-NA
  
  
  BLUPs[,j] <- tmp.BLUPs  
  dBLUPs[,j] <- tmp.dBLUPs  
  
  #H2 <- rbind(H2, pin(tmp.mod, h2 ~ V2 / ( V2 + V4)))
  H2<-rbind(H2,h2.fun(tmp.mod,data=data,gTerm="Crosses"))
  convInfo[j] <- tmp.mod$convergence
  }
  
} else if(WithinYear==TRUE){
  ### Within Year
  NumCrosses<-length(unique(dataNH$Crosses))
  dBLUPs<-BLUPs<-matrix(nrow=NumCrosses,ncol=length(traits))
  H2<-NULL
  convInfo<-vector()

  for (j in 1:length(traits)){
    Coltrait<-traits[j]
    dataNH$Trait<-dataNH[,Coltrait]  ### !!!!! TraitBLUE
    data<-dataNH
    # #for (Year in c(2019,2020)){
    #   i=Year-2018
    #   data<-dataNH[dataNH$Year==Year,]
    data<-droplevels(data)
    
    tmp.mod <- mmer(fixed = Trait ~ popChk, 
                    random = ~block + line + block:line+Crosses, 
                    data = data)
    
    PEV <- diag(tmp.mod$PevU$`Crosses`$Trait)
    varG <- as.numeric(tmp.mod$sigma$`Crosses`)
    tmp.BLUPs <- tmp.mod$U$`Crosses`$Trait  ##blups
    names(tmp.BLUPs)<-stringr::str_replace(names(tmp.BLUPs),"Crosses","")
    tmp.dBLUPs<-tmp.BLUPs / (1 - PEV/varG)  ##deregressed blups checks will have "inf" values
    print(identical(names(tmp.BLUPs),names(tmp.dBLUPs)))
    CheckCrosses<-unique(droplevels(data[!data$popChk=="ES",]$Crosses))  ## 4 Crosses
    CheckCrosses
    tmp.BLUPs[names(tmp.BLUPs)%in%CheckCrosses]<-NA
    tmp.dBLUPs[names(tmp.dBLUPs)%in%CheckCrosses]<-NA
    
    
    BLUPs[,j] <- tmp.BLUPs  
    dBLUPs[,j] <- tmp.dBLUPs  
    
    #H2 <- rbind(H2, pin(tmp.mod, h2 ~ V2 / ( V2 + V4)))
    H2<-rbind(H2,h2.fun(tmp.mod,data=data,gTerm="Crosses"))
    convInfo[j] <- tmp.mod$convergence
  }
  
}
################################

rownames(H2)<-traits
  H2
  convInfo
save(H2,convInfo,file=paste0(datafdr,Yr,"_Indiv_sommer_h2.fuc_convergence.Rdata")) #/Plots_

colnames(BLUPs)<-colnames(dBLUPs)<-traits
rownames(BLUPs)<-rownames(dBLUPs)<-names(tmp.dBLUPs)
  head(BLUPs)
  head(dBLUPs)
  dim(BLUPs)
  dim(dBLUPs)

BLUPs<-BLUPs[!rownames(BLUPs)%in%CheckCrosses,]     #RM the checks crosses
dBLUPs<-dBLUPs[!rownames(dBLUPs)%in%CheckCrosses,]    


####
Yr19_dBLUPs<-as.data.frame(dBLUPs)  ##!!!!
Yr19_dBLUPs$Year<-Yr
  dim(Yr19_dBLUPs)
Yr19_dBLUPs$plotNo<-expss::vlookup(rownames(Yr19_dBLUPs),dict=dataNH[dataNHpi$Year==Yr,],lookup_column = "Crosses",result_column = "plotNo")
Yr19_dBLUPs$Crosses<-rownames(Yr19_dBLUPs)
rownames(Yr19_dBLUPs)<-Yr19_dBLUPs$plotNo
  Yr19_dBLUPs[Yr19_dBLUPs$Crosses=="SL18-LD-13-FG2xSL18-LD-13-MG2",]
###
Yr20_dBLUPs<-as.data.frame(dBLUPs)  ##!!!!
Yr20_dBLUPs$Year<-Yr
  dim(Yr20_dBLUPs)
Yr20_dBLUPs$plotNo<-expss::vlookup(rownames(Yr20_dBLUPs),dict=dataNHpi[dataNHpi$Year==Yr,],lookup_column = "Crosses",result_column = "plotNo")
Yr20_dBLUPs$Crosses<-rownames(Yr20_dBLUPs)
rownames(Yr20_dBLUPs)<-Yr20_dBLUPs$plotNo  
  Yr20_dBLUPs[Yr20_dBLUPs$Crosses%in%comX,]
  Yr20_dBLUPs[Yr20_dBLUPs$Crosses=="SL18-LD-13-FG2xSL18-LD-13-MG2",]

#### Within Yr
WithinYr_Plot_dBLUPs<-rbind(Yr19_dBLUPs,Yr20_dBLUPs) ##!!!
save(WithinYr_Plot_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_plotlevel_WithinYear_AddBD.Rdata"))
  dim(WithinYr_Plot_dBLUPs)  

WithinYr_Indi_dBLUPs<-rbind(Yr19_dBLUPs,Yr20_dBLUPs)
  dim(WithinYr_Indi_dBLUPs)
save(WithinYr_Indi_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_Indivlevel_WithinYear.Rdata"))

WithinYr_Both_dBLUPs<-merge(WithinYr_Plot_dBLUPs,WithinYr_Indi_dBLUPs,by.x="row.names",by.y="row.names",all.x=TRUE)
save(WithinYr_Both_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata"))



#### Two Years
Plot_dBLUPs<-dBLUPs  ##!!!
save(BLUPs,dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_plotlevel_OverTwoYears_AddBD.Rdata"))

Indi_dBLUPs<-dBLUPs  ##!!!
save(BLUPs,dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_Individuallevel_OverTwoYears.Rdata"))

  dim(Indi_dBLUPs)
  dim(Plot_dBLUPs)
sum(rownames(Indi_dBLUPs)%in%rownames(Plot_dBLUPs))
Both_dBLUPs<-merge(Plot_dBLUPs,Indi_dBLUPs,by.x="row.names",by.y="row.names",all.x=TRUE)  
  dim(Both_dBLUPs)
  head(Both_dBLUPs)  
  
save(Both_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_overTwoYears_AddBD.Rdata"))

## Compare the BLUPs from BothYears vs WithinYear (cor ranged= 0.97-0.98)
both<-merge(Plot_dBLUPs,WithinYr_Plot_dBLUPs,by.x="row.names",by.y="row.names",all.x=TRUE)
for (i in 1:5){
  print(cor(both[,(i+1)],both[,(i+1+5)],use="complete"))
}

#######
#### Summarize the GS TP size--- BothYr
#datafdr<-paste0(WD,"OneTime1920/data/")  ## In terminal
datafdr<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/data/"
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_overTwoYears_AddBD.Rdata")) ##!!!

size<-NULL
traits<-colnames(Both_dBLUPs)[-1]
for (t in 1:length(traits)){
  size<-c(size,sum(!is.na(Both_dBLUPs[,traits[t]])|!is.nan(Both_dBLUPs[,traits[t]])))
}
names(size)<-traits
#######
#######

load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata"))


WithinYrSize<-NULL
for (Yr in c(2019,2020)){
  dataf<-droplevels(WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==Yr,])
  rmcols<-c("Row.names","Year.x","plotNo.x","Crosses.x","Year.y","plotNo.y","Crosses.y")
  traits<-colnames(dataf)[!colnames(dataf)%in%rmcols]
  size<-NULL
  for (t in 1:length(traits)){
    size<-c(size,sum(!is.na(dataf[,traits[t]])|!is.nan(dataf[,traits[t]])))
  }
  names(size)<-traits
  WithinYrSize<-rbind(WithinYrSize,size)
}

