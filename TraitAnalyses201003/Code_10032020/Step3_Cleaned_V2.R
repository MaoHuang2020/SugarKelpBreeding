###### Both Years
###################################
rm(list=ls())

library(here)
here()
datafdr<-here("TraitAnalyses201003/data/")
wd<-here("TraitAnalyses201003/UpdateAsh/") ## !!
setwd(wd)
  getwd()
  
#here::i_am("Code_10032020/Step3_Cleaned_V2.R") 
  ls()

library(rrBLUP)
# 
# load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHpi_withChk_3_sets_PhotoScore23.rdata")   ## Plot
# load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi
  

  ####
  #load(paste0("hMat_PedNH_CCmat_fndrMrkData_",yr,"_PhotoScore23.rdata"))
  #load(paste0("hMat_PedNH_CCmat_fndrMrkData_Both_PhotoScore23_NoSGP.rdata"))
  # biphasicPedNH<-read.csv(here("TraitAnalyses201003/ReorderPedigree","Ped_in_Order_866_Individuals_Fndr_New_Order_0116_2021.csv"),sep=",",header=TRUE,row.names=1)
  # biphasicPedNH[1:4,]
  # spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))    ### Progeny sps, col 2 and 3 are all values
  # nSp <- length(spRows)
  # nSp #245
  # 
  
load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23_UpdateAsh_0309_2021.rdata"))
load(paste0(datafdr,"dataNHim_withChk_3_sets_PhotoScore0123.rdata"))  

### !!!!Both
dataNHpi<-dataNHpiBoth_C  
dataNHim<-dataNHimboth_C
yr<-"Both"
  ls()
head(dataNHpiBoth_C)
bothYr<-TRUE


# ### !!!!2019
dataNHpi<-dataNHpi19_C
dataNHim<-dataNHim19_C
yr<-"2019"
bothYr<-FALSE  # use the within year model, not both years' data model
# 

### !!!!2020
dataNHpi<-dataNHpi20_C
dataNHim<-dataNHim20_C
yr<-"2020"
bothYr<-FALSE


dataNHpi$densityBlades<-ifelse(dataNHpi$densityBlades==0,NA,dataNHpi$densityBlades)  # densityblades as 0 then NA

dataNHpi$popChk <- ifelse(substr(dataNHpi$plotNo, 1, 1) == "Z", substr(dataNHpi$plotNo, 1, 2), "ES")  # Checks VS ES

dataNHpi$withinLoc <- ifelse(as.vector(dataNHpi$femaParLoc) == as.vector(dataNHpi$maleParLoc), 1, 0) # WithinLoc is 1
dataNHpi$development[dataNHpi$development=="#N/A"] <-NA

# Experimental design variables should be factors
for (col in c( "Year", "plotNo","Region","GrowDays","femaPar", "femaParLoc", "malePar", "maleParLoc", "block", "line","popChk")) 
  dataNHpi[,col] <- factor(dataNHpi[,col])
  nlevels(dataNHpi$plotNo)
  nlevels(dataNHpi$Year)
  unique(dataNHpi$Year)
# Make data cols numeric. Enforce data is numeric
for (col in c("wetWgtPlot", "lengthPlot", "wetWgtPerM","percDryWgt",  "dryWgtPerM","densityBlades","AshFreedryWgtPerM","AshFDwPM")) 
  dataNHpi[,col] <- as.numeric(dataNHpi[,col])

# If percentDryWeigth is NA, then the plot WetWeightPerM and DryWeightperM should be set at NA, WetWeigthPerM may be just rope
dataNHpi[is.na(dataNHpi$percDryWgt), c("wetWgtPerM", "dryWgtPerM")] <- NA
keepRows <- !is.na(dataNHpi$percDryWgt)
  sum(keepRows)  # 281

fndrF1<-strsplit(as.character(dataNHpi$femaPar), split="-", fixed=T)
fndrM1<-strsplit(as.character(dataNHpi$malePar),split="-",fixed=T)
fndrF<-sapply(fndrF1, function(vec) paste(vec[1:3], collapse="-"))
fndrM<-sapply(fndrM1, function(vec) paste(vec[1:3],collapse="-"))

isSelf<-ifelse(fndrF==fndrM,1,0)
dataNHpi$isSelf <- isSelf 
  str(dataNHpi)

# Experimental design variables should be factors
for (col in c("Year","plotNo","GrowDays", "femaPar", "femaParLoc", "malePar", "maleParLoc", "line", "block", "date", "popChk", "withinLoc", "isSelf")) 
  dataNHpi[,col] <- as.factor(dataNHpi[,col])

### aMat has fndr+GP+ all the plot SPs in it.
### sort 

dataNHpi_RMchk<-dataNHpi[!dataNHpi$crossID=="Check",] # Subset without chk levels
dataNHpi_RMchk$Crosses<-as.factor(as.character(dataNHpi_RMchk$Crosses)) # This RMed the check levels
  dim(dataNHpi_RMchk)   # 250
  str(dataNHpi_RMchk)
nrowplot<-nrow(dataNHpi_RMchk)
nrowchk<-nrow(dataNHpi) - nrowplot


##{{}}
##{{}}
#Dip
load(here("TraitAnalyses201003/data","outCovComb_dip_0116_2021.Rdata"))
  outCovComb<-outCovComb4_dipOrder
diphap<-"dip"


# #Hap
# load(here("TraitAnalyses201003/data","hMat_hap_0116_2021.Rdata"))
# #load(here("TraitAnalyses201003/Making_haploid_CovComb","outCovComb4_hap_Conden_0116_2021.Rdata")) # this outCovComb4_Hapconden=hMat_hap
# outCovComb<-hMat_hap
# diphap<-"hap"
## {}
## {}
  
  dim(outCovComb)
  outCovComb[1:4,1:4]
phenoNamesFact<-factor(dataNHpi_RMchk$Crosses,levels=rownames(outCovComb))   #colnames are sorted alphabetically for the outCovComb 
msZ0<-model.matrix(~-1 +phenoNamesFact,data=dataNHpi_RMchk)  

#### Does this matter??
identical(colnames(outCovComb)[544:787],levels(dataNHpi_RMchk$Crosses)) ####this= FALSE !!!Ordering of crosses in levels dataNHPi differs from outCovComb1

## 3. chks rows
chkSpMat<-matrix(0, nrowchk, nrow(outCovComb))
  dim(chkSpMat)   # 33 x 866  # 2019, 17 x 866  # 2020, 16 x 866
#  
msZ<-rbind(msZ0,chkSpMat)
  dim(msZ)   # 283 x 866   # 2019, 139 x 866  # 2020, 144x866

  
# Calculate heritability from the mixed.solve output
heritability <- function(msOut){
  return(msOut$Vu / (msOut$Vu + msOut$Ve))
}

Vu<-function(msOut){
 msOut$Vu
}

ErrVar<-function(msOut){
 msOut$Ve
}

#write.csv(dataNHpi,paste0("dataNHpi_Last_Used_in_Model_",yr,".csv"))

########### 6. use hMat as the relationship matrx
#### If or not chks could be included, not for ash, if this is within year or between year analysis

RunHeritability<-function(BothYear=TRUE,UseChk=FALSE,Trait=dataNHpi$AshFDwPM,dataNHpi=dataNHpi,hMat=outCovComb){
  dataNHpi<-droplevels(dataNHpi)

if (UseChk==FALSE){
  if (BothYear==TRUE){
    #Both Years !!!!!
    msX1 <- model.matrix( ~ line%in%Year+block%in%Year+Year, data=dataNHpi)  ### !!!! No popChk because no checks for ashes
    msX1 <- msX1[, apply(msX1, 2, function(v) !all(v == 0))]
  }else{
    #Within Year !!!! 2019 2020
    msX1 <- model.matrix( ~ line+block, data=dataNHpi)
    msX1 <- msX1[, apply(msX1, 2, function(v) !all(v == 0))]
    
  }
  library(Matrix)  
  print(rankMatrix(msX1))
  print(qr(msX1)$rank)  #Both, 27;  2019, 12
  
  msOutTrait<-mixed.solve(y=Trait,Z=msZ, K=hMat, X=msX1, SE=T)
  #msOutAshOnly<-mixed.solve(y=dataNHpi$Ash,Z=msZ,K=hMat,X=msX1,SE=T)
  #msOutAshFreeDW<-mixed.solve(y=dataNHpi$AshFreedryWgtPerM,Z=msZ, K=hMat, X=msX1, SE=T)
} else if(UseChk==TRUE){
  
  if (BothYear==TRUE){
    ###!!!!!!!!!!!!!!!!!! Both Years !!!!!!!!!!!
    msX <- model.matrix( ~ line%in%Year+block%in%Year+Year+popChk, data=dataNHpi)  ### !!!! has popChk for the traits
    msX <- msX[, apply(msX, 2, function(v) !all(v == 0))]
    
  }else{
    # Within Year !!!! 2019 2020
    msX <- model.matrix( ~ line+block+popChk, data=dataNHpi)  
    msX <- msX[, apply(msX, 2, function(v) !all(v == 0))]
  }
  library(Matrix)
   print(rankMatrix(msX))
   print(qr(msX)$rank)  #31
  
  msOutTrait <- mixed.solve(y=Trait, Z=msZ, K=hMat, X=msX, SE=T)
  #msOutWWPh <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msX, SE=T)
  #msOutDWPMh <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msX, SE=T)
  #msOutPDWh <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msX, SE=T)
}
  return(list(msOutTrait=msOutTrait))
}

#{}

### log transformed DwPM and WWpM

msOutAshFDwPM<-RunHeritability(BothYear=bothYr,UseChk=FALSE,Trait=dataNHpi$AshFDwPM,dataNHpi=dataNHpi,hMat=outCovComb)$msOutTrait
msOutAshOnly<-RunHeritability(BothYear=bothYr,UseChk=FALSE,Trait=dataNHpi$Ash,dataNHpi=dataNHpi,hMat=outCovComb)$msOutTrait

msOutDWPMh <-RunHeritability(BothYear=bothYr,UseChk=TRUE,Trait=log(dataNHpi$dryWgtPerM+1),dataNHpi=dataNHpi,hMat=outCovComb)$msOutTrait

msOutCN<-RunHeritability(BothYear=bothYr,UseChk=FALSE,Trait=dataNHpi$C_N_Combine,dataNHpi=dataNHpi,hMat=outCovComb)$msOutTrait

Scott<-c(heritability(msOutAshFDwPM),heritability(msOutAshOnly),heritability(msOutDWPMh))
names(Scott)<-c("AFDW/M","%Ash","DW/M")
  Scott

  ## Making the plot for Scott, 0309_2021

barplot(Scott,main="Trait heritability",cex.main=2,cex.names=2.0)
  
par(mar=par("mar")+c(4,0,0,0))
barplot(Scott, ylab="Heritability", las=2, cex.axis=1.0, cex.lab=1.3, cex.names=1.0)

msOutAshFreeDW<-RunHeritability(BothYear=bothYr,UseChk=FALSE,Trait=dataNHpi$AshFreedryWgtPerM,dataNHpi=dataNHpi,hMat=outCovComb)$msOutTrait
msOutWWPh <- RunHeritability(BothYear=bothYr,UseChk=TRUE,Trait=log(dataNHpi$wetWgtPlot+1),dataNHpi=dataNHpi,hMat=outCovComb)$msOutTrait
msOutPDWh <- RunHeritability(BothYear=bothYr,UseChk=TRUE,Trait=dataNHpi$percDryWgt,dataNHpi=dataNHpi,hMat=outCovComb)$msOutTrait


msOutDBh <- RunHeritability(BothYear=bothYr,UseChk=TRUE,Trait=dataNHpi$densityBlades,dataNHpi=dataNHpi,hMat=outCovComb)$msOutTrait


  heritability(msOutAshFreeDW) #both 0.561 #2019 #2020 0.7995067   # Dave cal
  heritability(msOutAshFDwPM) #both 0.461 #2019 #2020 0.5173875    # my cal
  heritability(msOutAshOnly)  #both 0.070 
  
traitsnames<-c("AshFDwPM","AshOnly", "dryWgtPerM","C:N ratio","wetWgtPlot", "percDryWgt","plotDensity")  
h2hMat <- c(heritability(msOutAshFDwPM),
            heritability(msOutAshOnly),
            heritability(msOutDWPMh),
            heritability(msOutCN),
            heritability(msOutWWPh),
            heritability(msOutPDWh),
            heritability(msOutDBh))
  options(scipen = 999) # turn off scientific
Vus<-c(Vu(msOutAshFDwPM),
       Vu(msOutAshOnly),
       Vu(msOutDWPMh),
       Vu(msOutCN),
       Vu(msOutWWPh),  
       Vu(msOutPDWh),
       Vu(msOutDBh)) 

ErrVars<-c(ErrVar(msOutAshFDwPM),
           ErrVar(msOutAshOnly),
           ErrVar(msOutCN),
           ErrVar(msOutDWPMh), 
           ErrVar(msOutWWPh),
           ErrVar(msOutPDWh),
           ErrVar(msOutDBh))
names(h2hMat) <-names(Vus)<-names(ErrVars)<-traitsnames
  Vus
  ErrVars
  round(h2hMat,3) 
  
H2Table<-rbind(h2hMat,Vus,ErrVars)
write.csv(H2Table,paste0(datafdr,"H2_Plot_Level_",yr,"years","_usingMatrix_",diphap,".csv"))
  

allBLUPs <- cbind(msOutAshFDwPM$u,msOutAshOnly$u,msOutDWPMh$u,msOutWWPh$u,  msOutPDWh$u,msOutDBh$u)
  dim(allBLUPs)
colnames(allBLUPs) <- c("AshFDwPM","AshOnly","DWpM","WWP","PDW","BD")
  head(allBLUPs)
#write.csv(allBLUPs,paste0(wd,"allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_0309_2021_dip.csv"))  
write.csv(allBLUPs,paste0(datafdr,"allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_","0315_2021_",diphap,".csv"))  


############# Individual BLUPs
######### Use these from plot level info
if (BothYear==TRUE){
  ###!!!!!!!!!!!!!!!!!! Both Years !!!!!!!!!!!
  msX <- model.matrix( ~ line%in%Year+block%in%Year+Year+popChk, data=dataNHpi)  ### !!!! has popChk for the traits
  msX <- msX[, apply(msX, 2, function(v) !all(v == 0))]
  
}else{
  # Within Year !!!! 2019 2020
  msX <- model.matrix( ~ line+block+popChk, data=dataNHpi)  
  msX <- msX[, apply(msX, 2, function(v) !all(v == 0))]
}

hMat=outCovComb
#############


for (col in c( "Year", "plotNo")) 
  dataNHim[,col] <- factor(dataNHim[,col])

  tail(dataNHim)
  dim(dataNHim)
  tail(dataNHpi)

for (col in c("bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter"))  
  dataNHim[,col] <- as.numeric(dataNHim[,col])

dataNHim2<-subset(dataNHim,dataNHim$plotNo%in%dataNHpi$plotNo) #!!!! Each plotNo is unique
  dim(dataNHim2)
  dim(dataNHim)
dataNHim2$plotNo<-factor(dataNHim2$plotNo)
  nlevels(dataNHim2$plotNo) 
  nlevels(dataNHpi$plotNo)
  nlevels(dataNHim$plotNo)

levels(dataNHim2$plotNo) <- levels(dataNHpi$plotNo)
dataNHim0<-dataNHim
dataNHim<-dataNHim2

### No Need to be  changed to Crosses, individual connects with plots
msXim <- model.matrix( ~ -1 + plotNo, data=dataNHim)  # nrow=nrow(dataNHim), ncol=number of plotNo=nlevels(plotNo)

#rownames(msXim)<-rownames(dataNHim)
  dim(msXim)
  dim(dataNHim)

msZim <- emZsp <- msXim %*% msZ
#msXim <- msXim %*% msXdb
msXim<-msXim%*%msX

  dim(msX)
  dim(msXim)
  dim(msZim)
  dim(dataNHim)

msXim<-msXim[, apply(msXim, 2, function(v) !all(v == 0))]

msOutBLh <- mixed.solve(y=log(dataNHim$bladeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBMWh <- mixed.solve(y=log(dataNHim$bladeMaxWidth+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBTh <- mixed.solve(y=log(dataNHim$bladeThickness+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutSLh <- mixed.solve(y=log(dataNHim$stipeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutSDh <- mixed.solve(y=log(dataNHim$stipeDiameter+1), Z=msZim, K=hMat, X=msXim, SE=T)
h2hMatim <- c(heritability(msOutBLh), heritability(msOutBMWh), heritability(msOutBTh), heritability(msOutSLh), heritability(msOutSDh)) #
names(h2hMatim) <- c("bladeLength", "bladeMaxWidth",  "bladeThickness", "stipeLength", "stipeDiameter") 

  round(h2hMatim,3)
  names(h2hMat) 
allh2h<-c(h2hMat,h2hMatim)
  allh2h
write.csv(allh2h,paste0(datafdr,"H2_Individual_Level_",yr,"years_usingMatrix_",diphap,".csv"))

allBLUPs_addIndiv <- cbind(msOutAshFDwPM$u,msOutAshOnly$u,msOutDWPMh$u,msOutWWPh$u,  msOutPDWh$u, msOutDBh$u, msOutBLh$u, msOutBMWh$u,  msOutBTh$u, msOutSLh$u, msOutSDh$u)
  head(allBLUPs_addIndiv)
colnames(allBLUPs_addIndiv) <- c("AshFDwPM","AshOnly","DWpM","WWP",  "PDW", "BDns", "BLen", "BMax",  "BThk", "SLen", "SDia")
write.csv(allBLUPs_addIndiv,paste0(datafdr,"allBLUPs_Plots+Individuals_withSGP_866_AddfndrsMrkData_","0315_2021_",diphap,".csv"))

Y<-dataNHpi
plot<-ggplot(data=Y,aes(densityBlades,Year))+
  geom_point(aes(color=as.factor(Year)))+ 
  geom_line(aes(group=as.factor(Crosses)))
print(plot)


# traitsherit1<-c("Ash","C_N_Combine","AshFDwPM")
# traitsherit2<-c("dryWgtPerM","wetWgtPlot","percDryWgt","densityBlades")
# 
# Vus<-NULL
# ErrVars<-NULL
# Herit<-NULL
# for (traits in traitsherit1){
#   Trait<-dataNHpi[,traits]
#   msOut<-RunHeritability(BothYear=TRUE,UseChk=FALSE,Trait=Trait,dataNHpi=dataNHpi,hMat=outCovComb)$msOutTrait
#   Vus<-c(Vus,Vu(msOut))
#   ErrVars<-c(ErrVars,ErrVar(msOut))
#   Herit<-c(Herit,heritability(msOut))
# }
# 
# H1all<-rbind(Herit,Vus,ErrVars)
# colnames(H1all)<-traitsherit1
#   H1all
# 
# Vus<-NULL
# ErrVars<-NULL
# Herit<-NULL
# for (traits in traitsherit2){
#   Trait<-dataNHpi[,traits]
#   msOut<-RunHeritability(BothYear=TRUE,UseChk=TRUE,Trait=Trait,dataNHpi=dataNHpi,hMat=outCovComb)$msOutTrait
#   Vus<-c(Vus,Vu(msOut))
#   ErrVars<-c(ErrVars,ErrVar(msOut))
#   Herit<-c(Herit,heritability(msOut))
# }

# H2all<-rbind(Herit,Vus,ErrVars)
# colnames(H2all)<-traitsherit2
#   H2all
# 
#   H2all





########### Plotting BVs correlation, and BVs by generations
#rownames(allBLUPs) is the same as that in u, and hMat

allBLUPs <- cbind(msOutAshFDwPM$u,msOutAshOnly$u,msOutDWPMh$u,msOutWWPh$u,  msOutPDWh$u, msOutDBh$u, msOutBLh$u, msOutBMWh$u,  msOutBTh$u, msOutSLh$u, msOutSDh$u)
colnames(allBLUPs) <- c("AshFDwPM","AshOnly","DWpM","WWP",  "PDW", "BDns", "BLen", "BMax",  "BThk", "SLen", "SDia")


### Combining all the GPs that have the GEBV
SProg_GPs<-allBLUPs[grep("UCONN-S",rownames(allBLUPs)),]  ### Grep only the SProg_GPs
allBLUPs0<-allBLUPs

### allBLUPs would remove all the SP_GP plots
allBLUPs<-allBLUPs[!rownames(allBLUPs)%in%rownames(SProg_GPs),]
###
  dim(allBLUPs)
  dim(allBLUPs0)
sampledSP <- which(nchar(rownames(allBLUPs)) < 11)  # The fndr SP is max of 10 characters
releaseGP <- nchar(rownames(allBLUPs))   
releaseGP <- which(10 < releaseGP & releaseGP <=15)  #FG and MG from fndr SP is max of 14 characters
progenySP <- which(nchar(rownames(allBLUPs)) > 15)   #the Cross name is larger than 15 characters

fndNames <- rownames(allBLUPs)[sampledSP]
hasDesc <- sapply(fndNames, function(n) grep(n, rownames(allBLUPs)[progenySP]))

nDesc <- sapply(hasDesc, length)
  nDesc

GP_allBLUPs<-as.data.frame(rbind(allBLUPs[releaseGP,],SProg_GPs)) ## !!! Added the SProg_GPs  
  dim(GP_allBLUPs)
#GP_allBLUPs<-as.data.frame(rbind(allBLUPs[releaseGP,]))
SP_allBLUPs<-as.data.frame(allBLUPs[progenySP,])
SP_allBLUPs<-SP_allBLUPs[order(-SP_allBLUPs$AshFDwPM),]

nrow(GP_allBLUPs)+length(fndNames)+nrow(SP_allBLUPs)
#517(=439 GPs_from_fndr+78 SProg_GPs) +104+ 245

pdf(paste0("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/","CorrPlot_SPplots",yr,"_PhotoScore23_WithHaploid.pdf"))
corrplot::corrplot.mixed(cor(SP_allBLUPs), diag="n", tl.cex=0.6, tl.col=1)
dev.off()


traitname<-"AshFDwPM" # !!!!!!!!!!
# Different colors for fndr, GP, progeny SP crosses
colSPGP <- c(rep("black", length(sampledSP)), rep("dark green", length(releaseGP)), rep("dark red", length(progenySP)),rep("red",nrow(SProg_GPs)))
colSPGP[which(nDesc == 0)] <- "grey"

OrderedBLUPs<-rbind(allBLUPs[sampledSP,],allBLUPs[releaseGP,],allBLUPs[progenySP,],SProg_GPs)
  dim(OrderedBLUPs)
  str(colSPGP)
#### Plot out each generation of BLUPs
pdf(paste0("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/",traitname,"vsIndividual",yr,"_PhotoScore23_withHaploid.pdf"))
plot(OrderedBLUPs[,colnames(OrderedBLUPs)%in%traitname], pch=16, col=colSPGP, xlab="Individual number", ylab=traitname) # !!The second col is DWpM
dev.off()


write.csv(OrderedBLUPs,paste0("allBLUPsDF_Ordred_",yr,"_Plots_WithHaploid.csv"))

# rnPos <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,colnames(allBLUPs)%in%traitname]>0))] #!! The second col is DWpM
# rnNeg <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,colnames(allBLUPs)%in%traitname]<0))] #!! The second col is DWpM
# locPos <- table(substring(rnPos, 6, 7))
# locNeg <- table(substring(rnNeg, 6, 7))

allBLUPsDF <- as.data.frame(allBLUPs)
allBLUPsDF <- cbind(pedigree=rownames(allBLUPs), allBLUPsDF)
head(allBLUPsDF)
tail(allBLUPsDF)

pdf(paste0("DWpMvsPDWsporophytes",yr,"_PhotoScore23.pdf"))
plot(allBLUPsDF$PDW[progenySP], allBLUPsDF$DWpM[progenySP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Sporophytes", cex.lab=1.3)
dev.off()
pdf(paste0("DWpMvsPDWgametophytes",yr,"_PhotoScore23.pdf"))
plot(allBLUPsDF$PDW[releaseGP], allBLUPsDF$DWpM[releaseGP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Gametophytes", cex.lab=1.3)
dev.off()

allBLUPsDF <- cbind(pedigree=allBLUPsDF$pedigree, index=2*allBLUPsDF$DWpM + allBLUPsDF$PDW, allBLUPsDF[,-1]) #Re

saveRDS(allBLUPsDF, file=paste0("allBLUPsDF_analyzeNH_",yr,"_PhotoScore23.rds"))
write.csv(allBLUPsDF,paste0("allBLUPs_DF_",yr,".csv"))





#######################################################
##{{}}
#########!!!!!!!!!!!!! Run this to get the outCovComb4 Only Once !!!!!!!!
source("is.square.matrix.R")
source("is.positive.definite.R")
source("is.symmetric.matrix.R")
load("GPsAmat_NA0.8_P1P2P3_09282020.Rdata")

#make positive definite: decimal points differ off-diagonal causes it being non-symmetric
# Further reduced the "SL18-LD-13-Mg-3" #GPsA<-GPsA[!rownames(GPsA)=="SL18-LD-13-Mg-3",!colnames(GPsA)=="SL18-LD-13-Mg-3"]
which(colnames(GPsA)=="SL18-PI-1-FG1")
IndiRM<-c("SL18-PI-1-FG1","3","SL18-CT1-FG3","SL18-CT1-MG2","SL18-CT1-MG3","SL18-OI-15-Female","SL18-ME-1-FG1","SL18-ME-1-MG1","SL18-ME-1-MG2")

GPsA<-GPsA[!rownames(GPsA)%in%IndiRM,!colnames(GPsA)%in%IndiRM] # 269 by 269
  dim(GPsA)

fndrsA<-round(mrkRelMat)     ### Added "shrink=TRUE" when estimating the mrkRelMat2
diag(fndrsA) <- diag(fndrsA) + 1e-5
is.positive.definite(fndrsA)   ### RelationshipMatrix_1. fnders A

GPsA<-round(GPsA,digits=5)   ### RelationshipMatrix_2.GPs A
diag(GPsA) <- diag(GPsA) + 1e-5
is.positive.definite(GPsA)  

#diag(mrkRelMat2) <- diag(mrkRelMat2) + 1e-5
diag(aMat) <- diag(aMat) + 1e-5
is.positive.definite(aMat)     ### RelationshipMatrix_3: pedigree based for ALL

#save(fndrsA,GPsA,aMat,file="CovList_3_As.Rdata")
#### Run these in terminal too !!!!!
#load("CovList_3_As.Rdata")
#install.packages("CovCombR")

library(CovCombR)

#### III. Make the combined RelationshipMatrix
### 1. NO initial, 3-list
CovList<-NULL
CovList[[1]]<-fndrsA ## fndrsA
CovList[[2]]<-GPsA  ## Further RMed the "SL18-LD-13-Mg-3"
CovList[[3]]<-aMat   

### outCovComb1<-CovComb(CovList,nu=1500)   

### 2. aMat initial, 2-list ------- Not correct, Only 239 invididuals
# CovList2<-NULL
# CovList2[[1]]<-fndrsA ## fndrsA
# CovList2[[2]]<-GPsA
# outCovComb2<-CovComb(CovList2,nu=1500,Kinit=aMat) 
####

### 3. aMat initial, 3-list
### outCovComb3<-CovComb(CovList,nu=1500,Kinit=aMat) 

### 4. aMat initial, 3-list, add weight
weights<-c(2,2,1)
outCovComb4<-CovComb(CovList,nu=1500,w=weights,Kinit=aMat) 

save(outCovComb4,file="outCovComb4_10012020_withSGP_866.Rdata")
#########!!!!!!!!!!!!! Run the above Only Once !!!!!!!!
##{{}}
########################################################




##### CV !!!!!!  Done only once 
#   #This is for setting CV, with a 20 folds scheme
#   #nreps<-round(nrow(dataNHpiYr)/20)+1
#   #sets<-rep(1:20,nreps)[-1:-(nreps*20-nrow(dataNHpiYr))]
#   #sets<-sets[order(runif(nrow(dataNHpiYr)))]
## Test pop: GPs used for making crosses (repeats + new combinations) and the ones genotyped
## Adding these union sets into a matrix


### random sampling to make 30 NAs on the dataNHpi
### test<-dataNHpi with NAs
### test rows of the u
### cor(test_u and the test_dataNHpi)

cycle=200
nSample=20
traitname<-"dryWgtPerM"

hMat <- calcHmatrix(mrkRelMat, aMat, aMatFounders=rownames(mrkRelMat))
hMat<-hMat[rownames(outCovComb4),colnames(outCovComb4)]
identical(colnames(hMat),colnames(outCovComb4))
aMat<-aMat[rownames(outCovComb4),colnames(outCovComb4)]

### cor=1 for outCovComb1 and outCovComb3 

CovCombList<-list(aMat,hMat,outCovComb1,outCovComb3,outCovComb4)

#corMeans<-matrix(nrow=length(CovCombList),ncol=2)


cor.r<-matrix(nrow=cycle,ncol=length(CovCombList))
colnames(cor.r)<-c("aMat","hMat","NoInit","aMatInit","aMatInit_withWeight")
rownames(cor.r)<-c(paste0("cor.r_",1:cycle))

cor.p<-matrix(nrow=cycle,ncol=length(CovCombList))
colnames(cor.p)<-c("aMat","hMat","NoInit","aMatInit","aMatInit_withWeight")
rownames(cor.p)<-c(paste0("cor.p_",1:cycle))

#corList<-NULL


for (i in 1:cycle){
  sampleYr<-sample(c(2019,2020),1,replace=TRUE)
  dataNHpiYr<-dataNHpi_RMchk[dataNHpi_RMchk$Year==sampleYr,]
  
  subsets<-order(runif(nrow(dataNHpiYr)))
  #print(subsets)
  pred.set<-dataNHpiYr[subsets<=nSample,]
  pred.set$Crosses<-as.factor(pred.set$Crosses)
  #print(str(pred.set[,1:15]))
  
  Y<-dataNHpi
  Y[rownames(Y)%in%rownames(pred.set),colnames(Y)==traitname]=NA
  
  for (j in 1:length(CovCombList)){
    hMat<-CovCombList[[j]]
    
    msOutDWPMh <- mixed.solve(y=log(Y[,traitname]+1), Z=msZ, K=hMat, X=msX, SE=T)
    AllhMatBLUPs<-msOutDWPMh$u
    pred.set.u<-AllhMatBLUPs[as.character(pred.set$Crosses)]
    pred.set.y<-pred.set[,colnames(pred.set)==traitname]
    
    cor.r[i,j]<-cor.test(pred.set.u,pred.set.y)$estimate
    cor.p[i,j]<-cor.test(pred.set.u,pred.set.y)$p.value
  }
}

#corList[[j]]<-cor
corMeans<-matrix(nrow=ncol(cor.r),ncol=2)
rownames(corMeans)<-colnames(cor.r)
colnames(corMeans)<-c("r","p")

corMeans[,1]<-t(colMeans(cor.r))
corMeans[,2]<-t(colMeans(cor.p))
corMeans

#r   p.value (20 cycles) (same sample set)
#aMat                0.2711160 0.3381301
#hMat                0.3094897 0.2670548
#NoInit              0.3030436 0.2853737
#aMatInit            0.3030436 0.2853737
#aMatInit_withWeight 0.3085474 0.2757711

#r   p.value (200 cycles) (same sample set)
#aMat                0.2295683 0.3912472
#hMat                0.2576408 0.3626239
#NoInit              0.2522064 0.3600591
#aMatInit            0.2522064 0.3600591
#aMatInit_withWeight 0.2544694 0.3559478

#r   p.value (200 cycles) (not the same sample set)
#NoInit              0.23319457 0.3716543
#aMatInit            0.24153700 0.3679468
#aMatInit_withWeight 0.23943776 0.3571202
##### CV ends here !!!!!




# Best sporophytes
nTop<-20  ### Pick the top 20
bestSP_DWpM <- allBLUPsDF[progenySP,][order(allBLUPsDF$DWpM[progenySP], decreasing=T)[1:nTop],] # based on DWpM
bestSP_DWpM[,-(1:2)] <- bestSP_DWpM[,-(1:2)] %>% round(3) # Only keep the numeric part
write.table(bestSP_DWpM, paste0("BestSPbyDryWgtPerM_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)

bestSP_idx <- allBLUPsDF[progenySP,][order(allBLUPsDF$index[progenySP], decreasing=T)[1:nTop],]
bestSP_idx[,-(1:2)] <- bestSP_idx[,-(1:2)] %>% round(3)
write.table(bestSP_idx, paste0("BestSPbyDWpMandPDW_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)

# Best gametophytes
bestGP_DWpM <- allBLUPsDF[releaseGP,][order(allBLUPsDF$DWpM[releaseGP], decreasing=T)[1:nTop],]
bestGP_DWpM[,-(1:2)] <- bestGP_DWpM[,-(1:2)] %>% round(3)
write.table(bestGP_DWpM[,-1], paste0("BestGPbyDryWgtPerM_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)
bestGP_idx <- allBLUPsDF[releaseGP,][order(allBLUPsDF$index[releaseGP], decreasing=T)[1:nTop],]
bestGP_idx[,-(1:2)] <- bestGP_idx[,-(1:2)] %>% round(3)
write.table(bestGP_idx[,-1], paste0("BestGPbyDWpMandPDW_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)

### save file to Diego
load("outCovComb4_10012020_withSGP.Rdata")
keep<-rownames(outCovComb4)[1:15]
dipkeep<-which(rownames(outCovComb4)%in%keep)
hapkeep<-grep(paste(c("OI","LD"), collapse="|"),rownames(outCovComb4))
outCovSample<-outCovComb4[c(dipkeep,hapkeep),c(dipkeep,hapkeep)]
other<-c("NL","CB","SF","CC","JS","NC","LL")
haprm<-grep(paste(other,collapse="|"),rownames(outCovSample))
outCovSample2<-outCovSample[-haprm,-haprm]

write.csv(outCovSample2,"outCovSample_TODiego.csv")

#### Datafile in the 2020_2019_Phenotypic_Data path
# #load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/CovComb/outCovComb4_and_Conden.Rdata")

