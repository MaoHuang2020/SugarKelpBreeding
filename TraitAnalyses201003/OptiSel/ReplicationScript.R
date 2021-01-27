### File1 PedKelp file:
### All fndrs, Sire/Dam  0, sex NA
### GPs, Sire is its fndrs, Dam NA, sex according to itself
### ProgSP, Sire is parental MG, Dam is parental FG, sex NA
### ProgSP_GP, Sire being NA, Dam being ProgSP, sex NA

### File 2 phen:
### Merged, to get all the pedigree info (equiGen, etc)

### File 3 kinship:
### I used our outCovComb output estimated at diploid level

rm(list=ls())
install.packages("here")
library(here)

here::i_am("TraitAnalyses201003/OptiSel/ReplicationScript.R")
setwd("TraitAnalyses201003/OptiSel/Examples")


### Input the outCovComb_dip
load(here("TraitAnalyses201003/ReorderPedigree","outCovComb_dip_0116_2021.Rdata"))
  ls()
  outCovComb4_dipOrder[1:4,1:4]
  
### Input Pedi
pedir<-here("TraitAnalyses201003/ReorderPedigree","Ped_in_Order_866_Individuals_Fndr_New_Order_0116_2021.csv")
biphasicPedNH<-read.csv(file=pedir,sep=",",header=TRUE,row.names=1)

# ### Input raw Pheno
# load(here("TraitAnalyses201003","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))
# dataNHpi<-dataNHpiBoth_C  

### Input GEBVs
phenoBV<-read.csv(here("allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_0116_2021_dip.csv"),sep=",",header=T,row.names=1)

library("optiSel")
library("data.table")


### FILE 1: 
### Reformat pedigree matrix to match the package format ###
cols<-c("Indiv","Sire","Dam","Sex","Born")

PedKelp<-matrix(nrow=nrow(biphasicPedNH),ncol=length(cols))
rownames(PedKelp)<-rownames(biphasicPedNH)

colnames(PedKelp)<-cols
PedKelp<-as.data.frame(PedKelp)

PedKelp[,1:3]<-biphasicPedNH[,c(1,3,2)]  ### !!!! biphasicPedNH: fndRow, F(Dame), M(Sire)
  head(PedKelp)

### Distinguish each individual to fndr, GP, SP, SP_GP in order to define generation # and sexes of Indiv
find_Individual<-function(pedNameList){
  fndrRow<-which(!grepl("FG",pedNameList) & !grepl("MG",pedNameList))
  GPRow<-which(!grepl("x",pedNameList) & !grepl("UCONN-S",pedNameList) & c(grepl("FG",pedNameList) | grepl("MG",pedNameList)))
  SPRow<-which(grepl("x",pedNameList))
  SP_GPRow<-which(grepl("UCONN-S",pedNameList))
  return(list(fndrRow=fndrRow,GPRow=GPRow,SPRow=SPRow,SP_GPRow=SP_GPRow))
}

find_Individual<-find_Individual(rownames(PedKelp))

fndrRow<-find_Individual$fndrRow
GPRow<-find_Individual$GPRow
SPRow<-find_Individual$SPRow
SP_GPRow<-find_Individual$SP_GPRow

### Add the generation number to PedKelp$Born !!!!!!!!!!
PedKelp$Born[fndrRow]<-1 # fndrs
PedKelp$Born[GPRow]<-2
PedKelp$Born[SPRow]<-3
PedKelp$Born[SP_GPRow]<-4

### Finding the FG and MG sexes
find_FGMG<-function(GPNameList){
  FGP<-GPNameList[grep("FG",GPNameList)]  # FGP list
  MGP<-GPNameList[grep("MG",GPNameList)]  # MGP list
  return(list(FGP=FGP,MGP=MGP))
}

GPNameList<-rownames(PedKelp)[c(GPRow,SP_GPRow)]  ### !!!!

FGP<-find_FGMG(GPNameList)$FGP
MGP<-find_FGMG(GPNameList)$MGP

for (i in 1:nrow(PedKelp)){
  if (rownames(PedKelp)[i]%in%FGP==TRUE){
    PedKelp$Sex[i]<-2                           # female
  }else if(rownames(PedKelp)[i]%in%MGP==TRUE){
    PedKelp$Sex[i]<-1                           # male
  }else{
    is.na<-PedKelp$Sex[i]
  }
}

### Prepare/Pruning the Pedigree file "PedKelp"
PedKelp$Indiv2<-rownames(PedKelp)
Pedig<-prePed(PedKelp)
rownames(Pedig)<-Pedig$Indiv2

PedKelp<-PedKelp[,!colnames(PedKelp)=="Indiv2"]
Pedig<-Pedig[,!colnames(Pedig)=="Indiv2"]
  
### The sexes turned out to be wrong, after checking their source code!!!
library(expss)
Pedig$Sex<-vlookup(rownames(Pedig),dict=PedKelp,result_column = "Sex",lookup_column = "row.names")


### FILE 2:  
### Merge phenotypic file to the pedigree info
phen<-merge(phenoBV,PedKelp,all.x=TRUE,by.x="row.names",by.y="row.names")
  head(phen)
phen<-phen[order(match(phen$Row.names,rownames(phenoBV))),]

### Compute the number of equivalent complete generations in the pedigree of each individual ###
### equiGen: Number of equivalent complete generations, defined as the sum of the proportion of known ancestors over all generations traced,
Sy    <- summary(Pedig) 
  head(Sy)
phen  <- merge(phen, Sy[, c("Indiv", "equiGen")], on="Indiv")

phen$Sex[phen$Sex==1]<-"male"
phen$Sex[phen$Sex==2]<-"female"
phen$isCandidate <- phen$Born==2 |phen$Born==4 


##### Checking on completeness
compl <- completeness(Pedig, keep=phen$Indiv, by="Indiv")
  head(compl)
  tail(compl)
compl <- completeness(Pedig, keep=phen$Indiv, by="Sex")
library("ggplot2")
ggplot(compl, aes(x=Generation, y=Completeness, col=Sex)) + 
  geom_line()
#######


### L is the generation interval in years
### approximately set to be 2
 L<-2
 
#### FILE 3:
fPED<-outCovComb4_dipOrder
  #fPED<-pedIBD(PedKelp)      # This one is identical to using Pedig
  #fPED2  <- pedIBD(Pedig)    # only r=0.737 with the cor(c(outCovComb4_dipOrder),c(fPED))

save(PedKelp,Pedig,phen,fPED,file="Prepared_PedKelp_Pedig_phen.Rdata") 

### Convert the Indiv from numeric back to character, so that it matches those in fPED
phen$IndivNum<-as.factor(phen$Indiv)    
phen$Indiv<-phen$Row.names    # has to be "Indiv", matching fPED names
  head(phen)
  str(phen)
  
cand <- candes(phen=phen, fPED=fPED, cont=NULL)
  cand$mean

### Traditional OCS with upper bounds for not just female, no care sex, contributions ###
n<-5 # number of progeny per individual
N0<-400 # evaluation population size
threshold<-n/(2*N0)  

### Define upper limits for kinships and native kinships ###
Ne<-60
ub.fPED  <- cand$mean$fPED  + (1-cand$mean$fPED)/(2*Ne*L)

females <- cand$phen$Sex=="female" & cand$phen$isCandidate
ub  <- setNames(rep(threshold, sum(females)), cand$phen$Indiv[females])
con <- list(ub         = ub, 
            ub.fPED    = ub.fPED,  
            lb.equiGen = cand$mean$equiGen)

### Calculate the optimum contribution
fit<-opticont("max.AshFDwPM",cand,con,solver="slsqp",trace=FALSE)  ### choose a trait!  "slsqp"
fit$info
fit$mean
Candidate <- fit$parent
  head(Candidate[Candidate$oc>0.01,c("Sex","AshFDwPM","oc")])   ### choose a trait!
Candidate_rank<-Candidate[order(-Candidate$AshFDwPM,-Candidate$oc),]
  head(Candidate)

####!!!!!!!!!!! Question, why did this not work???
### Estimate population means ###
cont <- agecont(Pedig, use=!is.na(Pedig$Born),maxAge = NA)  ### Did not work ?????
head(cont) 

L  <- 1/(4*cont$male[1]) + 1/(4*cont$female[1])  ### Cannot get L due to NAs in the cont step

cont<-genecont(Pedig,from=NULL,to=NULL)  ### Did not work??
str(cont)







# Ignore these below:
# ###########OLDER way of estimating the PedKelp
# 
# rownames(PedKelp)<-biphasicPedNH$fndRow
# PedKelp<-as.data.frame(PedKelp)
# PedKelp$Indiv<-c(rownames(biphasicPedNH))
# 
# Name_Order<-c(rownames(biphasicPedNH))
# names(Name_Order)<-c(biphasicPedNH[,1])
# head(Name_Order)
# ### Find the "sire" "dam" name according to the biphasicPedNH for PedKelp
# PedName<-function(row){
#   if(row==0 | is.na(row)){
#     value<-row
#   }else{
#     value<-Name_Order[row]
#   }
#   return(value)
# }
# 
# for (i in 1:nrow(PedKelp)){
#   PedKelp$Sire[i]<-PedName(biphasicPedNH[i,2])
#   PedKelp$Dam[i]<-PedName(biphasicPedNH[i,3])
# }
# 
# # if without FG, it is fndr, 2018 (generation 1)
# # if with FG, MG but no "x" or UCONN, it is GP, 2018 (generation 2)
# # if with "x", it is crosses made in 2018 and 2019 -> need to match to which year (generation 3)
# # if with UCONN, it is UCONN_GP, 2019 (generation 4)
# 
# find_Individual<-function(pedNameList){
#   fndrRow<-which(!grepl("FG",pedNameList) & !grepl("MG",pedNameList))
#   GPRow<-which(!grepl("x",pedNameList) & !grepl("UCONN-S",pedNameList) & c(grepl("FG",pedNameList) | grepl("MG",pedNameList)))
#   SPRow<-which(grepl("x",pedNameList))
#   SP_GPRow<-which(grepl("UCONN-S",pedNameList))
#   return(list(fndrRow=fndrRow,GPRow=GPRow,SPRow=SPRow,SP_GPRow=SP_GPRow))
# }
# 
# ### Add the generation number to PedKelp$Born
# find_Individual<-find_Individual(PedKelp$Indiv)
# fndrRow<-find_Individual$fndrRow
# GPRow<-find_Individual$GPRow
# SPRow<-find_Individual$SPRow
# SP_GPRow<-find_Individual$SP_GPRow
# 
# 
# 
# PedKelp$Born[fndrRow]<-1 # fndrs
# PedKelp$Born[GPRow]<-2
# PedKelp$Born[SPRow]<-3
# PedKelp$Born[SP_GPRow]<-4
# 
# is.na<-PedKelp$Sex
# str(PedKelp)
# 
# ### Finding the FG and MG rows
# find_FGMG<-function(GPNameList){
#   FGP<-GPNameList[grep("FG",GPNameList)]  # FGP list
#   MGP<-GPNameList[grep("MG",GPNameList)]  # MGP list
#   return(list(FGP=FGP,MGP=MGP))
# }
# 
# FGP<-find_FGMG(PedKelp$Indiv[c(GPRow,SP_GPRow)])$FGP
# MGP<-find_FGMG(PedKelp$Indiv[c(GPRow,SP_GPRow)])$MGP
# 
# 
# for (i in 1:nrow(PedKelp)){
#   if (rownames(PedKelp)[i]%in%FGP==TRUE){
#     PedKelp$Sex[i]<-2
#   }else if(rownames(PedKelp)[i]%in%MGP==TRUE){
#     PedKelp$Sex[i]<-1
#   }else{
#     is.na<-PedKelp$Sex[i]
#   }
# }
# 
# #### Sorry, Decided to convert the PedKelp back to numeric form
# rownames(PedKelp)<-rownames(biphasicPedNH)
# PedKelp$Indiv<-biphasicPedNH$fndRow
# PedKelp$Sire<-biphasicPedNH$X.1
# PedKelp$Dam<-biphasicPedNH$X.2
# PedKelp[103:105,]
# 
# head(PedKelp)  
# PedKelp[103:105,]





#### Other codes might be useful

# ### Estimate native contributions from pedigrees ###
# Pedig2 <- prePed(Ped, lastNative=1970, thisBreed="Angler", keep=animals)
# BC     <- pedBreedComp(Pedig2, thisBreed="Angler")
# phen   <- merge(phen, BC[, c("Indiv", "unknown", "native")], on="Indiv")
# setnames(phen, old="native", new="pedNC")
# BC[I, 1:4, on="Indiv"]

# ### Estimate native kinships from pedigrees ###
# fPEDN  <- pedIBDatN(Pedig2, thisBreed="Angler", keep.only=phen$Indiv)
# natKin <- fPEDN$Q1/fPEDN$Q2
# natKin[I, I]
# 

#con<-list(uniform="female",ub.sKin=0.010)
# This is assuming female contributions uniform. derive the ub.sKin/ub.fPED. But it needs L !!!! Did not quite work
# Ne <- 60
# con <- list(
#   uniform = "female",
#   ub.fPED = 1-(1-cand$mean$fPED)*(1-1/(2*Ne))^(1/L),
#   lb.equiGen = cand$mean$equiGen
# )

# ### Mate allocation
# 
# con <- list(ub=ub, ub.fSEG=ub.fSEG)
# fit <- opticont("max.EBV", cand, con, solver="cccp")
# 
# Candidate <- fit$parent
# Candidate$n <- noffspring(Candidate, N=200, random=FALSE)$nOff
# Mating      <- matings(Candidate, fSEG)
# head(Mating)
# 
# attributes(Mating)$objval

