rm(list=ls())
install.packages("here")
library(here)

here::i_am("TraitAnalyses201003/OptiSel/ReplicationScript.R")
#setwd("TraitAnalyses201003/OptiSel/Examples")


### Input the outCovComb_dip
load(here("TraitAnalyses201003/ReorderPedigree","outCovComb_dip_0116_2021.Rdata"))
  ls()
  outCovComb4_dipOrder[1:4,1:4]
  
### Pedi
pedir<-here("TraitAnalyses201003/ReorderPedigree","Ped_in_Order_866_Individuals_Fndr_New_Order_0116_2021.csv")
biphasicPedNH<-read.csv(file=pedir,sep=",",header=TRUE,row.names=1)

### raw Pheno
load(here("TraitAnalyses201003","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))
dataNHpi<-dataNHpiBoth_C  

### GEBVs
phenoBV<-read.csv(here("allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_0116_2021_dip.csv"),sep=",",header=T,row.names=1)

library("optiSel")
library("data.table")

#animals <- read.indiv("Population/Angler.Chr1.phased")
#I <- c("animal7396", "animal8713", "animal11514" )

### Estimate kinships from pedigrees ###
cols<-c("Indiv","Sire","Dam","Sex","Born")

PedKelp<-matrix(nrow=nrow(biphasicPedNH),ncol=length(cols))
colnames(PedKelp)<-cols

### All fndrs, Sire/Dam being NA, sex NA
### GPs, Sire being NA, Dam being fndrs, sex according to itself
### ProgSP, Sire being MG, Dam being MG, sex NA
### ProgSP_GP, Sire being NA, Dam being ProgSP, sex NA
rownames(PedKelp)<-biphasicPedNH$fndRow
PedKelp<-as.data.frame(PedKelp)
PedKelp$Indiv<-c(rownames(biphasicPedNH))

Name_Order<-c(rownames(biphasicPedNH))
names(Name_Order)<-c(biphasicPedNH[,1])
  head(Name_Order)
  
#PedKelp$Sire<-Name_Order[biphasicPedNH[,1]]

### Find the "sire" "dam" name according to the biphasicPedNH for PedKelp
PedName<-function(row){
  if(row==0 | is.na(row)){
    value<-row
  }else{
  value<-Name_Order[row]
  }
  return(value)
}

for (i in 1:nrow(PedKelp)){
  PedKelp$Sire[i]<-PedName(biphasicPedNH[i,2])
  PedKelp$Dam[i]<-PedName(biphasicPedNH[i,3])
}

# if without FG, it is fndr, 2018 (generation 1)
# if with FG, MG but no "x" or UCONN, it is GP, 2018 (generation 2)
# if with "x", it is crosses made in 2018 and 2019 -> need to match to which year (generation 3)
# if with UCONN, it is UCONN_GP, 2019 (generation 4)

fndrRow<-which(!grepl("FG",PedKelp$Indiv) & !grepl("MG",PedKelp$Indiv))
GPRow<-which(!grepl("x",PedKelp$Indiv) & !grepl("UCONN-S",PedKelp$Indiv) & c(grepl("FG",PedKelp$Indiv) | grepl("MG",PedKelp$Indiv)))
SPRow<-which(grepl("x",PedKelp$Indiv))
SP_GPRow<-which(grepl("UCONN-S",PedKelp$Indiv))

### Add the generation number to PedKelp$Born
PedKelp$Born[fndrRow]<-1 # fndrs
PedKelp$Born[GPRow]<-1
PedKelp$Born[SPRow]<-2
PedKelp$Born[SP_GPRow]<-2

is.na<-PedKelp$Sex
  str(PedKelp)

FGP<-PedKelp$Indiv[c(GPRow,SP_GPRow)][grep("FG",PedKelp$Indiv[c(GPRow,SP_GPRow)])]  # FGP list
MGP<-PedKelp$Indiv[c(GPRow,SP_GPRow)][grep("MG",PedKelp$Indiv[c(GPRow,SP_GPRow)])]  # MGP list

for (i in 1:nrow(PedKelp)){
  if (PedKelp$Indiv[i]%in%FGP==TRUE){
    PedKelp$Sex[i]<-"female"
  }else if(PedKelp$Indiv[i]%in%MGP==TRUE){
    PedKelp$Sex[i]<-"male"
  }else{
    is.na<-PedKelp$Sex[i]
  }
}

#Ped   <- fread("Population/Pedigree.txt")
Pedig<-prePed(PedKelp)
### The sexes turned out to be wrong !!!
library(expss)
Pedig$Sex<-vlookup(Pedig$Indiv,dict=PedKelp,result_column = "Sex",lookup_column = "Indiv")
  head(Pedig)
  tail(Pedig)
fPED  <- pedIBD(Pedig)    # only r=0.737 with the cor(c(outCovComb4_dipOrder),c(fPED))


phen<-merge(phenoBV,PedKelp,all.x=TRUE,by.x="row.names",by.y="Indiv")
colnames(phen)[colnames(phen)=="Row.names"]<-"Indiv"
  head(phen)
phen<-phen[order(match(phen$Indiv,rownames(phenoBV))),]


compl <- completeness(Pedig, keep=phen$Indiv, by="Indiv")
head(compl)

### Compute the number of equivalent complete     ###
### generations in the pedigree of each individual ###
# 
#   phen2  <- fread("Population/BreedingValues.txt")
#   phen2<-phen[1:6,]
#   

Sy    <- summary(Pedig)  #### How to get the equiGen, fullGen, maxGen, PCI in our case?
phen  <- merge(phen, Sy[, c("Indiv", "equiGen")], on="Indiv")

# ### Estimate kinships from phased marker data  ###
# 
# bfiles <- paste0("Population/Angler.Chr", 1:29, ".phased")
# map   <- fread("Population/map.txt")
# fSEG  <- segIBD(bfiles, map, minSNP=20, minL=2.5, keep=animals)
# fSEG[I,I]
# 
# ### Estimate native contributions from pedigrees ###
# 
# Pedig2 <- prePed(Ped, lastNative=1970, thisBreed="Angler", keep=animals)
# BC     <- pedBreedComp(Pedig2, thisBreed="Angler")
# phen   <- merge(phen, BC[, c("Indiv", "unknown", "native")], on="Indiv")
# setnames(phen, old="native", new="pedNC")
# BC[I, 1:4, on="Indiv"]
# 
# 
# ### Estimate native contributions from phased marker data ###
# 
# bfiles <- paste0("Population/Angler.Chr",     1:29, ".phased")
# rfiles <- paste0("refBreeds/OtherBreeds.Chr", 1:29, ".phased")
# files  <- list(hap.thisBreed = bfiles, hap.refBreeds=rfiles)
# Cattle <- fread("genotypedIndiv.txt")
# wfile  <- haplofreq(files, Cattle, map, thisBreed="Angler",
#                        minSNP=20,  minL=2.5, ubFreq=0.01,  
#                        what="match", w.dir="Population")
# Comp   <- segBreedComp(wfile$match, map)
# phen   <- merge(phen, Comp[, c("Indiv","native")], on="Indiv")
# setnames(phen, old="native", new="segNC")     
# 
# ### Estimate native kinships from pedigrees ###
# 
# fPEDN  <- pedIBDatN(Pedig2, thisBreed="Angler", keep.only=phen$Indiv)
# natKin <- fPEDN$Q1/fPEDN$Q2
# natKin[I, I]
# 
# ### Estimate native kinships from phased marker data ###
# 
# 
# fSEGN <- segIBDatN(files, Cattle, map, thisBreed="Angler", 
#                       minSNP=20, ubFreq=0.01, minL=2.5)
# natKin <- fSEGN$Q1/fSEGN$Q2
# natKin[I, I]

### Estimate population means ###

 cont <- agecont(Pedig, use=!is.na(Pedig$Born),maxAge = NA)  ### Did not work ?????
 head(cont) 
 
 Ne <- 60
 L  <- 1/(4*cont$male[1]) + 1/(4*cont$female[1])  ### Cannot get L due to NAs in the cont step
 
 cont<-genecont(Pedig,from=NULL,to=NULL)  ### Did not work??
 head(cont)
# 
# phen$isCandidate <- phen$Born <= 2026
 
 # cont<-pedBreedComp(Pedig)

cand <- candes(phen=phen, fPED=fPED, cont=NULL)
cand$mean

#con<-list(uniform="female",ub.sKin=0.047)
# derive the ub.sKin. But it needs L !!!!
Ne <- 60
con <- list(
  uniform = "female",
  ub.sKin = 1-(1-cand$mean$fPED)*(1-1/(2*Ne))^(1/L)
)


Offspring<-opticont("max.BV",cand,con,trace=FALSE)


### Define upper limits for kinships and native kinships ###


# ub.fSEG  <- cand$mean$fSEG  + (1-cand$mean$fSEG)/(2*Ne*L)
ub.fPED  <- cand$mean$fPED  + (1-cand$mean$fPED)/(2*Ne*L)

# ub.fSEGN <- cand$mean$fSEGN + (1-cand$mean$fSEGN)/(2*Ne*L)
# ub.fPEDN <- cand$mean$fPEDN + (1-cand$mean$fPEDN)/(2*Ne*L)

# ### Traditional OCS of male contributions ###
# 
# con <- list(uniform    = "female", 
#             ub.fPED    = ub.fPED, 
#             lb.equiGen = cand$mean$equiGen)
# 
# fit <- opticont("max.EBV", cand, con, solver="cccp")
# 
# fit$info
# 
# fit$mean
# 
# Candidate <- fit$parent
# Candidate[Candidate$oc>0.01,c("Sex","EBV","oc")]

### Traditional OCS with upper bounds for female contributions ###

females <- cand$phen$Sex=="female" & cand$phen$isCandidate
ub  <- setNames(rep(0.0125, sum(females)), cand$phen$Indiv[females])
con <- list(ub         = ub, 
            ub.fPED    = ub.fPED,  
            lb.equiGen = cand$mean$equiGen)

# ### Advanced OCS ###
# 
# con <- list(ub=ub, ub.fSEG=ub.fSEG)
# fit <- opticont("max.EBV", cand, con, solver="cccp")
# 
# con <- list(ub=ub)
# fit <- opticont("min.fSEG", cand, con, solver="cccp")
# 
# con <- list(ub=ub, lb.EBV=101)
# fit <- opticont("min.fSEG", cand, con, solver="cccp")
# 
# #con <- list(ub=ub, lb.EBV=101, lb.segNC=0.8)
# #fit <- opticont("min.fSEG", cand, con, solver="cccp")
# 
# con <- list(ub=ub, ub.fSEGN=ub.fSEGN)
# fit <- opticont("max.segNC", cand, con, solver="cccp")

##### Mate allocation

con <- list(ub=ub, ub.fSEG=ub.fSEG)
fit <- opticont("max.EBV", cand, con, solver="cccp")

Candidate <- fit$parent
Candidate$n <- noffspring(Candidate, N=200, random=FALSE)$nOff
Mating      <- matings(Candidate, fSEG)
head(Mating)

attributes(Mating)$objval

