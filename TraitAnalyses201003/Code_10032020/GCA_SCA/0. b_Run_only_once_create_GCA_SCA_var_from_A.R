# Run Only Once
# Generate GCA and SCA within Both Years'
# Generate GCA and SCA var within 2019, then within 2020

####### Both Years

rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
#load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi
####
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])

########### Run only Once !!!!!!!!!To get the GCA and SCA components
load(paste0("/local/workdir/mh865/outCovComb/outCovComb4_Mix_Conden_0527_2021.Rdata")) ##!!!
A<-outCovComb4_MixOrder
  IDs<-unique(droplevels(dataNHpiBoth_C$Crosses))

write.csv(A,paste0(WD,"OneTime1920/data/","A.csv"))
mm.file<- paste0(WD,"OneTime1920/data/","A.csv")          # path to covariates file

## This is to only get the calc GCA and SCA functions
source(paste0(WD,"OneTime1920/code/","BGLR_functions_noFMLoc.R")) # !!!terminal

#the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = which(colnames(Y)=="femaPar"),colNam = "P1",savefiledir = "OneTime1920/GP1/250Individual/")
CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = which(colnames(Y)=="malePar"),colNam = "P2",savefiledir = "OneTime1920/GP2/250Individual/")

CalSCA(G1.file=paste0(WD,"OneTime1920/GP1/250Individual/","G.rda"),
       G2.file=paste0(WD,"OneTime1920/GP2/250Individual/","G.rda"),
       savefileDir="OneTime1920/GP1P2/250Individual/")


############# Run it only once. To get the sampling CV
# # #### Run it once
#  sampleCV<-matrix(nrow=nrow(Y),ncol=500)
#   for (n in 1:500){
#     sets<-rep(1:10,25)  #250 ES ones
#     sampleCV[,n]<-sets[order(runif(nrow(Y)))]
#   }
#  save(sampleCV,file="sampleCV_250Indiv_0527_2021.Rdata")
 
# #  ###### One Time


##### 2019

rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
####
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])

library(BGLR)

########### Run only Once !!!!!!!!!To get the GCA and SCA components

# load(paste0(WD,"OneTime1920/data/","outCovComb4_Mix_Conden_0527_2021.Rdata"))
# Already run from above Both Year's section
# A<-outCovComb4_MixOrder
# write.csv(A,here("OneTime1920/data","A.csv"))
mm.file<- paste0(WD,"OneTime1920/data/","A.csv")          # path to covariates file

## This is to only get the calc GCA and SCA functions
source(paste0(WD,"OneTime1920/code/","BGLR_functions_noFMLoc.R")) # !!!terminal


############# Rerun the Y statement, then Run this separately for each Year
Y<-droplevels(Y[Y$Year==2019,])   ##!!!
Yr<-2019                          ##!!!
folder<-"OneTime1920/Yr19Only/"   ##!!!
nreps<-13
subtract<-c(1:2)
#############


############# Rerun the Y statement, then Run this separately for each Year
Y<-droplevels(Y[Y$Year==2020,])    ##!!!
Yr<-2020                           ##!!!
folder<-"OneTime1920/Yr20Only/"    ##!!!
nreps<-13
subtract<-c(1:2)
################

#the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = which(colnames(Y)=="femaPar"),colNam = "P1",savefiledir = paste0(folder,"GP1/"))
CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = which(colnames(Y)=="malePar"),colNam = "P2",savefiledir = paste0(folder,"GP2/"))

CalSCA(G1.file=paste0(WD,folder,"GP1/","G.rda"),
       G2.file=paste0(WD,folder,"GP2/","G.rda"),
       savefileDir=paste0(folder,"GP1P2/"))
############# Run it only once


#### Run it once
sampleCV<-matrix(nrow=nrow(Y),ncol=500)
for (n in 1:500){
  sets<-rep(1:10,nreps)[-subtract]  #122 ES ones
  sampleCV[,n]<-sets[order(runif(nrow(Y)))]
}
save(sampleCV,file=paste0(WD,"OneTime1920/data/sampleCV_Yr",Yr,"_",nrow(Y),"Indiv_0423_2021.Rdata"))
###### One Time
############## Run only Once !!!!!!!!!
