## ScenarioI_Yr19_Yr21 OneTime Running
rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/OneTime1921/"  # run in terminal
datafdr<-paste0(WD,"data/")

load(paste0(datafdr,"dataNHpi_withChk_3Yrs_PhotoScore23_07122021.rdata"))   ## Plot
#load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi
####
dataNHpi<-dataNHpi3yrs_C
dataNHpi$densityBlades<-ifelse(dataNHpi$densityBlades==0,NA,dataNHpi$densityBlades)  # densityblades as 0 then NA
dataNHpi$popChk <- ifelse(substr(dataNHpi$plotNo, 1, 1) == "Z", substr(dataNHpi$plotNo, 1, 2), "ES")  # Checks VS ES

dataNHpi$withinLoc <- ifelse(as.vector(dataNHpi$femaParLoc) == as.vector(dataNHpi$maleParLoc), 1, 0) # WithinLoc is 1
dataNHpi$development[dataNHpi$development=="#N/A"] <-NA


Y1<-dataNHpi  #dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM") #,"AshFreedryWgtPerM"
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])

########### Run only Once !!!!!!!!!To get the GCA and SCA components
load(paste0("/local/workdir/mh865/outCovComb/outCovComb4_Mix_Conden_0527_2021.Rdata")) ##!!!
A<-Pediconden
write.csv(A,paste0(datafdr,"A.csv"))
mm.file<- paste0(datafdr,"A.csv")          # path to covariates file

## This is to only get the calc GCA and SCA functions
source(paste0("/local/workdir/mh865/GCA_SCA/OneTime1920/code/","BGLR_functions_noFMLoc.R")) # !!!terminal

#the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = which(colnames(Y)=="femaPar"),colNam = "P1",savefiledir = "OneTime1920/GP1/")
CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = which(colnames(Y)=="malePar"),colNam = "P2",savefiledir = "OneTime1920/GP2/")

CalSCA(G1.file=paste0(WD,"OneTime1920/GP1/","G.rda"),
       G2.file=paste0(WD,"OneTime1920/GP2/","G.rda"),
       savefileDir="OneTime1920/GP1P2/")


############# Run it only once. To get the sampling CV
# # #### Run it once
# sampleCV<-matrix(nrow=nrow(Y),ncol=500)
#  for (n in 1:500){
#    sets<-rep(1:10,25)  #250 ES ones
#    sampleCV[,n]<-sets[order(runif(nrow(Y)))]
#  }
# save(sampleCV,file="sampleCV_250Indiv_0527_2021.Rdata")

# #  ###### One Time

### NOTE!!! Pediconden and outCovComb4_MixOrder is loosely correlated
# cultevo::mantel.test(dist(outCovComb4_MixOrder),dist(Pediconden),trials=99) 
# Mantel permutation test (method: spearman)
# r = 0.639, N = 374545
# 99 permutations, mean = -0.00607, sd = 0.0276
# p (empirical) = 0.01 (veridical correlation is highest found)


##### 2019
rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/PedigreeOnly_GCA_SCA/"  # run in terminal
datafdr<-paste0(WD,"data/")

load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
####
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])

library(BGLR)

########### Run only Once !!!!!!!!!To get the GCA and SCA components

# load(paste0(datafdr,"outCovComb4_Mix_Conden_0527_2021.Rdata"))
# Already run from above Both Year's section
# A<-Pediconden
# write.csv(A,here("OneTime1920/data","A.csv"))
mm.file<- paste0(datafdr,"A.csv")          # path to covariates file

## This is to only get the calc GCA and SCA functions
source(paste0("/local/workdir/mh865/GCA_SCA/OneTime1920/code/","BGLR_functions_noFMLoc.R")) # !!!terminal


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


# #### Run it once
# sampleCV<-matrix(nrow=nrow(Y),ncol=500)
# for (n in 1:500){
#   sets<-rep(1:10,nreps)[-subtract]  #122 ES ones
#   sampleCV[,n]<-sets[order(runif(nrow(Y)))]
# }
# save(sampleCV,file=paste0(datafdr,"sampleCV_Yr",Yr,"_",nrow(Y),"Indiv_0423_2021.Rdata"))
# ###### One Time
# ############## Run only Once !!!!!!!!!
# 















rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/PedigreeOnly_GCA_SCA/"  # run in terminal
datafdr<-paste0(WD,"data/")

load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
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
A<-Pediconden
write.csv(A,paste0(datafdr,"A.csv"))
mm.file<- paste0(datafdr,"A.csv")          # path to covariates file

## This is to only get the calc GCA and SCA functions
source(paste0("/local/workdir/mh865/GCA_SCA/OneTime1920/code/","BGLR_functions_noFMLoc.R")) # !!!terminal

#the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = which(colnames(Y)=="femaPar"),colNam = "P1",savefiledir = "OneTime1920/GP1/")
CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = which(colnames(Y)=="malePar"),colNam = "P2",savefiledir = "OneTime1920/GP2/")

CalSCA(G1.file=paste0(WD,"OneTime1920/GP1/","G.rda"),
       G2.file=paste0(WD,"OneTime1920/GP2/","G.rda"),
       savefileDir="OneTime1920/GP1P2/")


############# Run it only once. To get the sampling CV
# # #### Run it once
# sampleCV<-matrix(nrow=nrow(Y),ncol=500)
#  for (n in 1:500){
#    sets<-rep(1:10,25)  #250 ES ones
#    sampleCV[,n]<-sets[order(runif(nrow(Y)))]
#  }
# save(sampleCV,file="sampleCV_250Indiv_0527_2021.Rdata")

# #  ###### One Time

### NOTE!!! Pediconden and outCovComb4_MixOrder is loosely correlated
# cultevo::mantel.test(dist(outCovComb4_MixOrder),dist(Pediconden),trials=99) 
# Mantel permutation test (method: spearman)
# r = 0.639, N = 374545
# 99 permutations, mean = -0.00607, sd = 0.0276
# p (empirical) = 0.01 (veridical correlation is highest found)


##### 2019
rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/PedigreeOnly_GCA_SCA/"  # run in terminal
datafdr<-paste0(WD,"data/")

load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
####
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])

library(BGLR)

########### Run only Once !!!!!!!!!To get the GCA and SCA components

# load(paste0(datafdr,"outCovComb4_Mix_Conden_0527_2021.Rdata"))
# Already run from above Both Year's section
# A<-Pediconden
# write.csv(A,here("OneTime1920/data","A.csv"))
mm.file<- paste0(datafdr,"A.csv")          # path to covariates file

## This is to only get the calc GCA and SCA functions
source(paste0("/local/workdir/mh865/GCA_SCA/OneTime1920/code/","BGLR_functions_noFMLoc.R")) # !!!terminal


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


# #### Run it once
# sampleCV<-matrix(nrow=nrow(Y),ncol=500)
# for (n in 1:500){
#   sets<-rep(1:10,nreps)[-subtract]  #122 ES ones
#   sampleCV[,n]<-sets[order(runif(nrow(Y)))]
# }
# save(sampleCV,file=paste0(datafdr,"sampleCV_Yr",Yr,"_",nrow(Y),"Indiv_0423_2021.Rdata"))
# ###### One Time
# ############## Run only Once !!!!!!!!!
# 
