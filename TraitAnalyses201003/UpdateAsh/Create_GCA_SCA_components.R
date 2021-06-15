#### making the GCA and SCA variance components 

# In terminal

rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot

Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])


### Run it only once
### This is to only get the calc GCA and SCA functions
source(paste0(WD,"OneTime1920/code/","BGLR_functions_noFMLoc.R")) # !!!terminal

### Write and Load the A matrix
load(paste0(WD,"OneTime1920/data/","outCovComb4_Mix_Conden_0420_2021.Rdata"))
write.csv(outCovComb4_MixOrder,paste0(WD,"OneTime1920/data/","A.csv"))
mm.file<- paste0(WD,"OneTime1920/data/","A.csv")          # path to covariates file

# # #the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
 CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = 3,colNam = "P1",savefiledir = "OneTime1920/GP1/250Individual/")
 CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = 5,colNam = "P2",savefiledir = "OneTime1920/GP2/250Individual/")
# #
  CalSCA(G1.file=paste0(WD,"OneTime1920/GP1/250Individual/","G.rda"),
         G2.file=paste0(WD,"OneTime1920/GP2/250Individual/","G.rda"),
         savefileDir="OneTime1920/GP1P2/250Individual/")
# # ############# Run it only once
