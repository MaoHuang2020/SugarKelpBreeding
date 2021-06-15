# Yr20 to Yr19
rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
load(paste0(WD,"OneTime1920/data/","outCovComb_dip_0116_2021.Rdata"))
#write.csv(outCovComb4_dipOrder,here("OneTime1920/data","A.csv"))

mm.file<- paste0(WD,"OneTime1920/data/","A.csv")          # path to covariates file

source(paste0(WD,"OneTime1920/code/","BGLR_functions_noFMLoc.R")) # terminal

####
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2  #Both Years

### Get the "BLUE_Trait"
load(paste0(WD,"OneTime1920/data/","BLUE_DwPM_2vs1Year.rdata"))

library(expss)
Y$BLUE_Trait<-vlookup(Y$Crosses,dict=CrossBLUE,result_column = "BLUE_DwPM_2Yrs",lookup_column = "CrossName") ### This is the DwPM
head(Y)

Inputfiledir<-c("OneTime1920/GP1/","OneTime1920/GP2/","OneTime1920/GP1P2/") 
head(Y)
dim(Y)
y<-Y[,"dryWgtPerM"]  # phenotypes column  !!!
gid<-Y[,"Crosses"]   # Crosses
yBLUE<-Y[,"BLUE_Trait"] # This is the BLUE for DwPM


### 4. Yr20to19

yrTP<-2019
yrPP<-2020
RmCommonX<-droplevels(Y[Y$Year==yrPP,]$Crosses[which(Y[Y$Year==yrPP,]$Crosses%in%Y[Y$Year==yrTP,]$Crosses)])  # the CommonCrosses 

# yBLUE<-Y[!Y$Crosses%in%RmCommonX,"BLUE_Trait"]
# gid<-Y[!Y$Crosses%in%RmCommonX,"Crosses"]
# y<-Y[!Y$Crosses%in%RmCommonX,"dryWgtPerM"]

setwd(paste0(WD,"OneTime1920/Yr20to19_output/"))

r<-NULL
for (i in 1:20){
  dir.create(paste0("Rep",i))
  savepath<-paste0("OneTime1920/Yr20to19_output/Rep",i,"/")  # the path following WD
  
  yNA<-y
  testing<-which(Y$Year==2019)  ## PP year
  yNA[testing]<-NA
  
  BetweenYear<-RunBGLR(YearEffects=FALSE,
                       nIter=80000,burnIn=60000,
                       y=yNA,testing=testing,
                       Inputfiledir=Inputfiledir,
                       Outputfiledir=savepath)
  r<-c(r,predict(testing=testing,yBLUE=yBLUE,Y=Y,fmfiledir=savepath))
}
 print(r)

r2<-NULL  # cor ONLY the ES plots,RM checks. OK to keep the plots that are the same across years
for (i in 1:20){
  savepath<-paste0("OneTime1920/Yr20to19_output/Rep",i,"/")
  Ypred<-read.csv(paste0(WD,savepath,"TP_predict_PP_noFMloc.csv"),sep=",",header=TRUE,row.names=1)
  yBLUE2<-  #Only Cor the ES plots
    yPred2<-Ypred[Ypred$popChk=="ES","yPred"]  #Only Cor the ES plots
  
  r2<-c(r2,cor(Ypred[Ypred$popChk=="ES","yBLUE"],yPred2,use="complete"))
}

write.csv(r, paste0(WD,savepath,"cor_",i,"Reps.csv"))
print(mean(r))  # -0.02784699

