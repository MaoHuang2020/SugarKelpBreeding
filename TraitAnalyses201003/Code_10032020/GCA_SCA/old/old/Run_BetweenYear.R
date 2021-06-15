#Run_betweenYear

rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
load(paste0(WD,"OneTime1920/data/","outCovComb_dip_0116_2021.Rdata"))
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

### 3. Yr19to20
yrTP<-2019
yrPP<-2020
RmCommonX<-droplevels(Y[Y$Year==yrPP,]$Crosses[which(Y[Y$Year==yrPP,]$Crosses%in%Y[Y$Year==yrTP,]$Crosses)])  # the CommonCrosses 

setwd(paste0(WD,"OneTime1920/Yr19to20_output/"))

yBLUE<-Y[!Y$Crosses%in%RmCommonX,"BLUE_Trait"]
gid<-Y[!Y$Crosses%in%RmCommonX,"Crosses"]
y<-Y[!Y$Crosses%in%RmCommonX,"dryWgtPerM"]

r<-NULL 
for (i in 1:5){
  dir.create(paste0("Rep",i))
  savepath<-paste0("OneTime1920/Yr19to20_output/Rep",i,"/")  # the path following WD
  
  yNA<-y
  testing<-which(Y$Year==2020)  ## PP year
  yNA[testing]<-NA
  
  BetweenYear<-RunBGLR(YearEffects=FALSE,
                       nIter=80000,burnIn=60000,
                       y=yNA,testing=testing,
                       Inputfiledir=Inputfiledir,
                       Outputfiledir=savepath)
  
  r<-c(r,predict(testing=testing,gid=gid,yBLUE=yBLUE,Y=Y,fmfiledir=savepath))
}

write.csv(r, paste0(WD,savepath,"cor_",i,"Reps.csv"))
print(mean(r)) #-0.1159702



### 4. Yr20to19

yrTP<-2019
yrPP<-2020
RmCommonX<-droplevels(Y[Y$Year==yrPP,]$Crosses[which(Y[Y$Year==yrPP,]$Crosses%in%Y[Y$Year==yrTP,]$Crosses)])  # the CommonCrosses 
yBLUE<-Y[!Y$Crosses%in%RmCommonX,"BLUE_Trait"]
gid<-Y[!Y$Crosses%in%RmCommonX,"Crosses"]
y<-Y[!Y$Crosses%in%RmCommonX,"dryWgtPerM"]

setwd(paste0(WD,"OneTime1920/Yr20to19_output/"))

r<-NULL
for (i in 1:5){
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
}

for (i in 1:5){
  savepath<-paste0("OneTime1920/Yr20to19_output/Rep",i,"/")
  
  yBLUE2<-Y[!Y$Crosses%in%RmCommonX,"BLUE_Trait"]
  gid2<-Y[!Y$Crosses%in%RmCommonX,"Crosses"]
  #y2<-Y[!Y$Crosses%in%RmCommonX,"dryWgtPerM"]
  r<-c(r,predict(testing=testing,gid=gid2,yBLUE=yBLUE2,Y=Y,fmfiledir=savepath))
  
}

write.csv(r, paste0(WD,savepath,"cor_",i,"Reps.csv"))
print(mean(r))  # -0.02784699

