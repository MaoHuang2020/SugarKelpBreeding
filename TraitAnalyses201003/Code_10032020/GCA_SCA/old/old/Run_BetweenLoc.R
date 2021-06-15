# Between Loc Run
rm(list=ls())
#library(here) 
# I quit using this library. Its source code had some issues at one point and can't read in

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

locs<-c("CB","CC","JS","LD","NC","NL","OI","SF")

for (i in 1:5){
  setwd(paste0(WD,"OneTime1920/BetweenLoc_output/")) #where to create Rep# folder
  dir.create(paste0("Rep",i))
  WDloc<-paste0(WD,"OneTime1920/BetweenLoc_output/Rep",i,"/")  # the path following WD
  
  setwd(WDloc)   #where to create loc# folder
  r<-NULL
  for(loc in locs) {
    
    dir.create(paste0("loc",loc))
    savepath<-paste0("OneTime1920/BetweenLoc_output/Rep",i,"/loc",loc,"/")
    
    yNA<-y 
    testing<-which(Y$femaParLoc==loc)  # 156 plots
    yNA[testing]<-NA  # Other locations to predict this testing one
    
    BetweenLoc<-RunBGLR(YearEffects=TRUE,
                        nIter=80000,burnIn=60000,
                        y=yNA,testing=testing,
                        Inputfiledir=Inputfiledir,
                        Outputfiledir=savepath) 
    
    r<-c(r,predict(testing=testing,gid=gid,yBLUE=yBLUE,Y=Y,fmfiledir=savepath))
    
  }
  
  names(r)<-locs
  write.csv(r, paste0("r_1Loc_PP_Rep",i,".csv"))
}


rAll<-NULL
for (i in 1:5){
  
  WDloc<-paste0(WD,"OneTime1920/BetweenLoc_output/Rep",i,"/")  # the path following WD
  setwd(WDloc)
  r<-NULL
  
  for(loc in locs){
    #yNA<-y 
    testing<-which(Y$femaParLoc==loc)  # 156 plots
    #yNA[testing]<-NA
    savepath<-paste0("OneTime1920/BetweenLoc_output/Rep",i,"/loc",loc,"/")
    
    r<-c(r,predict(testing=testing,gid=gid,yBLUE=yBLUE,Y=Y,fmfiledir=savepath)) 
  }
  names(r)<-locs
  write.csv(r, paste0("r_1Loc_PP_Rep",i,".csv"))
  
  rAll[[i]]<-read.csv(paste0(WD,"OneTime1920/BetweenLoc_output/Rep",i,"/","/r_1Loc_PP_Rep",i,".csv"),row.names=1)
  #std error of the BLUEs??
}  

rowMeans(do.call(cbind.data.frame, rAll))

