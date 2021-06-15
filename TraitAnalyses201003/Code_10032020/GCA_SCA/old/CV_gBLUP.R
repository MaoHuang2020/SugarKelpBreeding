#CV rrBLUP
### Make EVD data
rm(list=ls())
#library(here) 
# I quit using this library. Its source code had some issues at one point and can't read in
# here()

#WD<-"/Users/maohuang/Desktop/Kelp/GCA_SCA/" # local

WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal
load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
#load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi
load(paste0(WD,"OneTime1920/data/","outCovComb_dip_0116_2021.Rdata"))
#write.csv(outCovComb4_dipOrder,here("OneTime1920/data","A.csv"))
mm.file<- paste0(WD,"OneTime1920/data/","A.csv")          # path to covariates file
#source(paste0(WD,"Code/","BGLR_functions.R"))  #local
source(paste0(WD,"OneTime1920/","BGLR_functions.R")) # terminal

####
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years

### Get the "BLUE_Trait"
load(paste0(WD,"OneTime1920/data/BLUE_DwPM_2vs1Year.rdata"))

library(expss)
Y$BLUE_Trait<-vlookup(Y$Crosses,dict=CrossBLUE,result_column = "BLUE_DwPM_2Yrs",lookup_column = "CrossName") ### This is the DwPM
head(Y)

plot<-ggplot(data=Y,aes(BLUE_Trait,Year))+
  geom_point(aes(color=as.factor(Year)))+ 
  geom_line(aes(group=as.factor(Crosses)))
print(plot)

###
Inputfiledir<-c("OneTime1920/GP1/","OneTime1920/GP2/","OneTime1920/GP1P2/") 
head(Y)
dim(Y)
y<-Y[,"dryWgtPerM"]  # phenotypes column  !!!
gid<-Y[,"Crosses"]   # Crosses
yBLUE<-Y[,"BLUE_Trait"] # This is the BLUE for DwPM



### 1.2 All data CV
# A. create sampling

#######
setwd(paste0(WD,"OneTime1920/Alldata_CV_output/"))

sampleCV<-matrix(nrow=nrow(Y),ncol=500)
for (n in 1:500){
  sets<-rep(1:10,29)[-c(1:7)]  #283 ones
  sampleCV[,n]<-sets[order(runif(nrow(Y)))]
}
save(sampleCV,file="sampleCV_283Indiv_0225_2021.Rdata")
###### One Time


load(paste0(WD,"OneTime1920/Alldata_CV_output/sampleCV_283Indiv_0225_2021.Rdata"))

sampleCV<-sampleCV
folds   <- 1:10 
reps<-20 # !!!
ntraits<-1  # !!!
cor<-matrix(nrow=reps,ncol=ntraits)
#          list(~factor(femaParLoc)+factor(maleParLoc),data=Y,model="BRR"),
load(paste0(WD,Inputfiledir[1],"EVD.rda"))  
EVD1<-EVD
rm(EVD)
load(paste0(WD,Inputfiledir[2],"EVD.rda"))       
EVD2<-EVD
rm(EVD)
load(paste0(WD,Inputfiledir[3],"EVD.rda"))       
EVD3<-EVD
rm(EVD)
ETA<-list(list(~factor(popChk)+factor(Year),data=Y,model="FIXED"),
          list(~factor(line)+factor(block),data=Y,model="BRR"),
          
          list(V=EVD1$vectors,d=EVD1$values,model="RKHS"),
          list(V=EVD2$vectors,d=EVD2$values,model="RKHS"),
          list(V=EVD3$vectors,d=EVD3$values,model="RKHS")
)

for (i in 1:reps){
  setwd(paste0(WD,"OneTime1920/Alldata_CV_output/"))
  dir.create(paste0("Rep",i))
  savepath<-paste0("OneTime1920/Alldata_CV_output/Rep",i,"/")  # the path following WD
  
  tmp<-NULL
  for (fold in folds){
    
    #dir.create(paste0('10folds_Cycle',i,"/"))
    #savepath<-paste0('10folds_Cycle',i,"/")   ### Creating the fold_# folder
    
    yNA<-y
    # print(fold)
    testing<-which(sampleCV[,i]==fold)
    yNA[testing]<-NA
    
    fm<-BGLR(y=yNA,
             ETA=ETA,
             nIter=80000,
             burnIn=60000,
             saveAt=paste0(WD,savepath,"CVData1920_",fold,"thfold_rep",i),
             verbose=TRUE)
    save(fm,file=paste0(WD,savepath,"fm_",fold,"thfold_rep",i,".rda"))
    
    yPred<-fm$ETA[[5]]$u+fm$ETA[[4]]$u+fm$ETA[[3]]$u 
    predict<-data.frame(testing,Crosses=gid[testing],yBLUE=yBLUE[testing],yPred=yPred[testing])
    predict$popChk<-expss::vlookup(predict$Crosses,dict=Y,lookup_column = "Crosses",result_column = "popChk")
    predict$Year<-expss::vlookup(predict$Crosses,dict=Y,lookup_column = "Crosses",result_column = "Year")
    
    tmp<-rbind(tmp,predict)
    # r<-c(r,predict(testing=testing,gid=gid,yBLUE=yBLUE,Y=Y,fmfiledir=savepath))
  }
  cor[i,]<-cor(tmp$yBLUE,tmp$yPred,use="complete")
  
  write.table(tmp,file=paste0(WD,savepath,"predictions_rep",i,".csv"),row.names=FALSE,sep=",") 
  rm(fm) 
  unlink("*.dat")
}
colMeans(cor)  # All 283 individuals from tmp   #noLOC 0.2696567
write.csv(cor,"10foldCV_cor2_using_283_plots_NoFMLoc.csv")