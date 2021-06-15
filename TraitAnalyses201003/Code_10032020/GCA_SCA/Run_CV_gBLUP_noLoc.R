### Make EVD data
rm(list=ls())

WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal
load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
#load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi
load(paste0(WD,"OneTime1920/data/","outCovComb_dip_0116_2021.Rdata"))
  #write.csv(outCovComb4_dipOrder,here("OneTime1920/data","A.csv"))
mm.file<- paste0(WD,"OneTime1920/data/","A.csv")          # path to covariates file
#source(paste0(WD,"Code/","BGLR_functions.R"))  #local
source(paste0(WD,"OneTime1920/","BGLR_functions_noFMLoc.R")) # terminal

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

    
# #############
# #the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
# CalGCA(Y=Y,mm.file=mm.file,colIDy = 3,colNam = "P1",savefiledir = "Onetime1920/Yr19_20/GP1/") 
# CalGCA(Y=Y,mm.file=mm.file,colIDy = 4,colNam = "P2",savefiledir = "Onetime1920/Yr19_20/GP2/") 
# 
# CalSCA(G1.file=paste0(WD,"Onetime1920/GP1/","G.rda"),
#        G2.file=paste0(WD,"Onetime1920/GP2/","G.rda"),
#        savefileDir="Onetime1920/GP1P2/")
# ############# Run it only once
# 


### 1. Both years data all used

r<-NULL 
for (i in 1:20){
  setwd(paste0(WD,"OneTime1920/Alldata_output/"))
  dir.create(paste0("noLocRep",i))
  savepath<-paste0("OneTime1920/Alldata_output/noLocRep",i,"/")  # the path following WD
  
y<-y  
testing<-which(!is.na(Y$crossID))  # all lines
  str(testing)
Alldata<-RunBGLR(YearEffects=TRUE,
              nIter=80000,burnIn=60000,
              y=y,testing=testing,
              Inputfiledir=Inputfiledir,
              Outputfiledir=savepath)  

  r<-c(r,predict(testing=testing,gid=gid,yBLUE=yBLUE,Y=Y,fmfiledir=savepath))
} 

write.csv(r, paste0(WD,savepath,"r_with_",i,"Reps.csv"))
mean(r)  # 0.936 over 5 reps

### 1.2 All data CV
# A. create sampling

#######
setwd(paste0(WD,"OneTime1920/Alldata_CV_output/"))

# sampleCV<-matrix(nrow=nrow(Y),ncol=500)
# for (n in 1:500){
#   sets<-rep(1:10,29)[-c(1:7)]  #283 ones
#   sampleCV[,n]<-sets[order(runif(nrow(Y)))]
# }
# save(sampleCV,file="sampleCV_283Indiv_0225_2021.Rdata")
# ###### One Time


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
  dir.create(paste0("noLocRep",i))
  savepath<-paste0("OneTime1920/Alldata_CV_output/noLocRep",i,"/")  # the path following WD
  
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

cor2<-matrix(nrow=reps,ncol=1)
for (i in 1:reps){
  savepath<-paste0("OneTime1920/Alldata_CV_output/noLocRep",i,"/")
  predict<-read.csv(paste0(WD,savepath,"predictions_rep",i,".csv"),sep=",",header=TRUE)
  cor2[i,]<-cor(predict[predict$popChk=="ES",]$yBLUE,predict[predict$popChk=="ES",]$yPred,use="complete")
}
  colMeans(cor2) #0.1827861     #noLOC 0.2699811
write.csv(cor2,"10foldCV_cor2_using_250ES_plots_NoFMLoc.csv")

##### GCA, SCA comp
# reps<-20
# varfm<-NULL #showed 5 varcomp
# for (i in 1:reps){
#   varcomp<-yHatVarMean_noloc(filedir = paste0("OneTime1920/Alldata_output/noLocRep",reps,"/"),vB=2,vGCA1=3,vGCA2=4,vSCA=5,filename="Model_noFMLoc")$varMean[2,]
#   varfm<-rbind(varfm,varcomp)
#   
#   varcomp2<-yHatVarMean_withloc(filedir = "OneTime1920/Alldata_output/withFMLocRep1/",vB=3,vGCA1=4,vGCA2=5,vSCA=6,filename="Model_withFMLoc")
# }
# colMeans(varfm)

