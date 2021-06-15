#Run_Between Yr_using_yBLUE_as_y.R
######## This is only for ONE Trait ---- OLD

#Run_CV_using_yBLUEs as predicted y. Work with 250 plots only
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

#load(paste0(WD,"OneTime1920/data/","outCovComb_dip_0116_2021.Rdata"))
#write.csv(outCovComb4_dipOrder,here("OneTime1920/data","A.csv"))
mm.file<- paste0(WD,"OneTime1920/data/","A.csv")          # path to covariates file

### This is to only get the calc GCA and SCA functions
#source(paste0(WD,"OneTime1920/code/","BGLR_functions_noFMLoc.R")) # !!!terminal
### Run it only once

### Get the "BLUE_Trait"
load(paste0(WD,"OneTime1920/data/","BLUE_DwPM_2vs1Year_Update03082021.rdata"))

library(expss)
##!!!!!!!!!! BLUE is the 2019 !

 for (Yr in c(2019)) {
    Y1<-Y[Y$Year==Yr,]
    Y1$BLUE_Trait<-vlookup(Y1$plotNo,dict=CrossBLUE,result_column = "BLUE_DwPM_2019",lookup_column = "plotNo")  # This is the DwPM
  } 

 for (Yr in c(2020)){
    Y2<-Y[Y$Year==Yr,]
    Y2$BLUE_Trait<-vlookup(Y2$plotNo,dict=CrossBLUE,result_column = "BLUE_DwPM_2020",lookup_column = "plotNo") 
}

Yrbind<-rbind(Y1,Y2)
Ydata<-Y
Y<-Yrbind



#plot<-ggplot(data=Y,aes(BLUE_Trait,Year))+
#  geom_point(aes(color=as.factor(Year)))+ 
#  geom_line(aes(group=as.factor(Crosses)))
#print(plot)

### Run it only once
# # #the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
# CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = 3,colNam = "P1",savefiledir = "OneTime1920/GP1/250Individual/") 
# CalGCA_noChk(Y=Y,mm.file=mm.file,colIDy = 5,colNam = "P2",savefiledir = "OneTime1920/GP2/250Individual/") 
# # 
#  CalSCA(G1.file=paste0(WD,"OneTime1920/GP1/250Individual/","G.rda"),
#         G2.file=paste0(WD,"OneTime1920/GP2/250Individual/","G.rda"),
#         savefileDir="OneTime1920/GP1P2/250Individual/")
# # ############# Run it only once
# # 
Inputfiledir<-c("OneTime1920/GP1/250Individual/","OneTime1920/GP2/250Individual/","OneTime1920/GP1P2/250Individual/") 
  head(Y)
  dim(Y)
  

load(paste0(WD,Inputfiledir[1],"EVD.rda"))  
EVD1<-EVD
  rm(EVD)
load(paste0(WD,Inputfiledir[2],"EVD.rda"))       
EVD2<-EVD
  rm(EVD)
load(paste0(WD,Inputfiledir[3],"EVD.rda"))       
EVD3<-EVD
  rm(EVD)
ETA<-list(
  list(V=EVD1$vectors,d=EVD1$values,model="RKHS"),
  list(V=EVD2$vectors,d=EVD2$values,model="RKHS"),
  list(V=EVD3$vectors,d=EVD3$values,model="RKHS")
)

reps<-20 # !!!
ntraits<-1  # !!!
 head(Y)  

# gid<-Y[,"Crosses"]   # Crosses
y<-yBLUE<-Y[,"BLUE_Trait"] # This is the BLUE for DwPM

cor<-matrix(nrow=reps,ncol=ntraits)
for (i in 1:reps){
  setwd(paste0(WD,"OneTime1920/Yr19to20_output/"))
  
  dir.create(paste0("yBLUEnoLocRep",i))
  savepath<-paste0("OneTime1920/Yr19to20_output/yBLUEnoLocRep",i,"/")
  
    yNA<-y
    # print(fold)
    testing<-which(Y$Year==2020)  ## PP year
    yNA[testing]<-NA
    
    fm<-BGLR::BGLR(y=yNA,
             ETA=ETA,
             nIter=80000,
             burnIn=60000,
             saveAt=paste0(WD,savepath,"Yr19To20_rep",i),
             verbose=TRUE)
    save(fm,file=paste0(WD,savepath,"fm_Yr19To20_rep",i,".rda"))
    
    yPred<-fm$ETA[[3]]$u+fm$ETA[[2]]$u+fm$ETA[[1]]$u  #SCA+GCA2+GCA1
    predict<-data.frame(testing,
                        Crosses=Y[c(testing),]$Crosses,
                        yBLUE=yBLUE[testing],
                        yPred=yPred[testing],
                        yHat=fm$yHat[testing],
                        popChk=Y[c(testing),"popChk"],
                        Year=Y[c(testing),]$Year)
    predict<-droplevels(predict)
    cor[i,]<-cor(predict$yBLUE,predict$yPred,use="complete")
    write.table(predict,file=paste0(WD,savepath,"predictions_rep",i,".csv"),row.names=FALSE,sep=",") 
  }
  
 
  rm(fm) 
  unlink("*.dat")

colMeans(cor)  # All 283 individuals from tmp   #noLOC 0.2696567
write.csv(cor,"cor_Yr19_predict_20_using_250ES_plots_NoFMLoc_yBLUEas_y.csv")



## This is to confirm the cors using "predict" data
# cor2<-matrix(nrow=reps,ncol=1)
# for (i in 1:reps){
#   savepath<-paste0("OneTime1920/Yr19to20_output/yBLUEnoLocRep",i,"/")
#   predict<-read.csv(paste0(WD,savepath,"predictions_rep",i,".csv"),sep=",",header=TRUE)
#   cor2[i,]<-cor(predict[predict$popChk=="ES",]$yBLUE,predict[predict$popChk=="ES",]$yPred,use="complete")
# }
# colMeans(cor2) 
# write.csv(cor2,"cor_Yr19_predict_20_using_250ES_plots_NoFMLoc_yBLUEas_y_confirm.csv")

