### ScenarioIII
#### Scenario with mixed=ploidy A matrix, but a simple GBLUP model
#### CV Within Yr2019

rm(list=ls())

WD<-"/local/workdir/mh865/GCA_SCA/GBLUPOnly_NoGCA_SCA/"  # run in terminal

datafdr<-"/local/workdir/mh865/GCA_SCA/OneTime1920/data/"

library(BGLR)

Yr<-2019

##!!! Load WithinYear dBLUPs Data
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata")) 
rownames(WithinYr_Both_dBLUPs)<-WithinYr_Both_dBLUPs$Row.names
WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,-1]

CrossBLUE<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==Yr,] ### Subset the Yr
CrossBLUE<-CrossBLUE[,!colnames(CrossBLUE)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")] ### RM extra cols

traits<-colnames(CrossBLUE)[!colnames(CrossBLUE)%in%c("Row.names","Crosses.x","Crosses.y","plotNo.x","plotNo.y","Year.x","Year.y")]
print(traits)

Y<-CrossBLUE

### Load the A matrix, outCovComb
load(paste0("/local/workdir/mh865/outCovComb/outCovComb4_Mix_Conden_0527_2021.Rdata")) ##!!!

A_866<-outCovComb4_MixOrder
IDs<-as.character(unique(Y$Crosses.x))  ### RM checks??
A<-as.matrix(A_866[IDs,IDs])

#
ETA<-list(list(K=A,model="RKHS"))

##### !!!


### Load Sample file
load(paste0(datafdr,"sampleCV_Yr",Yr,"_122Indiv_0312_2021.Rdata"))  
# 128 individuals in the sampleCV. But CrossBLUE only 127 lines
sampleCV<-sampleCV[-1,]  ##
folds<-1:10 
reps<-20 # !!!

ntraits<- length(traits) # !!!
cor<-matrix(nrow=reps,ncol=ntraits)
colnames(cor)<-traits


for (j in 1:length(traits)){
  Coltrait<-traits[j]
  
  y<-yBLUE<-Y[,colnames(Y)==Coltrait] 
  setwd(paste0(WD,"Alldata_CV_output/"))
  
  for (i in 1:reps){
    setwd(paste0(WD,"Alldata_CV_output/"))
    dir.create(paste0(Coltrait,"_OnlyYr",Yr,"_ydrBLUPsnoLocRep",i))
    savepath<-paste0("Alldata_CV_output/",Coltrait,"_OnlyYr",Yr,"_ydrBLUPsnoLocRep",i,"/")  # the path within WD!
    
    tmp<-NULL
    for (fold in folds){
      
      yNA<-y
      testing<-which(sampleCV[,i]==fold)
      yNA[testing]<-NA
      
      fm<-BGLR(y=yNA,
               ETA=ETA,
               nIter=80000,
               burnIn=60000,
               saveAt=paste0(WD,savepath,"CVData1920_",fold,"thfold_rep",i),
               verbose=TRUE)
      save(fm,file=paste0(WD,savepath,"fm_",fold,"thfold_rep",i,".rda"))
      
      yPred<-fm$ETA[[1]]$u  # equals to fm$yHat
      predict<-data.frame(testing,
                          Crosses=Y[c(testing),]$Crosses.x,
                          yBLUE=yBLUE[testing],
                          yPred=yPred[testing],
                          yHat=fm$yHat[testing],
                          plotNo=rownames(Y)[testing])
      # popChk=Y[c(testing),"popChk"],
      # Year=Y[c(testing),]$Year)
      predict<-droplevels(predict)
      
      tmp<-rbind(tmp,predict)
    }
    cor[i,j]<-cor(tmp$yBLUE,tmp$yPred,use="complete")
    
    write.table(tmp,file=paste0(WD,savepath,"predictions_rep",i,".csv"),row.names=FALSE,sep=",") 
    rm(fm) 
    unlink("*.dat")
  } 
  
}


standard_err<-function(x){sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))}
stderr<-apply(cor,2,standard_err)
cormean<-colMeans(cor)
cor_std<-rbind(cormean,stderr)
rownames(cor_std)<-c("corMean","StdErr")
colnames(cor_std)<-traits

write.csv(cor_std,paste0(paste0(WD,"Alldata_CV_output/","cor_CV_OnlyYr",Yr,"_ydrBLUPs_data_",length(traits),"Traits_Mean_05272021.csv")))
write.csv(cor,paste0(WD,"Alldata_CV_output/","cor_CV_OnlyYr",Yr,"_ydrBLUPs_data_",length(traits),"Traits_05272021.csv"))

### Warning
# Warning message:
#   In if (class(LT$K) != "matrix") stop("Kernel for linear term ",  :
#                                          the condition has length > 1 and only the first element will be used
#                                        
