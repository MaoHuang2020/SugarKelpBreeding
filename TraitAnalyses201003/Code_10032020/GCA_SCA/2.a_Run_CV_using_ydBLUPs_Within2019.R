# Run_CV_no_Loc_using_yBLUE_as_y_Within2019

rm(list=ls())

WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot

Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)

Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])

library(BGLR)



Y<-droplevels(Y[Y$Year==2019,])
Yr<-2019
folder<-"OneTime1920/Yr19Only/"
nreps<-13
subtract<-c(1:7)

Yr<-2019  ##!!!
folder<-"OneTime1920/Yr19Only/"

### Load Sample file
load(paste0(WD,"OneTime1920/data/sampleCV_Yr",Yr,"_122Indiv_0312_2021.Rdata"))

Inputfiledir<-c(paste0(folder,"GP1/"),paste0(folder,"GP2/"),paste0(folder,"GP1P2/")) 

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


#####!!!
datafdr<-paste0(WD,"OneTime1920/data/")
  ##!!! WithinYear Data
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata")) 
rownames(WithinYr_Both_dBLUPs)<-WithinYr_Both_dBLUPs$Row.names
WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,-1]

CrossBLUE<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==Yr,] ### Subset the Yr
CrossBLUE<-CrossBLUE[,!colnames(CrossBLUE)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")] ### RM extra cols

traits<-colnames(WithinYr_Both_dBLUPs)[!colnames(WithinYr_Both_dBLUPs)%in%c("Row.names","Crosses.x","Crosses.y","plotNo.x","plotNo.y","Year.x","Year.y")]
  print(traits)
##### !!!

sampleCV<-sampleCV
folds   <- 1:10 
reps<-20 # !!!
ntraits<- length(traits) # !!!
cor<-matrix(nrow=reps,ncol=ntraits)
colnames(cor)<-traits

for (j in 1:length(traits)){
  Coltrait<-traits[j]

 # Getting dBLUPs of 2019
    Y1<-Y[Y$Year==Yr,]
    
    CrossBLUEYr<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==Yr,]
    ## Look up the Crosses dBLUPs within Yr
    ## RM extra cols
    CrossBLUEYr<-CrossBLUEYr[,!colnames(CrossBLUEYr)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")]
    
    Y1$BLUE_Trait<-expss::vlookup(Y1$Crosses,dict=CrossBLUEYr,result_column = paste0(Coltrait),lookup_column = "Crosses.x")
 
 
  #Yrbind<-rbind(Y1,Y2)
  Ydata<-Y   # Save out the original Y, for in case use
  Y<-Y1  # Now Y is updated with its BLUEs from within 2019 Yr
  
  y<-yBLUE<-Y[,"BLUE_Trait"] 
  
  setwd(paste0(WD,"OneTime1920/Alldata_CV_output/"))

  for (i in 1:reps){
    setwd(paste0(WD,"OneTime1920/Alldata_CV_output/"))
    dir.create(paste0(Coltrait,"_OnlyYr",Yr,"_ydrBLUPsnoLocRep",i))
    savepath<-paste0("OneTime1920/Alldata_CV_output/",Coltrait,"_OnlyYr",Yr,"_ydrBLUPsnoLocRep",i,"/")  # the path within WD!
    
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
      
      yPred<-fm$ETA[[3]]$u+fm$ETA[[2]]$u+fm$ETA[[1]]$u  #SCA+GCA2+GCA1
      predict<-data.frame(testing,
                          Crosses=Y[c(testing),]$Crosses,
                          yBLUE=yBLUE[testing],
                          yPred=yPred[testing],
                          yHat=fm$yHat[testing],
                          popChk=Y[c(testing),"popChk"],
                          Year=Y[c(testing),]$Year)
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
  
write.csv(cor_std,paste0(paste0(WD,"OneTime1920/Alldata_CV_output/","cor_CV_OnlyYr",Yr,"_ydrBLUPs_data_",length(traits),"Traits_Mean_05272021.csv")))
write.csv(cor,paste0(WD,"OneTime1920/Alldata_CV_output/","cor_CV_OnlyYr",Yr,"_ydrBLUPs_data_",length(traits),"Traits_05272021.csv"))

