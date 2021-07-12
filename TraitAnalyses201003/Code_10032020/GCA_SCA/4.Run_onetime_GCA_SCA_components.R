### Use all data to run onetime GCA+SCA model, Get GCA and SCA var_component

##### GCA, SCA comp

rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot

Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])

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


#####!!!
datafdr<-paste0(WD,"OneTime1920/data/")

load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_overTwoYears_AddBD.Rdata")) ##!!!
rownames(Both_dBLUPs)<-Both_dBLUPs$Row.names
Both_dBLUPs<-Both_dBLUPs[,-1]
traits<-colnames(Both_dBLUPs)
  print(traits)
CrossBLUE<-Both_dBLUPs ##!!!
##### !!!

reps<-20
cor<-matrix(nrow=reps,ncol=length(traits)) ### i reps, j traits
  #colnames(cor)<-traits
corAll<-NULL
for (j in 1:length(traits)){
  Coltrait<-traits[j]
  Y$BLUE_Trait<-expss::vlookup(Y$Crosses,dict=CrossBLUE,result_column = paste0(Coltrait),lookup_column = "row.names") ### This is the DwPM
  head(Y)
  
  head(Y)
  dim(Y)
  
  Y<-Y   # Between Year !!!!
  y<-yBLUE<-Y[,"BLUE_Trait"]  # phenotypes column  !!!

#tmp<-NULL 

for (i in 1:reps){
  setwd(paste0(WD,"OneTime1920/Alldata_output/"))
  dir.create(paste0(Coltrait,"_OneTimeAll_Rep",i))
  savepath<-paste0("OneTime1920/Alldata_output/",Coltrait,"_OneTimeAll_Rep",i,"/")  # the path following WD
  
  y<-y  
  testing<-which(!is.na(Y$crossID))  # all lines
  str(testing)
  yNA<-y[testing]
  fm<-BGLR(y=yNA,
           ETA=ETA,
           nIter=80000,
           burnIn=60000,
           saveAt=paste0(WD,savepath,"OnetimeAll_rep",i),
           verbose=TRUE)
  save(fm,file=paste0(WD,savepath,"fm.rda"))
  
  yPred<-fm$ETA[[3]]$u+fm$ETA[[2]]$u+fm$ETA[[1]]$u  #SCA+GCA2+GCA1
  predict<-data.frame(testing,
                      Crosses=Y[c(testing),]$Crosses,
                      yBLUE=yBLUE[testing],
                      yPred=yPred[testing],
                      yHat=fm$yHat[testing],
                      popChk=Y[c(testing),"popChk"],
                      Year=Y[c(testing),]$Year)
  predict<-droplevels(predict)
  
  #tmp<-rbind(tmp,predict)
   
cor[i,j]<-cor(predict$yBLUE,predict$yPred,use="complete") 
  }
write.table(predict,file=paste0(WD,savepath,"Onetimepredictions_rep",i,".csv"),row.names=FALSE,sep=",")

corAll[[j]]<-cor
}

# standard_err<-function(x){sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))}
# stderr<-apply(cor,2,standard_err)
# cormean<-colMeans(cor)
# cor_std<-rbind(cormean,stderr)
# rownames(cor_std)<-c("corMean","StdErr")
# colnames(cor_std)<-traits


write.csv(cor,paste0(WD,"OneTime1920/Alldata_output/","Cor_predictItself_Alldata_",length(traits),"traits.csv")) 
write.csv(colMeans(cor),paste0(WD,"OneTime1920/Alldata_output/","Cor_predictItself_Alldata_",length(traits),"traits_Mean.csv")) 

####!!! AddBD






#### FUNCTION for finding 95% quantiles
Quan_fnc<-function(var=varU1){
  Uquan<-quantile(var,c(0.025,0.975))  #lower 0.025, upper 0.975
  Umean<-mean(var[var>Uquan[1]&var<Uquan[2]])
  Rangelow<-Uquan[1]
  Rangehigh<-Uquan[2]
  return(list(Umean=Umean,Rangelow=Rangelow,Rangehigh=Rangehigh))
}


##### Finding 95% CI
reps<-20
var95CIlow<-var95CIhigh<-GCASCA<-matrix(nrow=length(traits),ncol=4)

for (j in 1:length(traits)){
  Coltrait<-traits[j]
varfm<-NULL #showed 5 varcomp
varCIlow<-NULL
varCIhigh<-NULL
for (i in 1:reps){
# varcomp<-yHatVarMean_noloc(filedir = paste0("OneTime1920/Alldata_output/",Coltrait,"_OneTimeAll_Rep",i,"/"),vGCA1=1,vGCA2=2,vSCA=3,filename=paste0("OnetimeAll_rep",i))$varMean[2,]
# varfm<-rbind(varfm,varcomp)
  
  #varcomp2<-yHatVarMean_withloc(filedir = "OneTime1920/Alldata_output/withFMLocRep1/",vGCA1=1,vGCA2=2,vSCA=3,filename="Model_withFMLoc")
  
# 95% credible quantiles
  filedir = paste0("OneTime1920/Alldata_output/",Coltrait,"_OneTimeAll_Rep",i,"/");vGCA1=1;vGCA2=2;vSCA=3;filename=paste0("OnetimeAll_rep",i)
  load(paste0(WD,filedir,"fm.rda")) 
  fm<-fm
#  varfm<-c(fm$varE,fm$ETA[[vGCA1]]$varU,fm$ETA[[vGCA2]]$varU,fm$ETA[[vSCA]]$varU) # varcomp output by fm
  
  varE<-scan(paste0(WD,filedir,filename,"varE.dat"))
  # varB<-scan(paste0(WD,filedir,filename,"ETA_",vB,"_varB.dat"))
  varU1<-scan(paste0(WD,filedir,filename,"ETA_",vGCA1,"_varU.dat"))
  varU2<-scan(paste0(WD,filedir,filename,"ETA_",vGCA2,"_varU.dat"))
  varU3<-scan(paste0(WD,filedir,filename,"ETA_",vSCA,"_varU.dat"))
    #pdf(file="varU1_all_values_hist.pdf",width=1200,height = 800)
    #hist(varU1)
    #dev.off()
  Emean<-Quan_fnc(varE)$Umean
  U1mean<-Quan_fnc(varU1)$Umean
  U2mean<-Quan_fnc(varU2)$Umean
  U3mean<-Quan_fnc(varU3)$Umean
  
  Elow<-Quan_fnc(varE)$Rangelow
  U1low<-Quan_fnc(varU1)$Rangelow
  U2low<-Quan_fnc(varU2)$Rangelow
  U3low<-Quan_fnc(varU3)$Rangelow
    
  Ehigh<-Quan_fnc(varE)$Rangehigh
  U1high<-Quan_fnc(varU1)$Rangehigh
  U2high<-Quan_fnc(varU2)$Rangehigh
  U3high<-Quan_fnc(varU3)$Rangehigh
  
  varMean<-as.vector(c(Emean,U1mean,U2mean,U3mean)) #calculated varcomp by hand
  varlow<-as.vector(c(Elow,U1low,U2low,U3low))
  varhigh<-as.vector(c(Ehigh,U1high,U2high,U3high))
  names(varMean)<-names(varlow)<-names(varhigh)<-c("varE","varGCA1","varGCA2","varSCA")
  varfm<-rbind(varfm,varMean)
  varCIlow<-rbind(varCIlow,varlow)
  varCIhigh<-rbind(varCIhigh,varhigh)


}
  GCASCA[j,]<-colMeans(varfm)
  var95CIlow[j,]<-colMeans(varCIlow)
  var95CIhigh[j,]<-colMeans(varCIhigh)
}
colnames(GCASCA)<-colnames(var95CIlow)<-colnames(var95CIhigh)<-colnames(varfm)
rownames(GCASCA)<-rownames(var95CIlow)<-rownames(var95CIhigh)<-traits

write.csv(GCASCA,paste0(WD,"OneTime1920/Alldata_output/","OnetimeAll_",length(traits),"varcomp_95%CredibleInterval.csv"))

write.csv(var95CIlow,paste0(WD,"OneTime1920/Alldata_output/","OnetimeAll_",length(traits),"varcomp_95%CredibleInterval_Lowerbond.csv"))
write.csv(var95CIhigh,paste0(WD,"OneTime1920/Alldata_output/","OnetimeAll_",length(traits),"varcomp_95%CredibleInterval_Higherbond.csv"))




### Count the number of samples out of 16000 that's F>M
reps<-20
Fbig<-Mbig<-matrix(nrow=length(traits),ncol=1)
rownames(Fbig)<-rownames(Mbig)<-traits

for (j in 1:length(traits)){
  Coltrait<-traits[j]

FbiggerM<-NULL
MbiggerF<-NULL
  for (i in 1:reps){
    # varcomp<-yHatVarMean_noloc(filedir = paste0("OneTime1920/Alldata_output/",Coltrait,"_OneTimeAll_Rep",i,"/"),vGCA1=1,vGCA2=2,vSCA=3,filename=paste0("OnetimeAll_rep",i))$varMean[2,]
    # varfm<-rbind(varfm,varcomp)
    
    #varcomp2<-yHatVarMean_withloc(filedir = "OneTime1920/Alldata_output/withFMLocRep1/",vGCA1=1,vGCA2=2,vSCA=3,filename="Model_withFMLoc")
    
    # 95% credible quantiles
    filedir = paste0("OneTime1920/Alldata_output/",Coltrait,"_OneTimeAll_Rep",i,"/");vGCA1=1;vGCA2=2;vSCA=3;filename=paste0("OnetimeAll_rep",i)
    load(paste0(WD,filedir,"fm.rda")) 
    fm<-fm
    #  varfm<-c(fm$varE,fm$ETA[[vGCA1]]$varU,fm$ETA[[vGCA2]]$varU,fm$ETA[[vSCA]]$varU) # varcomp output by fm
    
    varE<-scan(paste0(WD,filedir,filename,"varE.dat"))
    # varB<-scan(paste0(WD,filedir,filename,"ETA_",vB,"_varB.dat"))
    varU1<-scan(paste0(WD,filedir,filename,"ETA_",vGCA1,"_varU.dat"))
    varU2<-scan(paste0(WD,filedir,filename,"ETA_",vGCA2,"_varU.dat"))
    varU3<-scan(paste0(WD,filedir,filename,"ETA_",vSCA,"_varU.dat"))
 
    FbiggerM<-rbind(FbiggerM,(sum(round(varU1,digits=4)>round(varU2,digits=4))/16000))  # 16000 samples
    MbiggerF<-rbind(MbiggerF,(sum(round(varU1,digits=4)<=round(varU2,digits=4))/16000))
  }
   Fbig[j,]<-colMeans(FbiggerM)
   Mbig[j,]<-colMeans(MbiggerF)
}

identical(rownames(Fbig),rownames(Mbig))
FMcompare<-cbind(Fbig,Mbig)
colnames(FMcompare)<-c("GCAF>M %","GCAM>F %")

write.csv(FMcompare,paste0(WD,"OneTime1920/Alldata_output/","GCA_varcomp_FM_compare%.csv"))





# 
# #### 4.FUNCTION 
# yHatVarMean_noloc<-function(filedir,vGCA1=1,vGCA2=2,vSCA=3,filename="OnetimeAll_rep1"){
#   load(paste0(WD,filedir,"fm.rda")) 
#   fm<-fm
#   varfm<-c(fm$varE,fm$ETA[[vGCA1]]$varU,fm$ETA[[vGCA2]]$varU,fm$ETA[[vSCA]]$varU) # varcomp output by fm
#   
#   varE<-scan(paste0(WD,filedir,filename,"varE.dat"))
#   # varB<-scan(paste0(WD,filedir,filename,"ETA_",vB,"_varB.dat"))
#   varU1<-scan(paste0(WD,filedir,filename,"ETA_",vGCA1,"_varU.dat"))
#   varU2<-scan(paste0(WD,filedir,filename,"ETA_",vGCA2,"_varU.dat"))
#   varU3<-scan(paste0(WD,filedir,filename,"ETA_",vSCA,"_varU.dat"))
#   
#   # plot(varU3,type='o',col=2,cex=.5)
#   #### Save the plot!!!!!
#   pdf(paste0(WD,filedir,"varE_varGCA1_2_SCA_OnetimeAll.pdf"))
#   par(mfrow=c(2,2))
#   plot(varE,type='o',col=2,cex=.5)
#   plot(varU1,type='o',col=1,cex=.5)
#   plot(varU2,type='o',col=1,cex=.5)
#   plot(varU3,type='o',col=1,cex=.5)
#   dev.off()
#   
#   varMedian<-c(median(varE),median(varU1),median(varU2),median(varU3))
#   varMean<-c(mean(varE),mean(varU1),mean(varU2),mean(varU3)) #calculated varcomp by hand
#   names(varMean)<-c("varE","varGCA1","varGCA2","varSCA")
#   
#   varMean<-rbind(varMean,varfm)     
#   write.csv(varMean,paste0(WD,filedir,filename,"_varcomp.csv"))
#   return(list(varMean=varMean))
# }