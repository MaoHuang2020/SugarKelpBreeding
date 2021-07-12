# Run Yr19 predict Yr 20

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


pop<-"Yr19to20"
TPYr<-2019      # TPYr   
NAYr<-2020 
#####!!!
datafdr<-paste0(WD,"OneTime1920/data/")
##!!!WithinYr Data
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata")) 
WithinYr_Both_dBLUPs<-droplevels(WithinYr_Both_dBLUPs)

identical(as.character(WithinYr_Both_dBLUPs$Row.names),as.character(WithinYr_Both_dBLUPs$plotNo.y))

rownames(WithinYr_Both_dBLUPs)<-WithinYr_Both_dBLUPs$Row.names
WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,-1]
head(WithinYr_Both_dBLUPs)

traits<-colnames(WithinYr_Both_dBLUPs)[!colnames(WithinYr_Both_dBLUPs)%in%c("Row.names","Crosses.x","Crosses.y","plotNo.x","plotNo.y","Year.x","Year.y")]
print(traits)
##### !!!

reps<-20 # !!!
ntraits<-length(traits)  # !!!
head(Y)  

cor<-matrix(nrow=reps,ncol=length(traits))
colnames(cor)<-traits

for (j in 1:length(traits)){
  Coltrait<-traits[j]
  
  #!!!!!!!!!! USE drBLUP estimated WithinYra 2019 and 2020 !
  for (Yr in c(2019)) {
    Y1<-Y[Y$Year==Yr,]
    
    CrossBLUEYr<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==Yr,]
    ## Look up the Crosses dBLUPs within Yr
    ## RM extra cols
    CrossBLUEYr<-CrossBLUEYr[,!colnames(CrossBLUEYr)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")]
    
    Y1$BLUE_Trait<-expss::vlookup(Y1$Crosses,dict=CrossBLUEYr,result_column = paste0(Coltrait),lookup_column = "Crosses.x")
  }
  
  for (Yr in c(2020)){
    Y2<-Y[Y$Year==Yr,]
    
    CrossBLUEYr<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==Yr,]
    
    CrossBLUEYr<-CrossBLUEYr[,!colnames(CrossBLUEYr)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")]
    
    Y2$BLUE_Trait<-expss::vlookup(Y2$Crosses,dict=CrossBLUEYr,result_column = paste0(Coltrait),lookup_column = "Crosses.x")
  }
  
  Yrbind<-rbind(Y1,Y2)
  Ydata<-Y   # Save out the original Y, for in case use
  Y<-Yrbind  # Now Y is updated with its BLUEs from within 2019 and 2020 WithinYr
  
  y<-yBLUE<-Y[,"BLUE_Trait"]   # phenotypes column  !!!
  
  for (i in 1:reps){
    setwd(paste0(WD,"OneTime1920/",pop,"_output/"))
    
    dir.create(paste0(Coltrait,"_ydrBLUP_noLocRep",i))
    savepath<-paste0("OneTime1920/",pop,"_output/",Coltrait,"_ydrBLUP_noLocRep",i,"/")
    
    yNA<-y
    
    testing<-which(Y$Year==NAYr)  ## PP year
    yNA[testing]<-NA
    
    fm<-BGLR::BGLR(y=yNA,
                   ETA=ETA,
                   nIter=80000,
                   burnIn=60000,
                   saveAt=paste0(WD,savepath,pop,"_rep",i),
                   verbose=TRUE)
    save(fm,file=paste0(WD,savepath,"fm_",pop,"_rep",i,".rda"))
    
    yPred<-fm$ETA[[3]]$u+fm$ETA[[2]]$u+fm$ETA[[1]]$u  #SCA+GCA2+GCA1
    predict<-data.frame(testing,
                        Crosses=Y[c(testing),]$Crosses,
                        yBLUE=yBLUE[testing],
                        yPred=yPred[testing],
                        yHat=fm$yHat[testing],
                        popChk=Y[c(testing),"popChk"],
                        Year=Y[c(testing),]$Year)
    predict<-droplevels(predict)
    cor[i,j]<-cor(predict$yBLUE,predict$yPred,use="complete")
    write.table(predict,file=paste0(WD,savepath,Coltrait,"_predictions_rep",i,".csv")) 
  }
  rm(fm) 
  unlink("*.dat")
  print(colMeans(cor))  
}

standard_err<-function(x){sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))}
stderr<-apply(cor,2,standard_err)
cormean<-colMeans(cor)
cor_std<-rbind(cormean,stderr)
rownames(cor_std)<-c("corMean","StdErr")
colnames(cor_std)<-traits


write.csv(cor,paste0(WD,"OneTime1920/",pop,"_output/","cor_",pop,"_",length(traits),"_traits_","_using_250ES_ydrBLUP_NoLoc_05272021.csv"))
write.csv(cor_std,paste0(WD,"OneTime1920/",pop,"_output/","cor_",pop,"_",length(traits),"_traits_","_using_250ES_ydrBLUP_NoLoc_Mean_05272021.csv"))






# This is to confirm the cors using "predict" data
# cor2<-matrix(nrow=reps,ncol=1)
# for (i in 1:reps){
#   savepath<-paste0("OneTime1920/Yr20to19_output/yBLUEnoLocRep",i,"/")
#   predict<-read.csv(paste0(WD,savepath,"predictions_rep",i,".csv"),sep=",",header=TRUE)
#   cor2[i,]<-cor(predict[predict$popChk=="ES",]$yBLUE,predict[predict$popChk=="ES",]$yPred,use="complete")
# }
# colMeans(cor2)
# write.csv(cor2,"cor_Yr19_predict_20_using_250ES_plots_NoFMLoc_yBLUEas_y_confirm.csv")



# ### Get the "BLUE_Trait"
# load(paste0(WD,"OneTime1920/data/","BLUE_DwPM_2vs1Year_Update03082021.rdata"))
# 
# ##!!!!!!!!!! BLUE is estimated within 2019 and 2020 !
# 
# for (Yr in c(2019)) {
#   Y1<-Y[Y$Year==Yr,]
#   Y1$BLUE_Trait<-expss::vlookup(Y1$Crosses,dict=CrossBLUE,result_column = "BLUE_DwPM_2019",lookup_column = "CrossName")  # This is the DwPM
# } 
# 
# for (Yr in c(2020)){
#   Y2<-Y[Y$Year==Yr,]
#   Y2$BLUE_Trait<-expss::vlookup(Y2$Crosses,dict=CrossBLUE,result_column = "BLUE_DwPM_2020",lookup_column = "CrossName") 
# }
# 
# Yrbind<-rbind(Y1,Y2)
# Ydata<-Y
# Y<-Yrbind



