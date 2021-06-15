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

    
#############
#the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
CalGCA(Y=Y,mm.file=mm.file,colIDy = 3,colNam = "P1",savefiledir = "Onetime1920/Yr19_20/GP1/") 
CalGCA(Y=Y,mm.file=mm.file,colIDy = 4,colNam = "P2",savefiledir = "Onetime1920/Yr19_20/GP2/") 

CalSCA(G1.file=paste0(WD,"Onetime1920/GP1/","G.rda"),
       G2.file=paste0(WD,"Onetime1920/GP2/","G.rda"),
       savefileDir="Onetime1920/GP1P2/")
############# Run it only once



### 1. Both years data all used

r<-NULL 
for (i in 1:5){
  setwd(paste0(WD,"OneTime1920/Alldata_output/"))
  dir.create(paste0("Rep",i))
  savepath<-paste0("OneTime1920/Alldata_output/Rep",i,"/")  # the path following WD
  
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
reps<-5 # !!!
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

cor2<-matrix(nrow=reps,ncol=1)
for (i in 1:reps){
  savepath<-paste0("OneTime1920/Alldata_CV_output/Rep",i,"/")
  predict<-read.csv(paste0(WD,savepath,"predictions_rep",i,".csv"),sep=",",header=TRUE)
  cor2[i,]<-cor(predict[predict$popChk=="ES",]$yBLUE,predict[predict$popChk=="ES",]$yPred,use="complete")
}
  colMeans(cor2) #0.1827861     #noLOC 0.2699811
write.csv(cor2,"10foldCV_cor2_using_250ES_plots_NoFMLoc.csv")

### 2. Between Location prediction  

### RM one loc at a time, as testing set
table(droplevels(Y)$femaParLoc)
locs<-levels(droplevels(Y)$femaParLoc)
v<-NULL
for (j in 1:length(locs))
{
  uniqueCrossesPerLoc<-count(unique(droplevels(Y)[Y$femaParLoc==locs[j],]$Crosses))
  c<-nrow(uniqueCrossesPerLoc)
  v<-c(v,c)
}
v
names(v)<-locs
#Total unique Crosses with data
#CB CC JS LD LL ME NC NL OD OI SF 
#14 46 14 40  3  2 15 38  2 54 20

#Total plots with data
#CB CC JS LD LL ME NC NL OD OI SF 
#14 47 14 42  3 17 15 45  2 64 20


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

#One rep only, r made on all plots (not on Crosses)
  # CB          CC          JS          LD          ME(check)          NC
  # 0.54064825 -0.08664727 -0.65768603  0.66493417  0.99995689 -0.33652471
  # NL          OI          SF
  # 0.14753220  0.57215122  0.64058045
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


### std error of the BLUEs??
stdE<-NULL
mean<-NULL
for (loc in locs){
  # Read in the BLUE/predict file in Rep 1
  i=1
  BLUEpath<-paste0(WD,"OneTime1920/BetweenLoc_output/Rep",i,"/loc",loc,"/")
  x<-read.csv(paste0(BLUEpath,"/TP_predict_PP.csv"),sep=",",header=T)$yBLUE
  std <- function(x) sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))  # sdev/sqrt(n)
  stdE<-c(stdE,std(x))
  mean<-c(mean,mean(x,na.rm=TRUE))
}
  names(stdE)<-names(mean)<-locs
  
  
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
    r<-c(r,predict(testing=testing,gid=gid2,yBLUE=yBLUE2,Y=Y,fmfiledir=savepath))
  }
  
  
  # for (i in 1:5){
  #   savepath<-paste0("OneTime1920/Yr20to19_output/Rep",i,"/")
  #   
  #   yBLUE2<-Y[!Y$Crosses%in%RmCommonX,"BLUE_Trait"]
  #   gid2<-Y[!Y$Crosses%in%RmCommonX,"Crosses"]
  #   #y2<-Y[!Y$Crosses%in%RmCommonX,"dryWgtPerM"]
  # 
  #   
  # }
  
  write.csv(r, paste0(WD,savepath,"cor_",i,"Reps.csv"))
  print(mean(r))  # -0.02784699
  
  #####
  yHatVarMean(filedir = "OneTime1920/Alldata_output/Rep1/") # the directory that follows here()
  yHatVarMean(filedir = "OneTime1920/Yr19to20_output/Rep1/") # the directory that follows here()
  yHatVarMean(filedir = "OneTime1920/Yr20to19_output/") # the directory that follows here()
  BetweenYear$fm$varE;BetweenYear$fm$ETA[[2]]$varB;BetweenYear$fm$ETA[[3]]$varU;BetweenYear$fm$ETA[[4]]$varU;BetweenYear$fm$ETA[[5]]$varU
  yHatVarMean(filedir = "OneTime1920/BetweenLoc_output/")
  yHatVarMean(filedir = "OneTime1920/BothGPgeno_output/")
  yHatVarMean(filedir = "OneTime1920/FGgeno_output/Rep1/")
  #yHatVarMean(filedir = "OneTime1920/test_output/")
  
  varU<-NULL
  for (i in 1:5){
    var<-yHatVarMean(filedir = paste0("OneTime1920/Alldata_output/Rep",i,"/"))
    varU<-rbind(varU,var[2,])
    colnames(varU)<-colnames(var)
  }
colMeans(varU)


varU<-NULL
for (i in 1:5){
  var<-yHatVarMean(filedir = paste0("OneTime1920/BothGPgeno_output/Rep",i,"/"))
  varU<-rbind(varU,var[2,])
  colnames(varU)<-colnames(var)
}
colMeans(varU)

  


### 5. Both GPs genotyped,predict others  

# Input the list of GPs being genotyped
GPsequenced<-read.csv(paste0(WD,"OneTime1920/data/","3_plates_sequenced_names.csv"),sep=",",header=TRUE)  
dim(GPsequenced)
head(GPsequenced)
GPsequenced<-GPsequenced$Name_InPheno 
str(GPsequenced)

r<-NULL
for (i in 1:5){
  setwd(paste0(WD,"OneTime1920/BothGPgeno_output/"))
  dir.create(paste0("Rep",i))
  savepath<-paste0("OneTime1920/BothGPgeno_output/Rep",i,"/")  # the path following WD
  
  yNA<-y
  BothGP<-which(Y[Y$femaPar%in%GPsequenced,]$malePar%in%GPsequenced)  # 166 plots
  #
  OneFG<-which(Y$femaPar%in%GPsequenced & !Y$malePar%in%GPsequenced) # 49
  OneMG<-which(Y$malePar%in%GPsequenced & !Y$femaPar%in%GPsequenced) # 54
  NeitherGP<-which(!Y$femaPar%in%GPsequenced & !Y$malePar%in%GPsequenced) #14
  
  Y[NeitherGP,][,1:10]
  
  testing<-setdiff(1:nrow(Y),BothGP) #117 plots
  yNA[testing]<-NA
  
  BothGPgeno<-RunBGLR(YearEffects=TRUE,
                      nIter=80000,burnIn=60000,
                      y=yNA,testing=testing,
                      Inputfiledir=Inputfiledir,
                      Outputfiledir=savepath) 
  
  r<-c(r,predict(testing=testing,gid=gid,yBLUE=yBLUE,Y=Y,fmfiledir=savepath))
}
write.csv(r, paste0(WD,savepath,"cor_",i,"Reps.csv"))
print(mean(r))  # -0.0700428


### 6. SPs, one genotyped, another not, predict those both are not
r<-NULL
for (i in 1:5){
  setwd(paste0(WD,"OneTime1920/FGgeno_output/"))
  dir.create(paste0("Rep",i))
  savepath<-paste0("OneTime1920/FGgeno_output/Rep",i,"/")  # the path following WD
  
  yNA<-y
  #
  OneFG<-which(Y$femaPar%in%GPsequenced & !Y$malePar%in%GPsequenced)  #49
  
  testing<-setdiff(1:nrow(Y),OneFG)   #234
  yNA[testing]<-NA
  
  OneGPgeno<-RunBGLR(YearEffects=TRUE,
                     nIter=80000,burnIn=60000,
                     y=yNA,testing=testing,
                     Inputfiledir=Inputfiledir,
                     Outputfiledir=savepath) 
  
  r<-c(r,predict(testing=testing,gid=gid,yBLUE=yBLUE,Y=Y,fmfiledir=savepath))
}
write.csv(r, paste0(WD,savepath,"cor_",i,"Reps.csv"))
print(mean(r)) #0.07600304


###############################################  
#### Making BLUEs within Year and between Years
  rm(list=ls())
  library(here) 
  # I quit using this library. Its source code had some issues at one point and can't read in
  #here()
  
  WD<-"/Users/maohuang/Desktop/Kelp/GCA_SCA/" # local

  load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
  dataNHpi<-dataNHpiBoth_C  ##!!!!!
  
  library(lme4)
  exptlSP <- as.character(dataNHpi$popChk) == "ES"
  dataNHpi$entry <- as.character(dataNHpi$popChk)
  dataNHpi$entry[exptlSP] <- as.character(dataNHpi$plotNo[exptlSP])
  dataNHpi$group <- as.factor(ifelse(exptlSP, 1, 0))
  
  library(lmerTest) # can give p-values
  
  for (col in c( "line", "block","popChk","group","entry","Year")) 
    dataNHpi[,col] <- factor(dataNHpi[,col])
    dataNHpi$Trait<-dataNHpi$dryWgtPerM  ### !!!!! TraitBLUE
  
  data<-dataNHpi
  fitAug <- lm(Trait ~ Year+ line*Year+ block*Year+ group +Crosses , data=data) # Crosses:group
    summary(fitAug)
  colnames(summary(fitAug)$coef)
  CrossBLUE<-summary(fitAug)$coef[grep("x",rownames(summary(fitAug)$coef)),"Estimate"]  # "x": xtract the Crosses
    str(CrossBLUE)
  CrossBLUE<-as.data.frame(CrossBLUE)
  library(stringr)
  CrossBLUE$CrossName<-str_replace(rownames(CrossBLUE),"Crosses","")
    head(CrossBLUE)
  colnames(CrossBLUE)[colnames(CrossBLUE)=="CrossBLUE"]<-"BLUE_DwPM_2Yrs"
  save(CrossBLUE,file=paste0(WD,"OneTime1920/data/","BLUE_DwPM_2vs1Year.rdata"))
  
  library(expss)
  CrossBLUE$plotNo<-vlookup(CrossBLUE$CrossName,dict=data,result_column = "plotNo",lookup_column="Crosses")
  CrossBLUE$Year<-vlookup(CrossBLUE$CrossName,dict=data,result_column = "Year",lookup_column="Crosses")
    head(CrossBLUE)
  which(!is.na(CrossMerge$TraitBLUE)) 
  
  which(!unique(droplevels(dataNHpiBoth_C)$Crosses)%in%CrossBLUE$CrossName)
  dataNHpiBoth_C[c(1,205,247),]  ####??????Why these three crosses got thrown out??????
  

  
### Within Year BLUEs, This becomes not necessary, because the Year effects are fixed 
  # library(stringr) 
  # CrossBLUE2<-NULL
  # for (Year in c(2019,2020)){
  #   i=Year-2018
  #   data<-droplevels(dataNHpi[dataNHpi$Year==Year,])
  #   fitAug1 <- lm(Trait ~ line+ block+ group+Crosses , data=data) # Crosses:group
  #   CrossBLUE2[[i]]<-summary(fitAug)$coef[grep("x",rownames(summary(fitAug)$coef)),"Estimate"]
  #   CrossBLUE2[[i]]<-as.data.frame(CrossBLUE2[[i]])
  #   CrossBLUE2[[i]]$CrossName<-str_replace(rownames(CrossBLUE2[[i]]),"Crosses","")
  #   head(CrossBLUE2[[i]])
  #   colnames(CrossBLUE2[[i]])[colnames(CrossBLUE2[[i]])=="CrossBLUE2[[i]]"]<-"BLUE_DwPM_1Year"
  # }
  # 
  # CrossBLUE$BLUE_DwPM_2019<-vlookup(CrossBLUE$CrossName,dict=CrossBLUE2[[1]],result_column = "BLUE_DwPM_1Year",lookup_column = "CrossName")
  # CrossBLUE$BLUE_DwPM_2020<-vlookup(CrossBLUE$CrossName,dict=CrossBLUE2[[2]],result_column = "BLUE_DwPM_1Year",lookup_column = "CrossName")
  # 
  # head(CrossBLUE)
  # ######Check the BLUEs
  # Y<-Y2
  # commonX<-Y[Y$Year==2019,]$Crosses[which(Y[Y$Year==2019,]$Crosses%in%Y[Y$Year==2020,]$Crosses)]
  # CrossBLUE[CrossBLUE$CrossName%in%commonX,]$BLUE_DwPM_2019
  # CrossBLUE[CrossBLUE$CrossName%in%commonX,]$BLUE_DwPM_2020
  # cor(round(CrossBLUE$BLUE_DwPM_2Yrs),round(CrossBLUE$BLUE_DwPM_2019),use="complete") #1
  # ######
 
  
  
  
  
  
  
  
  
  
  
  
### Setting the CV
########
#### make more random samples
setwd(here("CVDat1920"))
folds   <- 1:10    ### If no CV, then set a -999 for "folds" in parameters

sampleCV<-matrix(nrow=nrow(Y),ncol=500)
for (n in 1:500){
sets<-rep(1:10,26)[-c(1:2)]  # nrow(Y)=nrow(CrossMerge1920)=258 lines !!!!!!!!Add 2 fake Checks
sampleCV[,n]<-sets[order(runif(nrow(Y)))]
}

save(sampleCV,file="sampleCV_0215_2021.Rdata")

cycles<-2 # !!!
ntraits<-1  # !!!
cor<-matrix(nrow=cycles,ncol=ntraits)

for (i in 1:cycles){

	tmp<-NULL

for(fold in folds){

    yNA<-y
   # print(fold)
    colCV<-i
    testing=which(sampleCV[,colCV]==fold)
    yNA[testing]=NA
    
    fm=BGLR(y=yNA,
	ETA=ETA,
	nIter=nIter,
	burnIn=burnIn,
	saveAt=paste0("CVData1920_10fold_Cycle",i),
	verbose=TRUE)
    
    predictions=data.frame(testing,Individual=gid[testing], y=y[testing], yHat=fm$yHat[testing])
	tmp<-rbind(tmp,predictions)
}  # End folds
	cor[i,]<-cor(tmp$y,tmp$yHat,use="complete")

### Save the predictions from full 10 folds
   
dir.create(paste0('10folds_Cycle',i))
setwd(paste0('10folds_Cycle',i))   ### Creating the fold_# folder
write.table(tmp,file=paste("predictions_cycle",i,".csv",sep=""),row.names=FALSE,sep=",") 
             
  # print(str(fm))
  rm(fm) 
  
  unlink("*.dat")
  
  setwd('..')

} # End cycles

cor

## 
setwd("/local/workdir/mh865/GCA_SCA/CVData1920/output/GP1_P2_P1P2/10folds_Cycle1")
cor<-matrix(nrow=cycles,ncol=ntraits)

for (i in 1:cycles){
setwd(paste0("../10folds_Cycle",i))
predictionSaved<-read.csv(paste("predictions_cycle",i,".csv",sep=""),sep=",",header=TRUE)
cor[i,]<-cor(predictionSaved$y,predictionSaved$yHat,use="complete")
}

write.csv(cor,paste0("GP1_P2_P1P2_cor_cycles_0119_2020",cycles,".csv"))
# GP1+GP2+GP1P2 Model
colMeans(cor)
# 0.489798

######## Scan the variance
varUETA3<-scan("CVData1920_10fold_Cycle50ETA_3_varU.dat")
varUETA2<-scan("CVData1920_10fold_Cycle50ETA_2_varU.dat")
varUETA1<-scan("CVData1920_10fold_Cycle50ETA_1_varU.dat")
varE<-scan("CVData1920_10fold_Cycle50varE.dat")
par(mfrow=c(2,2))
plot(varE,type='o',col=2,cex=.5)
plot(varUETA3,type='o',col=1,cex=.5)
plot(varUETA2,type='o',col=1,cex=.5)
plot(varUETA1,type='o',col=1,cex=.5)

mean(varUETA3)
mean(varUETA2)
mean(varUETA1)

### Estimate all individuals together to get the var.component
### Using only GP1, then GP2 to do CV

####### GP1+GP2 only Model
# colMeans(cor)
# 0.481274
### CV scheme parameters.R

rm(list=ls())
setwd(here("CVData1920/output"))

folds   <- 1:10    ### If no CV, then set a -999 for "folds" in parameters
nIter  <- 50000
burnIn <- 40000
phenotype.file <- '../data/CrossMerge1920.csv' 
AB <- list() 
AB[[1]] <- '../GP1/EVD.rda'       # path to pehnotype file 
AB[[2]] <- '../GP2/EVD.rda'       # path to pehnotype file 

type <- c('RKHS','RKHS','RKHS') #!!!
colENV  <- NULL  # column in phenotype file that gives the id of the environment 
colVAR  <- 2  # column in phenotype file that gives the id of the variety 
colPhen <- 5  # phenotypes column
colCV <- 7   # CV column
CV0 <- FALSE 
ESC <- FALSE 
r <- 1
set.seed(1)


#### CV scheme main code

#Load the BGLR library
library(BGLR)

Y <- read.csv(phenotype.file,sep=',', header=TRUE, stringsAsFactors=FALSE )
##rownames(Y)<-Y$GID
	
#load(files)

y   = Y[,colPhen]
gid = Y[,colVAR]

if(ESC) { y=scale(y,center=TRUE,scale=TRUE) }


## This is to define what model do you want to use, in BGLR, setting up the ETA respectively
###########
#### WHY??? This is to allow each AB a separate ETA?

nk <-length(AB)  # what is the nk here? AB has GP1, GP2, GP1P2
ETA<-list(nk)

for(i in 1:nk){
  
  if(type[i]=='BRR')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=Z,model='BRR')
    rm(Z)
  }
  
  if(type[i]=='FIXED')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=Z,model='FIXED')
    rm(Z)
  }
  
  if(type[i]=='RKHS')
  {
    load(AB[[i]])
    ETA[[i]] <- list(V=EVD$vectors,d=EVD$values,model='RKHS')
    rm(EVD)
  }
  
  if(type[i]=='BayesA')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BayesA')
    rm(X)
  }
  
  if(type[i]=='BayesB')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BayesB')
    rm(X)
  }
  
  if(type[i]=='BayesC')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BayesC')
    rm(X)
  }
  
  if(type[i]=='BL')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BL')
    rm(X)
  }  
  print(i)
}
########

### Setting the CV
########



#### make more random samples
sampleCV<-matrix(nrow=nrow(Y),ncol=500)
for (n in 1:500){
sets<-rep(1:10,26)[-c(1:4)]  # nrow(Y)=nrow(CrossMerge1920)=256 lines !!!!!!!!
sampleCV[,n]<-sets[order(runif(nrow(Y)))]
}

save(sampleCV,file="sampleCV.Rdata")


load("sampleCV.Rdata")

cycles<-50 # !!!
ntraits<-1  # !!!
cor<-matrix(nrow=cycles,ncol=ntraits)

for (i in 1:cycles){

	tmp<-NULL

for(fold in folds){

    yNA<-y
   # print(fold)


    colCV<-i
    testing=which(sampleCV[,colCV]==fold)
    
   
 #   if(CV0)
 #   {         
 #     testing <- which(gid %in% gid[testing])   #CV0 is DJ paper, random remove one at a time
 #   }  
   
  

    yNA[testing]=NA
    
    fm=BGLR(y=yNA,
	ETA=ETA,
	nIter=nIter,
	burnIn=burnIn,
	saveAt=paste0("CVData1920_10fold_Cycle",i),
	verbose=TRUE)

    	fm$y=y     ### !!!! Not clear about this step
    
    predictions=data.frame(testing,Individual=gid[testing], y=y[testing], yHat=fm$yHat[testing])
	tmp<-rbind(tmp,predictions)
}  # End folds
	cor[i,]<-cor(tmp$y,tmp$yHat,use="complete")

### Save the predictions from full 10 folds
   
dir.create(paste0('10folds_Cycle',i))
setwd(paste0('10folds_Cycle',i))   ### Creating the fold_# folder
write.table(tmp,file=paste("predictions_cycle",i,".csv",sep=""),row.names=FALSE,sep=",") 
             
  # print(str(fm))
  rm(fm) 
  
  unlink("*.dat")
  
  setwd('..')

} # End cycles

colMeans(cor)
## GP1+GP2 Model
# colMeans(cor)
# 0.481274

