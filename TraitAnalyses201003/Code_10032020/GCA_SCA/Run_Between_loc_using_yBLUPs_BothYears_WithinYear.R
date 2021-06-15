#Run_between_loc_yBLUE_as_y
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


#CB=Farmed skinny kelp from giant staircase 
#CC: cape Cod Canal; JS:; LD:Lubec Dock; NC:Newcastle; NL: Nubble Light; OI: Orr's Island; SF:
#(NB=narrow blade=GS=Giant staircase. TI and CB are farmed)
# JS= Fort Stark, New Hampshire (SNE?)
# SF:Sullivan Falls
##From XW PCA: NL,JS,NC,IS is grouped, close to CC

locs<-c("CC","NC","NL","JS","OI","SF","LD","CB")


if(WithinYr==TRUE){
 
  pop<-"WithinYr_BetweenLoc"
  
  ########## WithinYr Data
  datafdr<-paste0(WD,"OneTime1920/data/")
  load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata")) 
  WithinYr_Both_dBLUPs<-droplevels(WithinYr_Both_dBLUPs)
  identical(as.character(WithinYr_Both_dBLUPs$Row.names),as.character(WithinYr_Both_dBLUPs$plotNo.y))
  
  rownames(WithinYr_Both_dBLUPs)<-WithinYr_Both_dBLUPs$Row.names
  WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,-1]
  head(WithinYr_Both_dBLUPs)
  
  traits<-colnames(WithinYr_Both_dBLUPs)[!colnames(WithinYr_Both_dBLUPs)%in%c("Row.names","Crosses.x","Crosses.y","plotNo.x","plotNo.y","Year.x","Year.y")]
  print(traits)
  
  ### Make NaNs to NAs
  for (t in 1:length(traits)){
    trait<-traits[t]
    WithinYr_Both_dBLUPs[,trait][is.nan(WithinYr_Both_dBLUPs[,trait])]<-NA
  }
  

reps<-20 # !!!
ntraits<-length(traits)  # !!!
  head(Y)  

cor<-matrix(nrow=reps,ncol=length(locs))
colnames(cor)<-locs

### Start the BGLR
for (j in 1:length(traits)){
  Coltrait<-traits[j]
 
  #!!!!!!!!!! USE drBLUP estimated WithinYr 2019 and 2020 !
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
  Y<-Yrbind
  
  y<-Y[,"BLUE_Trait"]  # phenotypes column  !!!
  yBLUE<-Y[,"BLUE_Trait"]
  
for (i in 3:reps){
  setwd(paste0(WD,"OneTime1920/BetweenLoc_output/")) #where to create Rep# folder
  dir.create(paste0(Coltrait,"_ydBLUPs_Rep",i))
  WDloc<-paste0(WD,"OneTime1920/BetweenLoc_output/",Coltrait,"_ydBLUPs_Rep",i,"/")  # the path following WD
  
  setwd(WDloc)   #where to create loc# folder # This is where the r is saved
  
  r<-NULL
  for(loc in locs) {
    
    dir.create(paste0("loc",loc))
    savepath<-paste0("OneTime1920/BetweenLoc_output/",Coltrait,"_ydBLUPs_Rep",i,"/loc",loc,"/")
    
    yNA<-y 
    testing<-which(Y$femaParLoc==loc)  # 156 plots
    yNA[testing]<-NA  # Other locations to predict this testing one
    
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
    r<-c(r,cor(predict$yBLUE,predict$yPred,use="complete"))
    write.table(predict,file=paste0(WD,savepath,"predictions_rep",i,".csv"),row.names=FALSE,sep=",") 
  }
  
  names(r)<-locs
  write.csv(r, paste0("r_1Loc_PP_Rep",i,".csv"))
   }
  print(r)
  }
}else if(WithinYear==FALSE){
  
  pop<-"BothYears_BetweenLoc"
  
  ## Both Year data
  datafdr<-paste0(WD,"OneTime1920/data/")
  load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_overTwoYears_AddBD.Rdata")) ##!!!
  rownames(Both_dBLUPs)<-Both_dBLUPs$Row.names
  Both_dBLUPs<-Both_dBLUPs[,-1]
  traits<-colnames(Both_dBLUPs)
    print(traits)
  CrossBLUE<-Both_dBLUPs ##!!!
  
  ### Make NaNs to NAs
  for (t in 1:length(traits)){
    trait<-traits[t]
    CrossBLUE[,trait][is.nan(CrossBLUE[,trait])]<-NA
  }
  

  reps<-20 # !!!
  ntraits<-length(traits)  # !!!
  head(Y)  
  
  cor<-matrix(nrow=reps,ncol=length(locs))
  colnames(cor)<-locs
  
  for (j in 1:length(traits)){
    Coltrait<-traits[j]
    
    Y$BLUE_Trait<-expss::vlookup(Y$Crosses,dict=CrossBLUE,result_column = paste0(Coltrait),lookup_column = "row.names") ### This is the DwPM
      head(Y)
    
    y<-Y[,"BLUE_Trait"]  # phenotypes column  !!!
    yBLUE<-Y[,"BLUE_Trait"]
    
    for (i in 1:reps){
      setwd(paste0(WD,"OneTime1920/BetweenLoc_output/ydBLUPs_BothYearData/")) #where to create Rep# folder
      dir.create(paste0(Coltrait,"_ydBLUPs_Rep",i))
      WDloc<-paste0(WD,"OneTime1920/BetweenLoc_output/ydBLUPs_BothYearData/",Coltrait,"_ydBLUPs_Rep",i,"/")  # the path following WD
      
      setwd(WDloc)   #where to create loc# folder # This is where the r is saved
      
      r<-NULL
      for(loc in locs) {
        
        dir.create(paste0("loc",loc))
        savepath<-paste0("OneTime1920/BetweenLoc_output/ydBLUPs_BothYearData/",Coltrait,"_ydBLUPs_Rep",i,"/loc",loc,"/")
        
        yNA<-y 
        testing<-which(Y$femaParLoc==loc)  # 156 plots
        yNA[testing]<-NA  # Other locations to predict this testing one
        
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
        r<-c(r,cor(predict$yBLUE,predict$yPred,use="complete"))
        write.table(predict,file=paste0(WD,savepath,"predictions_rep",i,".csv"),row.names=FALSE,sep=",") 
      }
      
      names(r)<-locs
      write.csv(r, paste0("r_1Loc_PP_Rep",i,".csv"))
    }
    print(r)
  }
  
}


#### Within Year Data
WDStart<-paste0(WD,"OneTime1920/BetweenLoc_output/")

#### Both Year data 
WDStart<-paste0(WD,"OneTime1920/BetweenLoc_output/ydBLUPs_BothYearData/")

rAll_traits_Loc<-NULL
for (j in 1:length(traits)){
  Coltrait<-traits[j]
rAll_Loc<-NULL
for (i in 1:reps){
  WDloc<-paste0(WDStart,Coltrait,"_ydBLUPs_Rep",i,"/") 
  
  rAll_Loc[[i]]<-read.csv(paste0(WDloc,"r_1Loc_PP_Rep",i,".csv"),row.names=1)
      }  
  r_allLoc<-rowMeans(do.call(cbind.data.frame, rAll_Loc))
  rAll_traits_Loc<-rbind(rAll_traits_Loc,r_allLoc)
}
rownames(rAll_traits_Loc)<-traits

write.csv(rAll_traits_Loc,paste0(WD,"OneTime1920/BetweenLoc_output/","cor_BetweenLoc_20Reps_ydBLUPs_",pop,"_",length(traits),"_traits.csv"))



