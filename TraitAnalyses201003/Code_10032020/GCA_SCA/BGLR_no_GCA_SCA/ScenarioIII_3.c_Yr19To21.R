#ScenarioIII_3.c_Yr19To21.R

### Yr2019 predict Yr2020

rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/OneTime192021/"  # run in terminal
datafdr<-paste0(WD,"data/")

####1. To get the experimental factors column Info
load(paste0(datafdr,"dataNHpi_withChk_3Yrs_PhotoScore23_07122021.rdata"))   ## Plot
dataNHpi<-dataNHpi3yrs_C
dataNHpi$densityBlades<-ifelse(dataNHpi$densityBlades==0,NA,dataNHpi$densityBlades)  # densityblades as 0 then NA
dataNHpi$popChk <- ifelse(substr(dataNHpi$plotNo, 1, 1) == "Z", substr(dataNHpi$plotNo, 1, 2), "ES")  # Checks VS ES

dataNHpi$withinLoc <- ifelse(as.vector(dataNHpi$femaParLoc) == as.vector(dataNHpi$maleParLoc), 1, 0) # WithinLoc is 1
dataNHpi$development[dataNHpi$development=="#N/A"] <-NA

Y1<-dataNHpi

colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM") #,"AshFreedryWgtPerM"
Y2<-droplevels(Y1[,colKeep])
  head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])
ExpCol<-Y[,colnames(Y)%in%c("crossID","Crosses","femalPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore")]

##### !!!
pop<-"Yr19to21"
TPYr<-2019      # TPYr   
NAYr<-2021 
#####!!!

#### 2. Upload the dBLUPs (estimated WithinYr 2019 and 2020, respectively)
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata"))  #2019 and 2020 
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_n_Indiv_WithinYr2021_0721_2021.Rdata")) # 2021

### Best to merge the Yr2021 data to the other two years data
### This is just to make the Yr2021 colnames the same as the other two years

WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,!colnames(WithinYr_Both_dBLUPs)%in%c("Ash","AshFDwPM")]
  colnames(WithinYr21_dBLUPs)[!colnames(WithinYr21_dBLUPs)%in%colnames(WithinYr_Both_dBLUPs)]
  colnames(WithinYr_Both_dBLUPs)[!colnames(WithinYr_Both_dBLUPs)%in%colnames(WithinYr21_dBLUPs)]

WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,colnames(WithinYr21_dBLUPs)]
identical(colnames(WithinYr_Both_dBLUPs),colnames(WithinYr21_dBLUPs))

WithinYr_Both_dBLUPs<-rbind(WithinYr_Both_dBLUPs,WithinYr21_dBLUPs)

rownames(WithinYr_Both_dBLUPs)<-WithinYr_Both_dBLUPs$Row.names
WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,!colnames(WithinYr_Both_dBLUPs)=="Row.names"]

WithinYr_Both_dBLUPs<-droplevels(WithinYr_Both_dBLUPs)
  dim(WithinYr_Both_dBLUPs) 
  head(WithinYr_Both_dBLUPs)
identical(as.character(rownames(WithinYr_Both_dBLUPs)),as.character(WithinYr_Both_dBLUPs$plotNo.x))

traits<-colnames(WithinYr_Both_dBLUPs)[!colnames(WithinYr_Both_dBLUPs)%in%c("Row.names","Crosses.x","Crosses.y","plotNo.x","plotNo.y","Year.x","Year.y")]
#2021 does not have Ash; Its percDryWgt is still 0s

traits<-traits[!traits%in%c("Ash","AshFDwPM","percDryWgt")]
  print(traits)
save(WithinYr_Both_dBLUPs,file=paste0(datafdr,"dBLUPs_3_yrs_cal_WithinYr.Rdata"))
  
##### 

### 3. Merge the experimental factors col to drBLUP
Y1.0<-droplevels(ExpCol[ExpCol$Year==2019|ExpCol$Year==2020,])
CrossBLUEYr<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==2019|WithinYr_Both_dBLUPs$Year.x==2020,]
CrossBLUEYr<-droplevels(CrossBLUEYr[,!colnames(CrossBLUEYr)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")])

Y1<-merge(CrossBLUEYr,Y1.0,by.x="row.names",by.y="plotNo",all.x=TRUE)  ### If by "Crosses', there will be dup in 2019,2020
  dim(Y1)
  head(Y1)
Y1$plotNo<-Y1$Row.names
Y1<-Y1[,!colnames(Y1)=="Row.names"]

Y2.0<-ExpCol[ExpCol$Year==2021,]
CrossBLUEYr<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==2021,]
CrossBLUEYr<-CrossBLUEYr[,!colnames(CrossBLUEYr)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")]
Y2<-merge(Y2.0,CrossBLUEYr,all.x=TRUE,by.x="Crosses",by.y="Crosses.x")

### Y1 has one extra col
Y1<-Y1[,colnames(Y1)%in%colnames(Y2)]
  colnames(Y2)%in%colnames(Y1)
  colnames(Y1)%in%colnames(Y2)
Y1<-Y1[,colnames(Y2)]  
  identical(colnames(Y1),colnames(Y2))
Yrbind<-rbind(Y1,Y2)
Ydata<-Y   # Save out the original Y, for in case use
Y<-Yrbind  # Now Y is updated with its BLUEs from within 2019 and 2020 WithinYr
  dim(Y)

  
  
# # ##########DO it ONCE !!!!
# # ### Load the A matrix, outCovComb
# load(paste0(datafdr,"outCovComb4_Mix_Conden_0712_2021.Rdata")) ##!!!
# 
# A_950<-outCovComb4_MixOrder
# IDs<-as.character(unique(Y$Crosses))  ### RM checks??
# A<-as.matrix(A_950[IDs,IDs])
# # 
# # ###### Some plots had data, but are the same duplicated cross between years.
# # ###### Need to expand the A matrix to make it fit all those 250 crosses (Y file (6 plots contained duplicated Crosses))
# ####### This is just to expand the A and match its order to Y!!!
# 
#  Ydup<-as.character(droplevels(Y$Crosses[duplicated(Y$Crosses)]))
# 
#  Y.y<-Y[!(Y$Year==NAYr&as.character(Y$Crosses)%in%Ydup),]
# 
#  rownames(A)[!rownames(A)%in%as.character(Y.y$Crosses)]  #2020 has 2 dup crosses
#  A.DupxDup<-A[rownames(A)%in%Ydup,colnames(A)%in%Ydup]
#  A.DupxAll<-A[rownames(A)%in%Ydup,]
#  A.AllxDup<-A[,colnames(A)%in%Ydup]
#  identical(rownames(A.DupxDup),rownames(A.DupxAll))
#  identical(rownames(A.DupxAll),colnames(A.AllxDup))
# 
#  Anew.part1<-cbind(A,A.AllxDup)
#  Anew.part2<-cbind(A.DupxAll,A.DupxDup)
#  identical(rownames(A),rownames(A.AllxDup))
#  identical(colnames(A),colnames(A.DupxAll))
#  identical(colnames(Anew.part1),colnames(Anew.part2))
# 
#  Anew<-rbind(Anew.part1,Anew.part2)
#  ToPlotNo<-expss::vlookup(lookup_value=rownames(Anew),dict=Y,result_column="plotNo",lookup_column="Crosses")
#  Anew<-Anew[match(as.character(Y$Crosses),rownames(Anew)),match(as.character(Y$Crosses),colnames(Anew))]
#     identical(rownames(Anew),as.character(Y$Crosses))
#     identical(colnames(Anew),as.character(Y$Crosses))
# 
#   scheme<-"Yr1920_Yr21"
# save(A,Anew,Y,file=paste0(datafdr,"A_n_Y_",scheme,"_allES_matchedtoY.Rdata"))
# 
# ########## DO it ONCE
# #  


#### Load Anew that's already matched order to Y(250 ES indi)
#load("/local/workdir/mh865/GCA_SCA/GBLUPOnly_NoGCA_SCA/A_250Indi_allES_matchedtoY.Rdata")
  WD<-"/local/workdir/mh865/GCA_SCA/OneTime192021/"  # run in terminal
  datafdr<-paste0(WD,"data/")
load(paste0(datafdr,"dBLUPs_3_yrs_cal_WithinYr.Rdata"))
traits<-colnames(WithinYr_Both_dBLUPs)[!colnames(WithinYr_Both_dBLUPs)%in%c("Year.x","plotNo.x","Crosses.x","Year.y","plotNo.y","Crosses.y","percDryWgt","wetWgtPerM")]


##### !!!
pop<-"Yr1920_Yr21"
#TPYr<-2019      # TPYr   

pop<-"Yr19_Yr21"

pop<-"Yr20_Yr21"
#####!!!


NAYr<-2021 
load(paste0(datafdr,"A_n_Y_",pop,"_allES_matchedtoY.Rdata"))

ETA<-list(list(K=Anew,model="RKHS"))
reps<-20 # !!!
ntraits<-length(traits)  # !!!
head(Y)  

cor<-matrix(nrow=reps,ncol=length(traits))
colnames(cor)<-traits

for (j in 1:length(traits)){
  Coltrait<-traits[j]
  
  for (i in 1:reps){
    setwd(paste0(WD,pop,"_output/"))
    
    dir.create(paste0(Coltrait,"_ydrBLUP_noLocRep",i))
    savepath<-paste0(pop,"_output/",Coltrait,"_ydrBLUP_noLocRep",i,"/")
    
    
    Y$BLUE_Trait<-Y[,colnames(Y)==Coltrait]
    
    
    y<-yBLUE<-Y[,"BLUE_Trait"]   # phenotypes column  !!!
    
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
    
    yPred<-fm$ETA[[1]]$u  # equals to fm$yHat
    predict<-data.frame(testing,
                        Crosses=Y$Crosses[testing],
                        yBLUE=yBLUE[testing],
                        yPred=yPred[testing],
                        yHat=fm$yHat[testing])
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

write.csv(cor,paste0(WD,pop,"_output/","cor_",pop,"_",length(traits),"_traits_","_using_250ES_ydrBLUP_NoLoc_0712_72021.csv"))
write.csv(cor_std,paste0(WD,pop,"_output/","cor_",pop,"_",length(traits),"_traits_","_using_250ES_ydrBLUP_NoLoc_Mean_0712_2021.csv"))


