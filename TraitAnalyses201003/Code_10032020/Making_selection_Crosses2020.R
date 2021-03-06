#### 2020 Selection process based on DwPM GEBVs
#{}
###########This is making selection process
### Combining all the GPs that have the GEBV
SProg_GPs<-allBLUPs[grep("UCONN-S",rownames(allBLUPs)),]  ### Grep only the SProg_GPs
allBLUPs0<-allBLUPs
allBLUPs<-allBLUPs[!rownames(allBLUPs)%in%rownames(SProg_GPs),]
####
dim(allBLUPs)
dim(allBLUPs0)
sampledSP <- which(nchar(rownames(allBLUPs)) < 11)  # The fndr SP is max of 10 characters
releaseGP <- nchar(rownames(allBLUPs))   
releaseGP <- which(10 < releaseGP & releaseGP <=15)  #FG and MG from fndr SP is max of 14 characters
progenySP <- which(nchar(rownames(allBLUPs)) > 15)   #the Cross name is larger than 15 characters


fndNames <- rownames(allBLUPs)[sampledSP]
hasDesc <- sapply(fndNames, function(n) grep(n, rownames(allBLUPs)[progenySP]))

nDesc <- sapply(hasDesc, length)
nDesc

GP_allBLUPs<-as.data.frame(rbind(allBLUPs[releaseGP,],SProg_GPs)) ## !!! Added the SProg_GPs  
dim(GP_allBLUPs)
#GP_allBLUPs<-as.data.frame(rbind(allBLUPs[releaseGP,]))
SP_allBLUPs<-as.data.frame(allBLUPs[progenySP,])
SP_allBLUPs<-SP_allBLUPs[order(-SP_allBLUPs$AshFDwPM),]

nrow(GP_allBLUPs)+length(fndNames)+nrow(SP_allBLUPs)

#write.csv(SP_allBLUPs,paste0("GOM_RepeatCrosses_",yr,"Yr.csv")) 
#RepeatCrosses<-read.csv("GOM_RepeatCrosses_BothYr.csv",sep=",",header=TRUE)
RepeatCrosses<-SP_allBLUPs

#### Adding sugar content
# This was what's saved
#save(dataRegion19_C,dataRegion20_C,dataRegionBoth_C,file="dataRegion_withChk_3_sets_PhotoScore23.rdata")

load("dataRegion_withChk_3_sets_PhotoScore23.rdata")    ### !!!!
dataRegion<-dataRegionBoth_C

Sugars<-dataRegion[,colnames(dataRegion)%in%c("Crosses","Total.Sugars","Mannuronic","Guluronic","Glucuronic","Mannitol","Fucose" ,"Glucan","Mannan","Galactan")]
Sugars<-Sugars[, apply(Sugars, 2, function(v) !all(is.na(v)))]
dim(Sugars)
RepeatCrosses2<-merge(RepeatCrosses,Sugars,by.x="row.names",by.y="Crosses",all.x=TRUE)
head(RepeatCrosses2)
colnames(RepeatCrosses2)[1]<-"Crosses"  

library(expss)
GPnames<-str_split_fixed(string=as.character(RepeatCrosses2$Crosses), "x", 2)
head(GPnames)
RepeatCrosses2$FG<-GPnames[,1]
RepeatCrosses2$MG<-GPnames[,2]
head(RepeatCrosses2)
RepeatCrosses2<-RepeatCrosses2[order(-RepeatCrosses2$AshFDwPM),]

write.csv(RepeatCrosses2,"GOM_RepeatCrosses_AddedSugar_withSGPs_866_10032020.csv")


### Generate new crosses for a trait
traitname<-c("AshFDwPM","DWpM","AshOnly")  #################!!!!!!!!!
ListGPs<-rownames(GP_allBLUPs) 
getGPsex <- function(gpName){
  sexNum <- strsplit(gpName, "-")[[1]][4]    ## split text and get the 4th element which is MG/FGs
  sex <- substring(sexNum, 1, 1)
  return(ifelse(sex %in% c("F", "M"), sex, NA))
}
gpSex <- sapply(ListGPs, getGPsex) ### apply to each element of the vector

gpSex2<-matrix(nrow=length(gpSex),ncol=2)
gpSex2[,1]<-names(gpSex)
gpSex2[,2]<-gpSex
gpSex2<-as.data.frame(gpSex2)

colnames(gpSex2)<-c("GPList","Sex")
gpSex2<-gpSex2[order(gpSex2$Sex),]

gpF<-gpSex2[gpSex2$Sex=="F",]  #data frame
gpM<-gpSex2[gpSex2$Sex=="M",]  #data frame

CrossF<-NULL
for (i in 1:nrow(gpF)){
  CrossFemal<-rep(as.character(gpF[i,1]),times=nrow(gpM))
  CrossF<-c(CrossF,CrossFemal)
}
CrossM<-rep(as.character(gpM$GPList),times=nrow(gpF))
Crosses<-as.data.frame(cbind(CrossF,CrossM))

for (col in c("CrossF","CrossM")){
  Crosses[,col]<-as.character(Crosses[,col])
}
head(Crosses)
head(allBLUPs)
for (i in 1:length(traitname)){
  Crosses[,3*i]<-vlookup(Crosses$CrossF,dict=allBLUPs0,result_column=i,lookup_column = "row.names") # Female BLUPs
  colnames(Crosses)[3*i]<-paste0("F_",colnames(allBLUPs0)[i])
  
  Crosses[,(3*i+1)]<-vlookup(Crosses$CrossM,dict=allBLUPs0,result_column=i,lookup_column = "row.names") # Male BLUPs
  colnames(Crosses)[3*i+1]<-paste0("M_",colnames(allBLUPs0)[i])
  
  Crosses[,(3*i+2)]<-rowMeans(cbind(Crosses[,3*i],Crosses[,(3*i+1)]))
  colnames(Crosses)[3*i+2]<-paste0("Predicted_",colnames(allBLUPs0)[i])
}  


Crosses$FGLoc<-str_split_fixed(string=Crosses$CrossF,"-",4)[,2]
Crosses$MGLoc<-str_split_fixed(string=Crosses$CrossM,"-",4)[,2]

SNELoc<-c("TI", "FW", "PI", "FI", "BL", "DB", "LR")
SNERegion<-cbind(rep("SNE",times=length(SNELoc)),SNELoc)

colnames(SNERegion)<-c("Region","Loc")
head(SNERegion)
str(SNERegion)
Crosses$FRegion<-vlookup(Crosses$FGLoc,dict=SNERegion,result_column ="Region",lookup_column = "Loc")
Crosses$MRegion<-vlookup(Crosses$MGLoc,dict=SNERegion,result_column ="Region",lookup_column = "Loc")
head(Crosses)
Crosses$FRegion<-ifelse(is.na(Crosses$FRegion),"GOM",Crosses$FRegion) 
Crosses$MRegion<-ifelse(is.na(Crosses$MRegion),"GOM",Crosses$MRegion) 

### Has both GOM and SNE
Crosses<-Crosses[order(-Crosses$Predicted_AshFDwPM),]
head(Crosses)
dim(Crosses)
write.csv(Crosses,paste0("GOMSNE_NewCrosses_AddedTraits_",yr,"Yr_10012020_withSGPs_866_10032020.csv"))     #### !!!!!!!!!!!!!

Crosses<-read.csv("GOMSNE_NewCrosses_AddedTraits_BothYr_10012020_withSGPs_866_10032020.csv",sep=",",header=T,row.names=1) 



##### Update checking on biomass
#### Checking on biomass !!!!!!!!!!!!!!!!!
summary<-read.csv("Active_culture_biomass_assesment_for_Mao_100120.csv",sep=",",header=T,na.strings=c("","NA"))  ###!!!!
dim(summary)  
head(summary)
tail(summary)
str(summary)
str(summary$GametophyteID)

### Onceagain confirm that the biomass is not listed  
summary[which(!as.character(summary$GametophyteID)%in%c(as.character(Crosses$CrossF),as.character(Crosses$CrossM))),]

summary$GametophyteID<-as.character(summary$GametophyteID)
str(summary)

######### Repeat Crosses
RepeatCrosses2$FGBio<-vlookup(RepeatCrosses2$FG,dict=summary,result_column="X..Crosses",lookup_column ="GametophyteID")
RepeatCrosses2$MGBio<-vlookup(RepeatCrosses2$MG,dict=summary,result_column="X..Crosses",lookup_column ="GametophyteID")

RepeatCrosses2$Manual<-ifelse(is.na(RepeatCrosses2$FGBio) |is.na(RepeatCrosses2$MGBio) | RepeatCrosses2$FGBio==0 | RepeatCrosses2$MGBio==0,"RM","Keep")
write.csv(RepeatCrosses2,"GOM_RepeatCrosses_AddedSugar_10032020_AddBiomass_withSGP_866.csv")


# RM SNE regions !!!

Crosses<-Crosses[Crosses$FRegion=="GOM"&Crosses$MRegion=="GOM",]
Crosses<-droplevels(Crosses)
dim(Crosses)
Crosses$FBio<-vlookup(Crosses$CrossF,dict=summary,result_column = "X..Crosses",lookup_column="GametophyteID")
Crosses$MBio<-vlookup(Crosses$CrossM,dict=summary,result_column = "X..Crosses",lookup_column="GametophyteID")
Crosses$Manual<-ifelse(is.na(Crosses$FBio) |is.na(Crosses$MBio) | Crosses$FBio==0 | Crosses$MBio==0,"RM","Keep")

write.csv(Crosses,paste0("GOM_NewCrosses_AddedTraits_BothYr_10032020_AddBiomass_withSGP_866.csv"))     #### !!!!!!!!!!!!!

head(Crosses)
write.csv(rownames(hMat),"rownameshMat_866.csv")
write.csv(cbind(summary$GametophyteID,summary$X..Crosses),"summaryGametophyteID.csv")
write.csv(rownames(GP_allBLUPs),"rownames_GP_allBLUPs_866.csv")
###

### Consumed biomass Repeats
# RepeatsCons<-read.csv("GOM_RepeatCrosses_AddedSugar_10012020_AddBiomass_withSGP_withNotesManual.csv",sep=",",header=T)
# RepeatsCons<-RepeatsCons[!is.na(RepeatsCons$Notes),]
# RepeatsCons<-droplevels(RepeatsCons)
#   write.csv(RepeatsConts,"")
# library("matrixStats")
#   rowCounts(RepeatsCons)

RepeatsCons<-read.csv("GOM_RepeatCrosses_UsedBiomass_Updated_10022020.csv",sep=",",header=TRUE) ##！！！
head(RepeatsCons) 


### Only calculate with the "Keep" ones   
data<-Crosses[Crosses$Manual=="Keep",]
data<-droplevels(data)
data$MBio[data$MBio==5]=10  # Males made into 10
#### UCONN-S all Biomass to be 2 to reserve some
#### UCONN-S if Biomass is 1,2 change to 0

data$CrossF<-as.character(data$CrossF)
data$CrossM<-as.character(data$CrossM)

UCONNS_F<-data$CrossF[grep("UCONN-S",data$CrossF)]
UCONNS_M<-data$CrossM[grep("UCONN-S",data$CrossM)]

data$FBio[which(data$CrossF%in%UCONNS_F)]<-ifelse(data$FBio[which(data$CrossF%in%UCONNS_F)]>2,2,0)
data$MBio[which(data$CrossM%in%UCONNS_M)]<-ifelse(data$MBio[which(data$CrossM%in%UCONNS_M)]>2,2,0)

dim(data)
head(data)
#How many were already used in Repeats
data$FBioUsed<-vlookup(data$CrossF,dict=RepeatsCons,result_column = "F_Used",lookup_column = "FG")
data$MBioUsed<-vlookup(data$CrossM,dict=RepeatsCons,result_column = "M_Used",lookup_column = "MG")

# If the FBio is used in Repeats, then biomass - used in Repeats, otherwise no change on F/MBio
data$FBioAdj<-ifelse(is.na(data$FBioUsed),data$FBio,data$FBio-data$FBioUsed) ##！！！
data$MBioAdj<-ifelse(is.na(data$MBioUsed),data$MBio,data$MBio-data$MBioUsed) ##！！！
head(data)
#num<-NULL
#nrow(data)
dim(data)

### !!!! To this point, the data is already subtracted the RepeatCrosses Biomass  
### !!!!! RM ones that are no longer any biomass left
data<-data[data$FBioAdj>0 & data$MBioAdj>0,]

data<-data[order(-data$Predicted_AshFDwPM),]  # order again by predicts ##!!! Ordering
head(data) 
dim(data)  # 18120 -> 14317 x 22
#19031 -> 14420
NewCross<-as.data.frame(matrix(NA,ncol=2))
#data.frame(Fsel=character(),Msel=character()) ##
colnames(NewCross)<-c("CrossF","CrossM")
NewCross
#NewCross[1,]<-as.vector(data[1,colnames(data)%in%c("CrossF","CrossM")])  # Start with sth

######!!! Function to pick based on biomass from top down
Pick<-function(data,NewCross){
  
  FCounts<-data[,colnames(data)%in%c("CrossF","FBioAdj")]  
  FCounts<-FCounts[!duplicated(FCounts$CrossF),]
  
  MCounts<-data[,colnames(data)%in%c("CrossM","MBioAdj")] 
  MCounts<-MCounts[!duplicated(MCounts$CrossM),]
  head(FCounts)
  head(MCounts)
  dim(MCounts)
  dim(FCounts)
  
  FCounts$FC<-NA
  MCounts$MC<-NA
  FCounts$SelF<-NA
  MCounts$SelM<-NA
  
  head(FCounts)
  head(MCounts)
  NewCross
  library(expss)
  FCounts0<-FCounts
  MCounts0<-MCounts
  
  data$SelF<-NA
  data$SelM<-NA
  data$Select<-NA
  
  for (i in 1:nrow(data)){
    
    Fname<-data$CrossF[i]
    Mname<-data$CrossM[i]
    Fname
    Mname
    FCounts<-FCounts0[which(FCounts0$CrossF==Fname),]
    FCounts$FC<-count_if(FCounts$CrossF,NewCross$CrossF)  #FC:Femalecount count how many F (used) in the NewCross list
    FCounts
    MCounts<-MCounts0[which(MCounts0$CrossM==Mname),]
    MCounts$MC<-count_if(MCounts$CrossM,NewCross$CrossM) #MC:Malecount count how many M (used) in the NewCross list
    head(FCounts)
    MCounts
    
    FCounts$SelF<-ifelse(FCounts$FC < FCounts$FBioAdj,"OK","Out") # select F
    MCounts$SelM<-ifelse(MCounts$MC < MCounts$MBioAdj,"OK","Out") # select M
    FCounts
    MCounts
    FCounts0[which(FCounts0$CrossF==Fname),]<-FCounts
    MCounts0[which(MCounts0$CrossM==Mname),]<-MCounts
    head(FCounts0)
    head(FCounts0)
    data$SelF[i]<-vlookup(data$CrossF[i],dict=FCounts0,result_column="SelF",lookup_column = "CrossF")
    data$SelM[i]<-vlookup(data$CrossM[i],dict=MCounts0,result_column="SelM",lookup_column = "CrossM")
    head(data)
    if(data$SelF[i]=="OK" & data$SelM[i]=="OK"){
      data$Select[i]<-"Select"
    }else{
      data$Select[i]<-"no"
    }
    if(data$Select[i]=="Select"){
      NewCross<-rbind(NewCross,as.vector(data[i,colnames(data)%in%c("CrossF","CrossM")]))
    }else{
      NewCross<-NewCross
    }
    head(data)
    NewCross
    print(nrow(NewCross))
  }
  return(list(data=data,NewCross=NewCross,FCounts=FCounts,MCounts=MCounts))
} 
###!! Function done

#### First all lines, find the cut off place for 160th
dataall<-Pick(data=data,NewCross = NewCross)

data1<-dataall$data
data160<-data1[data1$Select=="Select",] # select the top 160, and find where the next 170 starts on row
#data160<-data160[1:160,]
dim(data160)
data160<-data160[160,]   #### Pick top 160
data160
nrowdata<-which(rownames(data)==rownames(data160))
nrowdata
dim(data) 
data[nrowdata,]  

#### First the top 160  
data160Only<- Pick(data[1:nrowdata,],NewCross=NewCross) 
head(data160Only$data) 
head(data160Only$NewCross)
dim(data161)  
data161<- data160Only$data 
dim(data161_S)
data161_S<-data161[data161$Select=="Select",]   #### 1st Top part

data160NewCross<-data160Only$NewCross
str(data160NewCross)

dataBottom10Only<-Pick(data[(nrow(data)-10):nrow(data),],data160NewCross) #### Bottom to get just 10 picks
TopBottomNewCross<-dataBottom10Only$NewCross
str(TopBottomNewCross)
dataBottom10<-dataBottom10Only$data 
dim(dataBottom10)
head(dataBottom10)
dim(dataBottom10_S)
dataBottom10_S<-dataBottom10[dataBottom10$Select=="Select",]   #### 2nd Bottom part


dataMiddle<-data[(nrow(data161)+1):(nrow(data)-nrow(dataBottom10)),] #### 3rd middle part

dataMiddleOnly<-Pick(dataMiddle,TopBottomNewCross)

dataMiddlePart<-dataMiddleOnly$data

data161$Category<-"Top160"
dataBottom10$Category<-"Low10"
dataMiddlePart$Category<-"MiddleExtra"

dataOrder<-rbind(data161,dataBottom10,dataMiddlePart)
dim(dataOrder)
dim(data)
dataOrder$Cross<-paste0(dataOrder$CrossF,"x",dataOrder$CrossM)  
data$Cross<-paste0(data$CrossF,"x",data$CrossM)
which(!dataOrder$Cross%in%data$Cross)

write.csv(dataBottom10,"GOM_Bottom10_predicted_select.csv")
write.csv(data161,"GOM_160Crosses_from_Top.csv")
write.csv(dataMiddlePart,"GOM_Middle_predicted_extra.csv")

write.csv(data1,"GOM_AllCrosses_Ranked_basedonBiomass.csv")
write.csv(dataOrder,"GOM_Crosses_Available_Ordered_10022020.csv")


### RepeatsSum<-RepeatsCons
RepeatsFSum<-RepeatsSum[,colnames(RepeatsSum)%in%c("FG","F_Used")]
RepeatsFSum<-RepeatsFSum[!(is.na(RepeatsFSum$FG) |is.na(RepeatsFSum$F_Used)),]
RepeatsSP<-str_split_fixed(RepeatsFSum$FG,"-",4)[,1:3]
#RepeatsSP<-RepeatsSP[!RepeatsSP[,1]=="",]

RepeatsFSum$SP<-apply(RepeatsSP,1, function(x) paste(x[1:3], collapse="-"))  
RepeatsFSum
RepeatsFSum$Loc<-RepeatsSP[,2]

RepeatsMSum<-RepeatsSum[,colnames(RepeatsSum)%in%c("MG","M_Used")]
RepeatsMSum<-RepeatsMSum[!(is.na(RepeatsMSum$MG) |is.na(RepeatsMSum$M_Used)),]
RepeatsSPM<-str_split_fixed(RepeatsMSum$MG,"-",4)[,1:3]
#RepeatsSP<-RepeatsSP[!RepeatsSP[,1]=="",]

RepeatsMSum$SP<-apply(RepeatsSPM,1, function(x) paste(x[1:3], collapse="-")) 
RepeatsMSum$Loc<-RepeatsSPM[,2]
head(RepeatsFSum)
head(RepeatsMSum)

colnames(RepeatsFSum)[1]<-"GPID"
colnames(RepeatsFSum)[2]<-"Cross.num"

colnames(RepeatsMSum)[1]<-"GPID"
colnames(RepeatsMSum)[2]<-"Cross.num"

RepeatsF_Sum<-sum.num(RepeatsFSum)
RepeatsM_Sum<-sum.num(RepeatsMSum)

RepeatsF_Used<-rbind(RepeatsF_Sum$num.Loc,RepeatsF_Sum$num.SP) 
RepeatsM_Used<-rbind(RepeatsM_Sum$num.Loc,RepeatsM_Sum$num.SP)

write.csv(RepeatsF_Used,"RepeatsF_Used.csv")
write.csv(RepeatsM_Used,"RepeatsM_Used.csv")

###!! Function to sum the biomass: set needs to have "Cross.num", "Loc" and "SP"
sum.num<-function(set){
  set<-droplevels(set)
  num.Loc<-aggregate(set$Cross.num,by=list(set$Loc),FUN=sum,na.rm=TRUE)  ### Cross# per location
  num.SP<-aggregate(set$Cross.num,by=list(set$SP),FUN=sum,na.rm=TRUE)  
  return(list(num.Loc=num.Loc,num.SP=num.SP))
}    
### Summary on cultures biomass that's in stock and used by LOC
#Make the biomass
head(summary)
colnames(summary)[colnames(summary)=="X..Crosses"]<-"Cross.num"
colnames(summary)[colnames(summary)=="Location"]<-"Loc"
colnames(summary)[colnames(summary)=="Parent"]<-"SP"  

setF<-summary[summary$Sex=="F",]
setM<-summary[summary$Sex=="M",]
setF<-droplevels(setF)
setM<-droplevels(setM)
### setM all the Ms if =5, set to 10
### if UCONN, set to 2 or 0
setM$Cross.num[setM$Cross.num==5]=10

head(setM)
CultureF<-sum.num(setF)
CultureM<-sum.num(setM)
CultureF$num.Loc
CultureF$num.SP
CultureM$num.Loc
CultureM$num.SP
CultureF.num<-rbind(CultureF$num.Loc,CultureF$num.SP) 
CultureM.num<-rbind(CultureM$num.Loc,CultureM$num.SP) 

write.csv(CultureF.num,"Biomass_Culture.F.csv")
write.csv(CultureM.num,"Biomass_Culture.M.csv")


##### Summary biomass by Loc and SP for the top 170 list by LOC
set170<-dataOrder[dataOrder$Category=="Top160" | dataOrder$Category=="Low10",]
dim(set170)

###set170<-dataOrder[dataOrder$Category=="Top160",]
### Select the ones being used for making crosses
set170<-set170[set170$Select=="Select",]
set170<-droplevels(set170)
set170$Cross.num<-rep(1,nrow(set170))  # indicates the number of times each line being used
head(set170)
dim(set170)  
library(stringr)  
set170F<-set170M<-set170

set170F$Loc<-set170F$FGLoc
SPF<-str_split_fixed(set170F$CrossF,"-",4)
set170F$SP<-apply(SPF,1, function(x) paste(x[1:3], collapse="-"))
head(set170F)

set170M$Loc<-set170M$MGLoc
SPM<-str_split_fixed(set170M$CrossM,"-",4)  ##!!!
set170M$SP<-apply(SPM,1, function(x) paste(x[1:3], collapse="-"))

set170.Culture.numF<-sum.num(set170F)
set170.Culture.numM<-sum.num(set170M)

write.csv(rbind(set170.Culture.numF$num.Loc,set170.Culture.numF$num.SP),"Biomass_170_F.csv")
write.csv(rbind(set170.Culture.numM$num.Loc,set170.Culture.numM$num.SP),"Biomass_170_M.csv")  




###### Manually compare the Location distribution by Location in excel
###### What's the specific GPs being used in the 160+10 list
UsedNewCross<-dataBottom10Only$NewCross
dim(UsedNewCross)
head(UsedNewCross)

UsedFUniq<-unique(UsedNewCross$CrossF)    
UsedNewCrossF<-sapply(UsedFUniq,function(x) count_if(x,UsedNewCross$CrossF))  # count_if on each element of UsedFUniq

UsedMUniq<-unique(UsedNewCross$CrossM)    
UsedNewCrossM<-sapply(UsedMUniq,function(x) count_if(x,UsedNewCross$CrossM))  

UseNewCrossF2<-t(as.data.frame(t(UsedNewCrossF)))
UseNewCrossM2<-t(as.data.frame(t(UsedNewCrossM)))

dataOrder$FG170Used<-vlookup(dataOrder$CrossF,dict=UseNewCrossF2,result_column=1,lookup_column = "row.names")
dataOrder$MG170Used<-vlookup(dataOrder$CrossM,dict=UseNewCrossM2,result_column=1,lookup_column = "row.names")

write.csv(dataOrder,"dataOrder.csv")




####### Checking on locations summary to get diversity
set170LocF<-set170.Culture.numF$num.Loc
colnames(set170LocF)<-c("Loc","Used160_10")

set170LocM<-set170.Culture.numM$num.Loc
colnames(set170LocM)<-c("Loc","Used160_10")

setCultureF<-CultureF$num.Loc
setCultureM<-CultureM$num.Loc

colnames(setCultureF)<-c("Loc","CultureSum_F")
colnames(setCultureM)<-c("Loc","CultureSum_M")

setCultureF$Used_170<-vlookup(setCultureF$Loc,dict=set170LocF,result_column = "Used160_10",lookup_column = "Loc")
setCultureM$Used_170<-vlookup(setCultureM$Loc,dict=set170LocM,result_column = "Used160_10",lookup_column = "Loc")

SNELoc<-c("TI", "FW", "PI", "FI", "BL", "DB", "LR")
SNERegion<-cbind(rep("SNE",times=length(SNELoc)),SNELoc)

colnames(SNERegion)<-c("Region","Loc")
head(SNERegion)
str(SNERegion)

setCultureF$Region<-vlookup(as.character(setCultureF$Loc),dict=SNERegion,result_column = "Region",lookup_column = "Loc")
setCultureM$Region<-vlookup(as.character(setCultureM$Loc),dict=SNERegion,result_column = "Region",lookup_column = "Loc")

setCultureF$Region[is.na(setCultureF$Region)]="GOM"
setCultureM$Region[is.na(setCultureM$Region)]="GOM"

Culture_sum<-merge(setCultureF,setCultureM,by="Loc",all=TRUE)
Culture_sum$TotalLoc<-Culture_sum$Used_170.x+Culture_sum$Used_170.y

# number of repeats were not summed properly
head(RepeatsCons)
uniqLsF<-unique(RepeatsCons$FG)
uniqLsM<-unique(RepeatsCons$MG)
#
RepeatsCons$LocF<-str_split_fixed(string=RepeatsCons$FG,"-",4)[,2]
RepeatsCons$LocM<-str_split_fixed(string=RepeatsCons$MG,"-",4)[,2]

#RepeatsCons$TotalF<- aggregate(RepeatsCons$LocF,by=list(RepeatsCons$LocF),FUN=sum,na.rm=TRUE)
#RepeatsCons$TotalM<- aggregate(RepeatsCons$LocM,by=list(RepeatsCons$LocM),FUN=sum,na.rm=TRUE)

Culture_sum$UsedRepeatsF<-vlookup(Culture_sum$Loc,dict=RepeatsCons,result_column = "FUsed",lookup_column = "Floc")
Culture_sum$UsedRepeatsM<-vlookup(Culture_sum$Loc,dict=RepeatsCons,result_column = "Mused",lookup_column = "Mloc")

Culture_sum$LeftFBio<-Culture_sum$CultureSum_F-Culture_sum$Used_170.x-Culture_sum$UsedRepeatsF
Culture_sum$LeftMBio<-Culture_sum$CultureSum_M-Culture_sum$Used_170.y-Culture_sum$UsedRepeatsM

dim(Culture_sum)
head(Culture_sum)
write.csv(Culture_sum,"Culture_sum.csv")
#Diversity<-read.csv("Diversity_Interested_Location.csv",sep=",",header=TRUE)

#dataOrder$Select=="Select"& 
CrossLeft<-dataOrder[!(dataOrder$Category%in%c("Top160","Low10")),]
dim(CrossLeft)


##### Manually check on the locations
### Interested:   
Diversity<-CrossLeft[CrossLeft$FGLoc%in%c("LD","NC","NL","OD","CB"),]
dim(Diversity)
Diversity<-Diversity[Diversity$MGLoc%in%c("NL","CC","OD","SF"),]

write.csv(Diversity,"TO_pick_20DiversityCrosses.csv")


### Picked Diversity

Diver20<-read.csv("Picked_20DiversityCrosses.csv",sep=",",header=T)
Diver20
Diver20$Cross<-paste0(Diver20$FG,"x",Diver20$MG)
dim(dataMiddle)  
dataMiddle$Cross<-paste0(dataMiddle$CrossF,"x",dataMiddle$CrossM)

dataMiddle_ExtractDiver20<-dataMiddle[which(dataMiddle$Cross%in%Diver20$Cross),]

Diver20_bio<-read.csv("Picked_20Diversity_Biomass.csv",sep=",",header=TRUE)

#subtract the biomass on the BioAdj column in dataMiddle
dataMiddle_RMDiver20<-dataMiddle[which(!dataMiddle$Cross%in%Diver20$Cross),]
dim(dataMiddle_RMDiver20) #9254
dataMiddle_RMDiver20$FUseDiv<-vlookup(dataMiddle_RMDiver20$CrossF,dict=Diver20_bio,result_column = "FUsed",lookup_column = "FG")
head(Diver20_bio)
dataMiddle_RMDiver20$MUseDiv<-vlookup(dataMiddle_RMDiver20$CrossM,dict=Diver20_bio,result_column = "MUsed",lookup_column = "MG")
head(dataMiddle_RMDiver20)

dataMiddle_RMDiver20$FBioAdj2<-ifelse(is.na(dataMiddle_RMDiver20$FUseDiv),dataMiddle_RMDiver20$FBioAdj,dataMiddle_RMDiver20$FBioAdj-dataMiddle_RMDiver20$FUseDiv)    
dataMiddle_RMDiver20$MBioAdj2<-ifelse(is.na(dataMiddle_RMDiver20$MUseDiv),dataMiddle_RMDiver20$MBioAdj,dataMiddle_RMDiver20$MBioAdj-dataMiddle_RMDiver20$MUseDiv)       

dataMiddle_RMDiver20$FBioAdj<-dataMiddle_RMDiver20$FBioAdj2
dataMiddle_RMDiver20$MBioAdj<-dataMiddle_RMDiver20$MBioAdj2

dim(dataMiddle)
dim(dataMiddle_ExtractDiver20)
dim(dataMiddle_RMDiver20)

# Use this list added Div20 to start with the pick function 
str(TopBottomNewCross)   
Diver20Cross<-Diver20[,colnames(Diver20)%in%c("FG","MG")]
colnames(Diver20Cross)<-c("CrossF","CrossM")
Diver20Cross

TopBottomNewCross_add20<-rbind(TopBottomNewCross,Diver20Cross)    
dim(TopBottomNewCross_add20)
head(dataMiddle_RMDiver20)
dataMiddle_RMDiver20<-dataMiddle_RMDiver20[order(-dataMiddle_RMDiver20$Predicted_AshFDwPM),]
head(dataMiddle_RMDiver20)  

# Pick with updated files      
Middle_AfterDiv20<-Pick(data=dataMiddle_RMDiver20,NewCross = TopBottomNewCross_add20)   

dataMiddle_AfterDiver20<-Middle_AfterDiv20$data  ###!!
write.csv(dataMiddle_AfterDiver20,"dataMiddle_AfterDiver20.csv")
dim(dataMiddle_AfterDiver20)
head(dataMiddle_AfterDiver20)
head(data161)
dim(data161)

dataMiddle_AfterDiver20<-dataMiddle_AfterDiver20[,!colnames(dataMiddle_AfterDiver20)%in%c("Cross","FUseDiv","MUseDiv","FBioAdj2","MBioAdj2")]    
dataMiddle_AfterDiver20$Category<-"Extra"

#Make the names the same    
dataMiddle_ExtractDiver20<-dataMiddle_ExtractDiver20[,!colnames(dataMiddle_ExtractDiver20)%in%c("Cross")]
dataMiddle_ExtractDiver20$SelF<-"OK"
dataMiddle_ExtractDiver20$SelM<-"OK"
dataMiddle_ExtractDiver20$Select<-"Select"
dataMiddle_ExtractDiver20$Category<-"Diversity"

dim(dataMiddle_ExtractDiver20)
dim(data161)
colnames(data161)
colnames(dataMiddle_ExtractDiver20)
colnames(dataMiddle_AfterDiver20)
colnames(data161)
dim(dataMiddle_AfterDiver20)
head(dataMiddle_ExtractDiver20)

dataOrder2<-rbind(data161,dataBottom10,dataMiddle_ExtractDiver20,dataMiddle_AfterDiver20)

write.csv(dataOrder2,"GOM_New_Crosses_Ranked_Add20_Diversity.csv")
# for (i in 1:nrow(data)){
#     Fname<-data$CrossF[i]
#     Mname<-data$CrossM[i]
#     Fname
#     Mname
#     FCounts<-FCounts0[which(FCounts0$CrossF==Fname),]
#     FCounts$FC<-count_if(FCounts$CrossF,NewCross$CrossF)  #FC:Femalecount count how many F (used) in the NewCross list
#     FCounts
#     MCounts<-MCounts0[which(MCounts0$CrossM==Mname),]
#     MCounts$MC<-count_if(MCounts$CrossM,NewCross$CrossM) #MC:Malecount count how many M (used) in the NewCross list
#     head(FCounts)
#     MCounts
#     
#     FCounts$SelF<-ifelse(FCounts$FC < FCounts$FBioAdj,"OK","Out") # select F
#     MCounts$SelM<-ifelse(MCounts$MC < MCounts$MBioAdj,"OK","Out") # select M
#     FCounts
#     MCounts
#     FCounts0[which(FCounts0$CrossF==Fname),]<-FCounts
#     MCounts0[which(MCounts0$CrossM==Mname),]<-MCounts
#     head(FCounts0)
#     head(FCounts0)
#     data$SelF[i]<-vlookup(data$CrossF[i],dict=FCounts0,result_column="SelF",lookup_column = "CrossF")
#     data$SelM[i]<-vlookup(data$CrossM[i],dict=MCounts0,result_column="SelM",lookup_column = "CrossM")
#     head(data)
#     if(data$SelF[i]=="OK" & data$SelM[i]=="OK"){
#       data$Select[i]<-"Select"
#     }else{
#       data$Select[i]<-"no"
#     }
#     if(data$Select[i]=="Select"){
#       NewCross<-rbind(NewCross,as.vector(data[i,colnames(data)%in%c("CrossF","CrossM")]))
#     }else{
#       NewCross<-NewCross
#     }
#     head(data)
#     NewCross
#     print(nrow(NewCross))
# }


# data$Fnum<-NA
#   for (i in (1:nrow(data))){
#     tmp<-data[1:i,]
#     data$Fnum[i]<-sum(tmp$CrossF==data$CrossF[i]) # countif in excel. sum of 1 TRUE = 1
#     #data$Mnum[i]<-sum(tmp$CrossM==data$CrossM[i])
#   }
# 
# # for each levels of Flist
# # creat a subset of file
#     head(data)
#     dim(data)
# data$F_keep<-ifelse(data$Fnum<=data$FBioAdj,"F_y","F_no")   
# 
# data2<-data[data$F_keep=="F_y",]
#   dim(data2)
# data2<-droplevels(data2)  
# 
# data2$Mnum<-NA
# for (i in (1:nrow(data2))){
#       tmp<-data2[1:i,]
#       data2$Mnum[i]<-sum(tmp$CrossM==data2$CrossM[i]) # countif in excel. sum of 1 TRUE = 1
#     }
#   head(data2)
#   dim(data2)
# data$Mnum_look<-vlookup(data$CrossM,dict=data2,result_column = "Mnum",lookup_column = "CrossM")    
# 
# data$M_keep<-ifelse(data$Mnum_look<=data$MBioAdj,"M_y","M_no")
#       
# 
# # If Fnum needed < FBioAdj, then y, otherwise no  
# data$Notes<-ifelse(data$Fnum<=data$FBioAdj & data$Mnum_look<=data$MBioAdj,"y","no")
# data$CrossMade<-paste0(data$CrossF,"x",data$CrossM) # Create unique variable
#   dim(data)


#write.csv(data,file="GOM_NewCrosses_AddedTraits_BothYr_10012020_AddBiomass_withSGP_HasBiomassOnly_866.csv")

# #### Merge these info back to Crosses
# data0<-data[,colnames(data)%in%c("FBioUsed","MBioUsed", "FBioAdj", "MBioAdj", "Fnum","Mnum", "Notes","CrossMade")]
#   head(data0)
#   dim(data0)
#   tail(data0)
# 
# Crosses$CrossMade<-paste0(Crosses$CrossF,"x",Crosses$CrossM) # Create unique variable
#   dim(Crosses)
#   head(Crosses)
# Crosses<-Crosses[order(-Crosses$Predicted_AshFDwPM),]   ##!!! Ordering
# 
# Crosses_AddNotes<-merge(Crosses,data0,by="CrossMade",all.x=TRUE) ## merge back to long list
#   dim(Crosses_AddNotes)
#   head(Crosses_AddNotes)
#   tail(Crosses_AddNotes)
# Crosses_AddNotes<-Crosses_AddNotes[order(-Crosses_AddNotes$Predicted_AshFDwPM),]    ## !! Ordering
# 
# #library("xlsx")
# #write.xlsx(Crosses_AddNotes,file="GOM_NewCrosses_AddedTraits_BothYr_10012020_AddBiomass_withSGP_IndicateBiomass.xlsx",sheetName="Sheet1",append=FALSE)
# 
# write.csv(Crosses_AddNotes,file="GOM_NewCrosses_AddedTraits_BothYr_10012020_AddBiomass_withSGP_IndicateBiomass_866.csv")
# 


