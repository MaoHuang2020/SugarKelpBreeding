### First check and see if all GPs from 2021SP plots that are having photo score of 2,3 are inside the 866 matrix.
### If not, add their fndrs, themselves
### Then update the pedigree again.
### 1.Compare the mean of Top160 vs Low10
### 2. A bivariate analysis of the Yes or No boulation

WD<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003"
datafdr<-paste0(WD,"/data/")
outputfdr<-paste0(WD,"/2021PhenotypicAnalysis/")

### I. PlotLevel
Yr2021<-read.csv(paste0(datafdr,"2021_PhenotypingDatasheet_PlotLevel.csv"),sep=",",header=TRUE)
  dim(Yr2021)
  head(Yr2021)
  
#RM PhotoScore
Yr2021_sub1<-Yr2021[Yr2021$Photo.Score%in%c(1,2,3),]   # Included PS 2 and 3

#RM other crosses categories
RMCat<-c("","1000","2000","250","4000","500","Industry Cross","Paint_replace","Spray","Spray_replace")
    dim(Yr2021_sub2)
Yr2021_sub2<-Yr2021_sub1[!Yr2021_sub1$Category%in%RMCat,]    

GPlist<-c(as.character(Yr2021_sub2$FemaleGP),as.character(Yr2021_sub2$MaleGP))
  str(GPlist)
GPlist<-unique(GPlist) # 98 unique ones

### Load the A matrix
load(paste0(datafdr,"outCovComb4_Mix_Conden_0527_2021.Rdata"))
  Pedi866<-Pediconden
  All866<-rownames(Pedi866)
NotIn866<-GPlist[!GPlist%in%All866] # 0; All the GPs are in the Pedi866

### So then add Yr 2021 data, go to the multiple Year model

load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23_UpdateAsh_0309_2021.rdata"))
load(paste0(datafdr,"dataNHim_withChk_3_sets_PhotoScore0123.rdata"))  
  
dataNHpi21_C<-Yr2021_sub2
    colnames(dataNHpiBoth_C)
    colnames(dataNHpi21_C)
    
###!!!!!!!!!!! CHANGE here, if more traits are going to be added at the plot level    
Names2021<-c("Cross.ID","FemaleGP","Female.Location","MaleGP","Male.Location","Cross.Location","Site","Painting.Density","Development.Ranking","Line","Position","Block","Photo.Score")
Names2021<-c(Names2021,"Wet.Weight.of.Plot..kg.","Length.of.Plot..m.","Wet.Weight.per.meter..kg.","X..Dry.Weight","Dry.Yield..kg.m.","Sporphyte.Density.per.Meter")
Names1920<-c("crossID","femaPar","femaParLoc","malePar","maleParLoc","CrossLoc","Region","PlantingDens","development","line","position","block","PhotoScore")
Names1920<-c(Names1920,"wetWgtPlot","lengthPlot","wetWgtPerM","percDryWgt","dryWgtPerM","densityBlades")

colnames(dataNHpi21_C)[colnames(dataNHpi21_C)%in%Names2021]<-Names1920
    dataNHpi21_C[1:3,]

## Renamed the plot level trait measurements 
dataNHpi21_C$Order_1<-substring(dataNHpi21_C$Plot..,first=2,last=5)
dataNHpi21_C$plotNo<-paste0("2021_",substring(dataNHpi21_C$Plot..,first=2,last=5))
dataNHpi21_C$Year<-2021

## Cols keep for plot level traits /// Note, not yet the "Category"
  Colkeep<-c(Names1920,"Order_1","plotNo","Year")
Merged3yrs<-rbind(dataNHpiBoth_C[,Colkeep],dataNHpi21_C[,Colkeep])

Yr19<-droplevels(Merged3yrs[Merged3yrs$Year==2019,]) #139
Yr20<-droplevels(Merged3yrs[Merged3yrs$Year==2020,]) #144
Yr21<-droplevels(Merged3yrs[Merged3yrs$Year==2021,]) #105
dataNHpi3yrs_C<-Merged3yrs
save(dataNHpi3yrs_C,Yr19,Yr20,Yr21,file=paste0(datafdr,"dataNHpi_withChk_3Yrs_PhotoScore23_07122021.rdata"))


### II. Individual measurements
Yr21_Indiv<-read.csv(paste0(datafdr,"2021_PhenotypingDatasheet_IndividualLevel.csv"),sep=",",header=TRUE)

### RM rows with NA in blade length
Yr21_Indiv<-droplevels(Yr21_Indiv[!is.na(Yr21_Indiv$Blade.Length..cm.),])

###Cols keep for individual level traits  #not the "Category" col
Names21_Indiv<-c("Blade.Length..cm.","Blade.Width.Max","Blade.Thickness..mm.","Stipe.Length..cm.","Stipe.Diameter..mm.")
Names1920_Indiv<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")

colnames(Yr21_Indiv)[colnames(Yr21_Indiv)%in%Names21_Indiv]<-Names1920_Indiv
  
### Get the "plotNo","crossID","Region","Year"
ExpCol<-c("Plot..","Crosses","CrossLoc","plotNo","Year","crossID","femaPar","femaParLoc","malePar","maleParLoc","CrossLoc","Region","PlantingDens","development","line","position","block","PhotoScore")
Exp.df<-dataNHpi21_C[,ExpCol]
dataNHim21_C<-merge(Yr21_Indiv,Exp.df,by.x="Plot..",by.y="Plot..",all.x=TRUE)

### RM rows with NA in plotNo, which means they are in the "RMCat" Category!
dataNHim21_C<-droplevels(dataNHim21_C[!is.na(dataNHim21_C$plotNo),])   ### 522 rows

### If you want extra col, add here !!!!!!!!! BUT, also need to rename the Yr21_Indiv for thos cols
ColkeepIndiv<-c("plotNo","crossID","Region","Year" ,"bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter") 

Yr1920_Indiv<-dataNHimboth_C[,ColkeepIndiv]
  str(Yr1920_Indiv)
for (col in c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")){
  Yr1920_Indiv[,col]<-as.numeric(Yr1920_Indiv[,col])
}
Yr1920_Indiv$Year<-as.numeric(Yr1920_Indiv$Year)
  
Yr2021_Indiv<-droplevels(dataNHim21_C[,ColkeepIndiv])
  str(Yr2021_Indiv)
Yr2021_Indiv$plotNo<-as.factor(Yr2021_Indiv$plotNo)

  identical(colnames(Yr1920_Indiv),colnames(Yr2021_Indiv))
dataNHim3yrs_C<-rbind(Yr1920_Indiv,Yr2021_Indiv)
  dim(dataNHim3yrs_C)
 
Yr19_Ind<-droplevels(dataNHim3yrs_C[dataNHim3yrs_C$Year==2019,]) 
  str(Yr19_Ind)
Yr20_Ind<-droplevels(dataNHim3yrs_C[dataNHim3yrs_C$Year==2020,])  
  dim(Yr20_Ind)
Yr21_Ind<-droplevels(dataNHim3yrs_C[dataNHim3yrs_C$Year==2021,])
  dim(Yr21_Ind)  
  
save(Yr19_Ind,Yr20_Ind,Yr21_Ind,dataNHim3yrs_C,file="dataNHim_withChk_3Yrs_PhotoScore0123_07122021.rdata")  
  