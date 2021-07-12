### ScenarioIII Setting up the A (expand, match to Y crosses)
rm(list=ls())
WD<-"/local/workdir/mh865/GCA_SCA/GBLUPOnly_NoGCA_SCA/"  # run in terminal
datafdr<-"/local/workdir/mh865/GCA_SCA/OneTime1920/data/"

####1. To get the experimental factors column Info
load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23.rdata"))   ## Plot
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])
ExpCol<-Y[,colnames(Y)%in%c("crossID","Crosses","femalPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore")]


#### 2. Upload the dBLUPs (estimated WithinYr 2019 and 2020, respectively)
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata")) 
WithinYr_Both_dBLUPs<-droplevels(WithinYr_Both_dBLUPs)

identical(as.character(WithinYr_Both_dBLUPs$Row.names),as.character(WithinYr_Both_dBLUPs$plotNo.y))

rownames(WithinYr_Both_dBLUPs)<-WithinYr_Both_dBLUPs$Row.names
WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,-1]
head(WithinYr_Both_dBLUPs)

traits<-colnames(WithinYr_Both_dBLUPs)[!colnames(WithinYr_Both_dBLUPs)%in%c("Row.names","Crosses.x","Crosses.y","plotNo.x","plotNo.y","Year.x","Year.y")]
print(traits)
##### 

#### 3. Merge the experimental factors col to drBLUP
Y1.0<-ExpCol[ExpCol$Year==2019,]
CrossBLUEYr<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==2019,]
CrossBLUEYr<-CrossBLUEYr[,!colnames(CrossBLUEYr)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")]
Y1<-merge(Y1.0,CrossBLUEYr,all.x=TRUE,by.x="Crosses",by.y="Crosses.x")

Y2.0<-ExpCol[ExpCol$Year==2020,]
CrossBLUEYr<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==2020,]
CrossBLUEYr<-CrossBLUEYr[,!colnames(CrossBLUEYr)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")]
Y2<-merge(Y2.0,CrossBLUEYr,all.x=TRUE,by.x="Crosses",by.y="Crosses.x")

Yrbind<-rbind(Y1,Y2)
Ydata<-Y   # Save out the original Y, for in case use
Y<-Yrbind  # Now Y is updated with its BLUEs from within 2019 and 2020 WithinYr


##########DO it ONCE
### Load the A matrix, outCovComb
load(paste0("/local/workdir/mh865/outCovComb/outCovComb4_Mix_Conden_0527_2021.Rdata")) ##!!!

A_866<-outCovComb4_MixOrder
IDs<-as.character(unique(Y$Crosses))  ### RM checks??
A<-as.matrix(A_866[IDs,IDs])

###### Some plots had data, but are the same duplicated cross between years.
###### Need to expand the A matrix to make it fit all those 250 crosses (Y file (6 plots contained duplicated Crosses))
####### This is just to expand the A and match its order to Y!!!

Ydup<-as.character(droplevels(Y$Crosses[duplicated(Y$Crosses)]))

Y.y<-Y[!(Y$Year==NAYr&as.character(Y$Crosses)%in%Ydup),]

rownames(A)[!rownames(A)%in%as.character(Y.y$Crosses)]  #2020 has 2 dup crosses
A.DupxDup<-A[rownames(A)%in%Ydup,colnames(A)%in%Ydup]
A.DupxAll<-A[rownames(A)%in%Ydup,]
A.AllxDup<-A[,colnames(A)%in%Ydup]
identical(rownames(A.DupxDup),rownames(A.DupxAll))
identical(rownames(A.DupxAll),colnames(A.AllxDup))

Anew.part1<-cbind(A,A.AllxDup)
Anew.part2<-cbind(A.DupxAll,A.DupxDup)
identical(rownames(A),rownames(A.AllxDup))
identical(colnames(A),colnames(A.DupxAll))
identical(colnames(Anew.part1),colnames(Anew.part2))

Anew<-rbind(Anew.part1,Anew.part2)
ToPlotNo<-expss:vlookup(lookup_value=rownames(Anew),dict=Y,results_column="plotNo",lookup_column="Crosses")
Anew<-Anew[match(as.character(Y$Crosses),rownames(Anew)),match(as.character(Y$Crosses),colnames(Anew))]
identical(rownames(Anew),as.character(Y$Crosses))
identical(colnames(Anew),as.character(Y$Crosses))

save(A,Anew,file="/local/workdir/mh865/GCA_SCA/GBLUPOnly_NoGCA_SCA/A_250Indi_allES_matchedtoY.Rdata")
########## DO it ONCE
