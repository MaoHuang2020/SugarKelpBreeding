#Assess between Loc 

# Amatrix by location, heat-map
# Off-diagonal means

WD<-"/Users/maohuang/Desktop/Kelp/GCA_SCA/"  # run in terminal

load(paste0(WD,"OneTime1920/data/","outCovComb_dip_0116_2021.Rdata"))
#write.csv(outCovComb4_dipOrder,here("OneTime1920/data","A.csv"))
Amat<- read.csv(paste0(WD,"OneTime1920/data/","A.csv"),sep=",",header=TRUE,row.names=1)          # path to covariates file

Link<-read.csv("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/ReorderPedigree/3_plates_sequenced_names_Update.csv",sep=",",header=T)
  head(Link)
  dim(Amat)
  Amat[1:4,1:5]
locs<-levels(Link[Link$Region=="GOM",]$Location)  #Only the GOM locations
  locs<-c("CC","NC","NL","JS","OI","SF","LD","CB")

  tmp<-NULL
  tmp2<-NULL
  for (i in 1:length(locs)){
    Amatsub<-as.matrix(Amat[rownames(Amat)%in%Link[Link$Location==locs[i],]$Name_InGeno,rownames(Amat)%in%Link[Link$Location==locs[i],]$Name_InGeno]) 
    tmp<-c(tmp,mean(c(Amatsub)))
    diag(Amatsub)<-NA
    tmp2<-c(tmp2,mean(c(Amatsub),na.rm =TRUE))
  }
  names(tmp)<-locs            
tmp
tmp2
# > tmp
# [1] 0.3179973 0.2420760 0.2574035 0.2675357 0.1970107 0.4476933 0.3681782 0.3280601
# > tmp2
# [1] 0.2432661 0.1503068 0.1968030 0.1892893 0.1403637 0.3248789 0.3062829 0.2760229
# 
find_Loc<-function(v=rownames(Amat)) {
Locations<-substring(v,first=6,last=7)
order<-1:length(v)
Group<-data.frame(order,Locations)
Group<-Group[order(Group$Locations),]
return(Group)
}

CrossIndx<-grepl("x",rownames(Amat),fixed=TRUE) #Crosses vs GPs


  dim(Group)
Amat1<-Amat[find_Loc(rownames(Amat))$order,find_Loc(rownames(Amat))$order]

Amat2<-Amat[CrossIndx,CrossIndx]

Amat3<-Amat2[find_Loc(v=rownames(Amat2))$order,find_Loc(v=rownames(Amat2))$order]
heatmap(Amat)
heatmap(Amat1)
heatmap(Amat2)
par(oma=c(5,0,1,0.5))
heatmap(Amat3,Rowv=NA,Colv=NA,cexRow=0.8,cexCol=0.7)  #Order with already ordered col/rownames


