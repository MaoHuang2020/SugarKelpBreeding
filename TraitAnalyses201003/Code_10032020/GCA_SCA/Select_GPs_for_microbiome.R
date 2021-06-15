### Select GPs for microbiome

### Subset the list of 245 GPs genotyped
### Subset out that had biomass

### Merge to OC file (generated 0204_2021 on DwPM)

### Re-generate BLUPs for all 866 individuals, because Ash is updated --- to get blade density BVs
### How to deal with the density? If it is 0, then remove that plot??
### Or should I average them?

### Re-generate OC for blade density too
### Need to organize and sort a file with all BLUEs, then all BLUPs, and then all OCs

# 1. select top and bottom DwPM OC
# 2. Rank blade density raw or BV? Maybe BV
rm(list=ls())
wd<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding"
datafdr<-paste0(wd,"/TraitAnalyses201003/data/")

load(paste0(datafdr,"outCovComb_dip_0116_2021.Rdata"))
ls()

GPsequenced<-read.csv(paste0(datafdr,"3_plates_sequenced_names_Update.csv"),sep=",",header=TRUE)
  dim(GPsequenced)
head(GPsequenced)

## GPs genotyped, used in pedigree
GPs270<-droplevels(GPsequenced$Name_InGeno[GPsequenced$Name_InGeno%in%rownames(outCovComb4_dipOrder)])
  str(GPs270)

GPs8<-droplevels(GPsequenced$Name_InGeno[!GPsequenced$Name_InGeno%in%rownames(outCovComb4_dipOrder)])
  GPs8    # checks GPs

head(GPsequenced)   

## Upload the OC_DwPM
OC_DwPM<-read.csv(paste0(datafdr,"Candidates_TwoTraits_oc_0204_2021.csv"),sep=",",header = TRUE)
  head(OC_DwPM)
  dim(OC_DwPM)

## GPsequenced with OC_DwPM
GPsequenced$OC_DwPM<-expss::vlookup(GPsequenced$Name_InGeno,dict=OC_DwPM,"oc_using_DWpM","Indiv")  

  head(GPsequenced)
  
## Upload the BD data
allBLUPs<-read.csv(paste0(datafdr,"allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_0315_2021_dip.csv"),sep=",",header=TRUE)  
  head(allBLUPs)  
plot(allBLUPs$BD,allBLUPs$DWpM)  

GPsequenced$BV_DwPM<-expss::vlookup(GPsequenced$Name_InGeno,dict=allBLUPs,"DWpM","X")
GPsequenced$BV_BD<-expss::vlookup(GPsequenced$Name_InGeno,dict=allBLUPs,"BD","X")

GPsequenced<-GPsequenced[order(-GPsequenced$OC_DwPM),]

# ## Merge the biomass
# Biomass<-read.csv(paste0(datafdr,"UCONN Culture update 3_3_2021.csv"),sep=",",header=TRUE)
#   head(Biomass)                  
  
### Grep the BLUPs for crosses
SP_BLUPs<-droplevels(allBLUPs[grepl("x",allBLUPs$X),])

plot(SP_BLUPs$DWpM,SP_BLUPs$BD)
abline(v=0.05,col="red")
abline(h=0,col="red")

abline(v=0,col="blue")
abline(h=10,col="blue")

Subset1<-SP_BLUPs[SP_BLUPs$BD<0 &SP_BLUPs$DWpM>0.05,]
Subset1<-Subset1[order(-Subset1$DWpM,Subset1$BD),]

Subset2<-SP_BLUPs[SP_BLUPs$DWpM<0 & SP_BLUPs$BD>10,]
  dim(Subset1)
Subset2<-Subset2[order(-Subset2$BD,Subset2$DWpM),]

Subset1$Category<-"High_DWpM_Low_BD_BVs"
Subset2$Category<-"High_BD_Low_DWpM_BVs"
  dim(Subset2)

Subsets<-rbind(Subset1,Subset2)
Subsets$FG<-stringr::str_split_fixed(string=as.character(Subsets$X),"x",2)[,1]
Subsets$MG<-stringr::str_split_fixed(string=as.character(Subsets$X),"x",2)[,2]

write.csv(Subsets,paste0(datafdr,"Subsets_DWpM_contraversial_BD.csv"))
SP_BLUPs$DWpM_to_BD<-abs(SP_BLUPs$DWpM)/abs(SP_BLUPs$BD)
  SP_BLUPs<-SP_BLUPs[order(-SP_BLUPs$DWpM_to_BD),]
GPsequenced$SP_Category_used_FG<-expss::vlookup(GPsequenced$Name_InGeno,dict=Subsets,lookup_column = "FG",result_column="Category")
GPsequenced$SP_Category_used_MG<-expss::vlookup(GPsequenced$Name_InGeno,dict=Subsets,lookup_column = "MG",result_column="Category")

plot(GPsequenced$BV_DwPM,GPsequenced$BV_BD)
abline(v=0.16,col="red")
abline(h=20,col="red")
abline(v=0.12,col="blue")
abline(h=65,col="blue")

GPsequenced$GP_highDWpM_lowBD<-ifelse(GPsequenced$BV_DwPM>=0.16 & GPsequenced$BV_BD<20,"GP_highDWpM_lowBD",NA)
GPsequenced$GP_lowDWpM_highBD<-ifelse(GPsequenced$BV_DwPM<=0.12 & GPsequenced$BV_BD>65,"GP_lowDWpM_highBD",NA)

write.csv(GPsequenced,paste0(datafdr,"Select_for_microbiome_GPsequenced_add_OC_DwPM_BV_BladeDensity.csv"))

### Manual selection

GPselect_microbiome<-read.csv(paste0(datafdr,"Select_for_microbiome_GPsequenced_add_OC_DwPM_BV_BladeDensity.csv"),sep=",",header=TRUE)

GPselect_microbiome$Select[GPselect_microbiome$Select==""]<-NA
### Load GP marker data
load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/SugarKelpBreeding_NoGitPush/GenotypicData_for_SugarKelpBreeding/GPsAmat_NA0.8_P1P2P3_09282020.Rdata")
source(paste0(wd,"/TraitAnalyses201003/Code_10032020/PCA_fcn_codes.R"))
ls()
GPSNPs<-imputedFromAmat

GPSNP<-GPSNPs[rownames(GPSNPs)%in%as.character(GPselect_microbiome[!is.na(GPselect_microbiome$Select),]$Name_InGeno),]
  dim(GPSNP)
  
pop<-"GP_select_microbiome"
PCA<-PCA_fcn(GPSNP,pop = pop)
Prop<-PCA$Portion
xlim=c(-450,600)
ylim=c(-450,600)

PCA2<-PCA_lim_fcn(xlim=xlim,ylim=ylim)


PCA.link<-read.csv(paste0(getwd(),"/",pop,"_PCA_scores.csv"),sep=",",header=TRUE,row.names=1)
  head(PCA.link)
  
PCA.link$Entry<-rownames(PCA.link) #!!
PCA.link$Program<-expss::vlookup(rownames(PCA.link),dict=droplevels(GPselect_microbiome),
                                 lookup_column = "Name_InGeno",
                                 result_column = "Select")

Levels<-levels(PCA.link$Program)
  Levels
Levelscol<-c("green","black","royalblue","skyblue","mediumorchid","mediumpurple","red")
  Levelscol
color_fnc<-function (x){
  if (x=="Top20_OC_DwPM"){
    return("red")
  }else if(x=="Bottom20_OC_DwPM"){
    return("green")
  }else if(x=="GOM_checks"){
    return("black")
  }else if (x=="GP_highDWpM_lowBD")
  { return ("royalblue")} 
  else if(x=="GP_lowDWpM_highBD"){
    return("skyblue")
  }else if(x=="ProgSP_highDWpM_lowBD"){
      return("mediumorchid")
  }else if(x=="ProgSP_lowDWpM_highBD"){
    return("mediumpurple")
  }
}


#PCA.link<-droplevels(PCA.link[,c("PC1","PC2","Entry","Program")])
  str(PCA.link)
#PCA.link$Program<-as.character(PCA.link$Program)

dev.set(3) 	
tiff(file=paste(pop,"_PCA_x,y,lim_with_color.tiff",sep=""),width=1200,height=800,units="px",pointsize=12,res=150)
par(cex=0.5)

plot(PCA.link$PC1,PCA.link$PC2,pch=1,main=paste("PCA for",pop,sep=" "),xlim=xlim,ylim=ylim,xlab=paste("PC1 (",Prop[1],"%)",sep=""),ylab=paste("PC2 (",Prop[2],"%)",sep=""),col=sapply(PCA.link$Program, color_fnc))

legend(370,615,legend=Levels,col=Levelscol,pch=1,cex=0.55)

dev.off()

