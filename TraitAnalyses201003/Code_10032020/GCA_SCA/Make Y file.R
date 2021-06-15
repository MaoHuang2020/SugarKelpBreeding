# Format the mega file with all possible combinations
# Get the list of GPs individuals, from that 866 hMat matrix
# These include the GPs being used in making crosses; not used for crosses but were genotyped and had biomass
### I think later the hMat matrix will be removed with its founders????
rm(list=ls())
setwd("/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction_tst")
load("hMat_PedNH_CCmat_fndrMrkData_Both_PhotoScore23_WithSGP_866.rdata")
  ls()
  dim(hMat)
  hMat[1:4,1:5]
  write.csv(hMat,"/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/data/A.csv")

all_individuals<-rownames(hMat)
all_FGs<-all_individuals[grep("FG",all_individuals)]
all_FGs<-all_FGs[-grep("x",all_FGs)]
all_MGs<-all_individuals[grep("MG",all_individuals)] # keep only with MG
all_MGs<-all_MGs[-grep("x",all_MGs)] # RM crosses with "x"

ListGPs<-c(all_FGs,all_MGs)
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

gpF<-droplevels(gpSex2[gpSex2$Sex=="F",] ) #data frame
gpM<-droplevels(gpSex2[gpSex2$Sex=="M",] ) #data frame

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

#### Now, here are all the FG and MG combinations
#### Need to bring in the phenotypic data for the ones that are tested in the field in year 2019 and 2020

Crosses$Cross<-paste0(Crosses$CrossF,"x",Crosses$CrossM)
colnames(Crosses)[1:3]<-c("P1","P2","GID")

# Load the phenotypic data
load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHpi_withChk_3_sets_PhotoScore23.rdata")   ## Plot
load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi
  ls()
  head(dataNHpiBoth_C)

### Get the LSmeans on the DwPM first. 
### So that the checks were taken into consideration and being removed
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

fitAug <- lm(Trait ~ Year+ line*Year+ block*Year+ Crosses , data=dataNHpi) # Crosses:group
       summary(fitAug)
       
##### ANOVA 
     anova(fitAug)
       #Response: Trait
       #             Df    SumSq MeanSq   Fvalue    Pr(>F)    
       # Year         1  2.6616 2.66155 257.1232 6.311e-08 ***
       #  line         4  0.1102 0.02754   2.6606 0.1025722    
       # block       10  3.2172 0.32172  31.0803 9.104e-06 ***
       #  Crosses    246 30.2028 0.12278  11.8609 0.0001887 ***
       #  Year:line    3  0.0326 0.01087   1.0501 0.4168238    
       # Year:block   7  0.0536 0.00766   0.7398 0.6466632    
       # Residuals    9  0.0932 0.01035 
       
       
       colnames(summary(fitAug)$coef)
CrossBLUE<-summary(fitAug)$coef[grep("x",rownames(summary(fitAug)$coef)),"Estimate"]
  str(CrossBLUE)
CrossBLUE<-as.data.frame(CrossBLUE)
library(stringr)
CrossBLUE$CrossName<-str_replace(rownames(CrossBLUE),"Crosses","")
  head(CrossBLUE)
library(expss)
 
##### Plot out data   
df <-dataNHpi[,colnames(dataNHpi)%in%c("crossID","Year","plotNo","Region","plotNo_Numeric","femaPar","femaParLoc","malePar","maleParLoc","Crosses","CrossLoc","line","position","block","PhotoScore","Trait")]
  #colnames(df)[1] <- "Name"  # crossID
  df<-df[order(df$Crosses,df$Year),] 
  df$Year<-as.factor(df$Year)
  
  its_fndr<-function(GPstring){
    fndr<-strsplit(GPstring,split="-", fixed=T) 
    its_fndr<-sapply(fndr,function(vec) paste(vec[1:3], collapse="-"))
  }
  
  df$FGfndr<-its_fndr(as.character(df$femaPar))
  df$MGfndr<-its_fndr(as.character(df$malePar))
  
  df$fndrPair<-paste0(df$FGfndr,"x",df$MGfndr)
  df$LocPair<-paste0(df$femaParLoc,"x",df$maleParLoc) 
  
  plot(df$Year,df$Trait, pch=16)
  pdf("Trait_2020vs2019_Boxplot.pdf", height=8, width=8)
  
  plot(df$Trait,df$Year, pch=16)
  
  library(ggplot2) 
  
  plot<-ggplot(data=df,aes(Trait,Year))+
    geom_point(aes(color=Year))+ 
    geom_line(aes(group=as.factor(Crosses)))
  print(plot)


 
  
  
### Add the TraitBLUE back to the Crosses file
Crosses$TraitBLUE<-vlookup(Crosses$GID,dict=CrossBLUE,result_column = "CrossBLUE",lookup_column="CrossName")
  head(Crosses)
  which(!is.na(CrossMerge$TraitBLUE)) 

### This is to merge the Cross pairs to raw phenotypic data w/ photo score 2 and 3 !!!!!!!!Could CHANGE!!
CrossMerge<-merge(Crosses,dataNHpiBoth_C,by.x="GID",by.y="Crosses",all.x=TRUE)
    head(CrossMerge)
    tail(CrossMerge)
    which(!is.na(CrossMerge$dryWgtPerM))

CrossMerge$Year[is.na(CrossMerge$TraitBLUE)&is.na(CrossMerge$Year)]<-2021  #### Make the new combinations as 2021
  head(CrossMerge)
CrossMerge<-CrossMerge[order(CrossMerge$Year),]

colkeep<-c("GID","P1","P2","TraitBLUE","Year")
write.csv(CrossMerge[,colkeep],"/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/data/CrossMergeDwPMTrial1.csv")
  

# Set up the CV scheme
CrossMerge<-read.csv("/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/data/CrossMergeDwPMTrial1.csv",sep=",",header=TRUE)

# Subset only the 2019 and 2020 samples
CrossMerge1920<-droplevels(CrossMerge[CrossMerge$Year==2019 |CrossMerge$Year==2020,])


# random sampling
sets<-rep(1:10,26)[-c(1:4)]  # nrow(CrossMerge1920)=256 lines
sample<-sets
CrossMerge1920$CV<-sets[order(runif(nrow(CrossMerge1920)))]  # column 7
write.csv(CrossMerge1920,"/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/data/CrossMerge1920.csv")

  
#fitAug <- lmer(dryWgtPerM ~ popChk + line*block + (1|entry:group), data=dataNHpi)
# print(aov <- anova(fitAug)) 
# anova(fitAug)

#2019_C ANOVA # Photot score only 2 and 3
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq   Mean Sq NumDF   DenDF F value  Pr(>F)  
# popChk     0.01043 0.0052155     2   2.000  0.4595 0.68518  
# line       0.06803 0.0226775     3 100.342  1.9978 0.11914  
# block      0.16775 0.0209688     8  97.945  1.8473 0.07726  # this was not sig when using all data, instead of photoscore=2,3 data.
# line:block 0.32593 0.0135803    24  95.000  1.1964 0.26555 

# Both_C ANOVA
# npar    Sum Sq   Mean Sq F value
# popChk        4 0.0067382 0.0016845  0.1749
# line          1 0.0072646 0.0072646  0.7543
# block         1 0.0006827 0.0006827  0.0709
# line:block    1 0.0003968 0.0003968  0.0412


  