### Multi-trait analysis
rm(list=ls())
library(asreml)
asreml.license.activate()
#enter this code CDEA-HECC-CDAH-FIED


WD<-"/local/workdir/mh865/GCA_SCA/"
datafdr<-paste0(WD,"OneTime1920/data/")

###1. Input outCovComb4_dipOrder
load(paste0(WD,"OneTime1920/data/","outCovComb4_Mix_Conden_0420_2021.Rdata"))
All_grm<-outCovComb4_MixOrder  # 866

###2. Input dataNH

#### Plot level
load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23_UpdateAsh_0309_2021.rdata"))  ## Plot -- Updated Ash
load(paste0(datafdr,"dataNHim_withChk_3_sets_PhotoScore0123.rdata"))  ## Individual

dataNHpi<-dataNHpiBoth_C  ##!!!!!

exptlSP <- as.character(dataNHpi$popChk) == "ES"
dataNHpi$entry <- as.character(dataNHpi$popChk)
dataNHpi$entry[exptlSP] <- as.character(dataNHpi$plotNo[exptlSP])
dataNHpi$group <- as.factor(ifelse(exptlSP, 1, 0))

for (col in c( "line", "block","popChk","group","entry","Year")){ 
  dataNHpi[,col] <- factor(dataNHpi[,col])}
dataNHpi<-droplevels(dataNHpi)

### Plot Level----Add blades!!!!
traits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","Ash","AshFDwPM","densityBlades") 
dataNH<-dataNHpi   ### !!!!! Only the plot traits


#### Individual level to get their experimental factors.
#### The dataNHpi is already filtered for their phenotypic PHOTO SCORE (>2)
dataNHim<-dataNHimboth_C
dataNHim$line<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "line")
dataNHim$block<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "block")
dataNHim$Year<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Year")
dataNHim$popChk<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "popChk")
dataNHim$Crosses<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Crosses")
dataNHim$PhotoScore<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "PhotoScore")
dataNHim<-dataNHim[which(dataNHim$PhotoScore >1),]  # 3969 rows with PhotoScore >1
str(dataNHim)
dim(dataNHim)

traits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
for (col in traits){
  dataNHim[,col]<-as.numeric(dataNHim[,col]) }

####### Averaging the individual measurements to a plot level
library(dplyr)
dataNHim_avg<-aggregate(cbind(bladeLength,bladeMaxWidth,bladeThickness,stipeLength,stipeDiameter) ~ Year+line+block+Crosses, data=dataNHim, FUN=mean, na.action = na.omit)
### Plot and individual traits
dataNHpi$MergeVar<-paste(dataNHpi$Year,dataNHpi$line,dataNHpi$block,dataNHpi$Crosses,sep="_")
dataNHim_avg$MergeVar<-paste(dataNHim_avg$Year,dataNHim_avg$line,dataNHim_avg$block,dataNHim_avg$Crosses,sep="_")

dataNH1<-merge(dataNHpi,dataNHim_avg,by.x="MergeVar",by.y="MergeVar",all.x=TRUE)
dataNH1[is.na(as.character(dataNH1$Crosses.y)==as.character(dataNH1$Crosses.y)),]
colnames(dataNH1)[colnames(dataNH1)=="Year.x"]<-"Year"
colnames(dataNH1)[colnames(dataNH1)=="line.x"]<-"line"
colnames(dataNH1)[colnames(dataNH1)=="block.x"]<-"block"
colnames(dataNH1)[colnames(dataNH1)=="Crosses.x"]<-"Crosses"
colnames(dataNH1)
#######


########### Multi-trait for plot-level +Individual level traits
################################## Adding the checks into Relationship matrix

### Transform the traits value by multiply sqrt(10)
dataNH<-dataNHpi
Trait_Crosses<-as.character(dataNH$Crosses)
Trait_grm<-All_grm[rownames(All_grm)%in%Trait_Crosses,colnames(All_grm)%in%Trait_Crosses]
print(dim(dataNH))
print(sum(is.na(dataNH$Trait)))
print(dim(Trait_grm))

### Adding the checks into the Trait_grm, all 1s in diagonal,all 0s for others
data<-dataNH1 ### data has both of the Plot level and Individual level traits

data[data$popChk=="ES",]$Crosses
droplevels(data[data$popChk=="ES",])$Crosses
ChkCross<-unique(droplevels(data[!data$popChk=="ES",])$Crosses) # 33 plots of checks, 4 unique ones
Col0<-matrix(0,nrow=nrow(Trait_grm),ncol=length(ChkCross))
colnames(Col0)<-ChkCross

Trait_grm2<-cbind(Trait_grm,Col0)

Row0<-matrix(0,nrow=length(ChkCross),ncol=ncol(Trait_grm))
Chk1<-diag(x=1,nrow=length(ChkCross),ncol=length(ChkCross))
Row0_Chk1<-cbind(Row0,Chk1)
rownames(Row0_Chk1)<-ChkCross

Trait_grm3<-rbind(Trait_grm2,Row0_Chk1)

print(nrow(Trait_grm3)==length(unique(droplevels(data$Crosses))))
Trait_grm4<-Trait_grm3
Cross_link<-as.data.frame(rownames(Trait_grm3))
colnames(Cross_link)<-"Cross_chr"

Cross_link$num<-1:nrow(Trait_grm3)
ChkCross

#Cross_link$num<-ifelse(as.character(Cross_link$Cross_chr)%in%ChkCross,999,Cross_link$num0)

rownames(Trait_grm3)<-expss::vlookup(rownames(Trait_grm3),dict=Cross_link,result_column = "num",lookup_column ="Cross_chr" )
colnames(Trait_grm3)<-expss::vlookup(colnames(Trait_grm3),dict=Cross_link,result_column = "num",lookup_column ="Cross_chr" )

data$CrossesNum<-expss::vlookup(data$Crosses,dict=Cross_link,result_column = "num",lookup_column ="Cross_chr" )

data$CrossesNum<-as.factor(data$CrossesNum)
### Order the data$CrossesNum to be matching with the Trait_grm3 
data<-data[order(data$Year,data$CrossesNum),]

### Inverse the Trait_grm3
snpRelMat<-Trait_grm3
Gsnp=solve(snpRelMat+diag(1e-6,length(snpRelMat[,1]),length(snpRelMat[,1])))

# Map elements in the relationship matrix to the phenotypes

nrow(Gsnp)==sum(rownames(Gsnp)==levels(data$CrossesNum))
nrow(Gsnp)==sum(colnames(Gsnp)==levels(data$CrossesNum))


attr(Gsnp, "INVERSE")=TRUE

asreml.options(extra=30)

#Within plot level traits and Individual traits
modMT1<-asreml(cbind(percDryWgt,dryWgtPerM,AshFDwPM,densityBlades,bladeLength,bladeMaxWidth,bladeThickness,stipeLength,stipeDiameter) ~ trait+trait:line+trait:block+trait:popChk+trait:Year,
               random= ~ us(trait):vm(CrossesNum,Gsnp),
               residual = ~id(units):us(trait),
               data = data, maxiter=100, trace=TRUE)

GenVar1<-summary(modMT1)$varcomp
write.csv(GenVar1,"Plot5traits_Indivial_MultiTrait_Genetic_var_04222021.csv") #and_

colnames(GenVar1)
## Extract rows with extra parts

Gen_varcov<-GenVar1[grepl("trait:vm",rownames(GenVar1)),]
## RM names with extra parts
rownames(Gen_varcov)<-stringr::str_split_fixed(rownames(Gen_varcov),"_",2)[,2]

#Mtraits<-c("wetWgtPlot","dryWgtPerM","AshFDwPM","densityBlades")

### The varcomp is in long format, convert to wide format
Gen_varcov$var1<-stringr::str_split_fixed(rownames(Gen_varcov),":",2)[,1]
Gen_varcov$var2<-stringr::str_split_fixed(rownames(Gen_varcov),":",2)[,2]
Gen_varcov$cor<-Gen_varcov$component
Gen_varcov<-Gen_varcov[order(rownames(Gen_varcov)),]
Gen_varcov_wide<-tidyr::spread(Gen_varcov[,c("var1","var2","component")],var2,component)

write.csv(Gen_varcov_wide,"Plot6traits_Individual_Gen_varcov_wide_04222021.csv")

### Reorganize: Fill in NAs
rownames(Gen_varcov_wide)<-Gen_varcov_wide$var1
Gen_varcov_wide<-Gen_varcov_wide[,-1]

Gen_varcov_imp<-matrix(nrow=nrow(Gen_varcov_wide),ncol=ncol(Gen_varcov_wide))
rownames(Gen_varcov_imp)<-rownames(Gen_varcov_wide)
colnames(Gen_varcov_imp)<-colnames(Gen_varcov_wide)

for (i in 1:nrow(Gen_varcov_wide)){
  for (j in 1:ncol(Gen_varcov_wide)){
    Gen_varcov_imp[i,j]<-ifelse(is.na(Gen_varcov_wide[i,j]),Gen_varcov_wide[j,i],Gen_varcov_wide[i,j])
  }
}

isSymmetric(Gen_varcov_imp)
Gen.cor<-cov2cor(Gen_varcov_imp)
orders<-c("percDryWgt","dryWgtPerM","AshFDwPM","densityBlades","bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")

Gen.cor<-Gen.cor[match(orders,rownames(Gen.cor)),match(orders,colnames(Gen.cor))]

write.csv(Gen.cor,"Plot4traits_Individual_5traits_Multi_trait_Genetic_cor_04262021.csv") #
