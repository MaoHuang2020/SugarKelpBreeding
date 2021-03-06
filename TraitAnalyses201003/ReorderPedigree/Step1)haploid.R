rm(list=ls())
###### This is to work the  the haploid level of matrices
#### 1. Founders marker data, make haploid level, one sample 2 rows: 0->0,0 1->0,1 2->1,1
#### Cal relationship matrix for the fnders
#### 2. Cal hapccMatrix for the pedigree
#### 3. GPs relationship matrix, 0->0,2->1
#### 4. Use CovComb() to combine this Founder relationship amat, and the GP sequenced amat and the PedNH ccmatrix
#### 5. condense to diploid level
#######

#### function to make diploid marker data to be haploid
#!!!!!! This is for haploid
######### Diploid split into 2 rows (onerow at a time) format
ConvertOneRow<-function(fndrMrkDataImp){
onerow<-fndrMrkDataImp
onerow[ onerow< -0.5]<--1
onerow[onerow>= -0.5 & onerow<= 0.5]<-0
onerow[onerow> 0.5]<-1

### Two row format for haploid
geno3<-onerow
geno3[geno3==-1]=0

geno4<-onerow
geno4[geno4==1]=1
geno4[geno4==0]=1
geno4[geno4==-1]=0
  
geno4<-as.matrix(geno4)
geno3<-as.matrix(geno3)
genoC<-matrix(nrow=nrow(geno4)*2,ncol=ncol(geno4))
colnames(genoC)<-colnames(geno4)

for (i in 1:nrow(geno4)){
  genoC[((2*i)-1),]=as.vector(geno3[i,])
  genoC[(2*i),]=geno4[i,]
  genoC<-as.data.frame(genoC)
  rownames(genoC)[(2*i)-1]<-rownames(geno3)[i]
  rownames(genoC)[(2*i)]<-paste0(rownames(geno4),"_2")[i]
  }

return(list(OneRowG=genoC))
}
#####

InputFile<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/ReorderPedigree"

#### Input files:
load(paste0(InputFile,"/CovList_3_As_0116_2021.Rdata"))  #####!!!!! Now fndrs has 58 individuals
  # To get the list of fndrs from diploid level: "fndrsA"
fndrsA_diploid<-rownames(fndrsA)

load(paste0(InputFile,"/FarmCPU_GAPIT.Rdata"))  
  # To get the fndrs marker data: "geno2"
geno2<-geno[,-1]
rownames(geno2)<-geno$taxa
rownames(geno2) <- paste0(substring(rownames(geno2), first=1, last=2), "18", substring(rownames(geno2), first=3))
geno2[geno2==0]=-1
geno2[geno2==1]=0
geno2[geno2==2]=1

biphasicPedNH<-read.csv(here("TraitAnalyses201003/ReorderPedigree","Ped_in_Order_866_Individuals_Fndr_New_Order_0116_2021.csv"),sep=",",header=TRUE,row.names=1)
  # To get the pedigree at diploid level": "biphasicPedNH"
source(here("TraitAnalyses201003/Making_haploid_CovComb","calcCCmatrixBiphasic.R")) 
  # To get the function calculating haploid level ccMatrix

load(here("../2020_2019_Phenotypic_Data/SugarKelpBreeding_NoGitPush/GenotypicData_for_SugarKelpBreeding","GPs_mainGenome_NA0.8_P1P2P3.Rdata"))
  # To get the GPs the raw SNPs data: "GPSNP" is the raw SNPs"

###IF to RUN in Terminal
# load("/local/workdir/mh865/SNPCalling/Saccharina_latissima/bamAll/bam/mpileup/Filtering/GPs_mainGenome_NA0.8.Rdata")


#### 1. fnder marker data
#######Did not do this step
#geno2<-geno2[rownames(geno2)%in%CrossedSP,]  # 93 fnders made crosses, only 56 were genotyped
#so some fnders were not in the pedigree file but were included to estimate aMat, they were removed out later
#though this should not be an issue?????
#######

############### make this list only the same as the diploid level 47 individuals used

geno2<-geno2[rownames(geno2)%in%fndrsA_diploid,]

## Impute fndr Mrk at diploid level, convert to haploid format
fndrMrkData<-geno2
library(rrBLUP)
mrkRelMat <- A.mat(fndrMrkData, impute.method="EM", return.imputed=T,shrink=TRUE) ## Add shrink per Deniz
fndrMrkDataImp <- mrkRelMat$imputed

fndOneRow<-ConvertOneRow(fndrMrkDataImp)$OneRowG

# Calculate GRM as (W %*% W^T) / sum(locusDosageVariance)
# To calculate p, the function expects dosage coding to be 0, 1, 2 for diploid
# and 0, 1 for haploid
  calcGenomicRelationshipMatrix <- function(locusMat, ploidy=2){
    if (!any(ploidy == 1:2)) stop("Ploidy must be 1 or 2")
    freq <- colMeans(locusMat) / ploidy
    locusMat <- scale(locusMat, center=T, scale=F)
    return(tcrossprod(locusMat) / sum(ploidy*freq*(1-freq)))
  }
  
fndrA<-calcGenomicRelationshipMatrix(fndOneRow,ploidy=1)  ###!!! 1. fnder marker Rel matrix, haploid

##Another way of estimating the relationship Matrix. BUT fndrA and fndrA2 are not the same????
#M<-fndOneRow  
#fndrA2<-cov(t(M))/mean(diag(cov(t(M))))  #assuming M is (nxp) GenotypesxMarkers allele frequencies matrix (p>>n)

  
#### 2. Read in the diploid Pedigree file, ccMatrix change rows names    

biphasicPedNH<-as.matrix(biphasicPedNH)  

HaploidCCcal <- calcCCmatrixHaploid(biphasicPedNH)

biphasichapccMat<-HaploidCCcal$hccMat   ###!!! 2. CCmatrix/aMat at haploid level for pedigree
hapOnePointer<-HaploidCCcal$hapOnePointer  ###!!! For condensing
dipTwoPointers<-HaploidCCcal$dipTwoPointers  ###!!! For condensing

### This is to rename the rownames out of the biphasichapccMat
pedinames<-NULL
for (i in 1:nrow(biphasicPedNH)){
  if(!is.na(biphasicPedNH[i,2]) & !is.na(biphasicPedNH[i,3])){
    # this is diploid, both parents have numbers
    pedinames<-c(pedinames,rownames(biphasicPedNH)[i],paste0(rownames(biphasicPedNH)[i],"_2"))
  }else{
    # this is GP, one parent is missing
    pedinames<-c(pedinames,rownames(biphasicPedNH)[i])
  }
}

rownames(biphasichapccMat)<-colnames(biphasichapccMat)<-pedinames


#### 3. GPs data  # did not redo this part on 0116_2021

#### Need to redo the formating of the raw marker data again  
  GPs0<-GPSNP   # 
  rownames(GPs0)<-paste0(GPs0$chrom,"_",GPs0$pos)  
  map<-GPs0[,1:2]
  GPs<-GPs0[,-c(1:2)]
  
  ######{{}} Did this in the Terminal to estimate the RelationshipMatrix2 
  #### Correct the names of GPs
  
  library(stringr)
  names<-colnames(GPs)
  names<-str_replace_all(names,"SA-","SA18-")
  names<-str_replace_all(names,"SL-","SL18-")
  names<-str_replace_all(names,"FG-","FG")
  names<-str_replace_all(names,"MG-","MG")
  
  colnames(GPs)<-names
  
  ## RM checks, the mistake one (RM SL-PI-1-FG-1,SL-CT1-FG-3)
  IndiRM<-c("SL18-PI-1-FG1","3","SL18-CT1-FG3","SL18-CT1-MG2","SL18-CT1-MG3","SL18-OI-15-Female","SL18-ME-1-FG1","SL18-ME-1-MG1","SL18-ME-1-MG2")
  GPs1<-GPs[,!colnames(GPs)%in%IndiRM]
  GPs<-GPs1
  
  #Impute GP marker matrix with A.mat()
  tSNP<-t(GPs)
  tSNP[tSNP==2]=1  # Change to haploid level dosage
    
#mrkMat needs to be m xn dimension
  
  imputeToMean <- function(mrkMat){
    imputeVec <- function(mrkVec){
      meanScore <- mean(mrkVec, na.rm=T)
      mrkVec[is.na(mrkVec)] <- meanScore
      return(mrkVec)
    }
    imputedMrkMat <- apply(mrkMat, 2, imputeVec)
    return(imputedMrkMat)
  }
  
GPsMrkImp<-imputeToMean(tSNP)  
  dim(GPsMrkImp)  #269 x 909749
  
#Change GPs marker to 0, 1 format
  
onerow<-GPsMrkImp 
  onerow[onerow> 0 & onerow< 0.5]=0
  onerow[onerow>= 0.5 & onerow< 1]=1
  onerow[onerow> 1]=1  
  
### Different ways of estimating GPsAs
GPsA<-calcGenomicRelationshipMatrix(GPsMrkImp,ploidy=1)

GPsA2<-calcGenomicRelationshipMatrix(onerow,ploidy=1) ##Use This !!!## This function assumes data to be 0 and 1 for haploid

cultevo::mantel.test(dist(GPsA),dist(GPsA2),trials=99) # 0.969
  
M<-GPsMrkImp
GPsA3<-cov(t(M))/mean(diag(cov(t(M))))
mantel.test(dist(GPsA),dist(GPsA3),trials=99) # 0.984

M2<-onerow
GPsA4<-cov(t(M2))/mean(diag(cov(t(M))))
mantel.test(dist(GPsA2),dist(GPsA4),trials=99) # 0.9

### Choose the onerow, calcGenomicRelationshipMatrix(), so GPsA2

### This is the GPsA2

### 4. ###

source("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Code_10032020/is.square.matrix.R")
source("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Code_10032020/is.positive.definite.R")
source("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Code_10032020/is.symmetric.matrix.R")

GPsA<-GPsA2   #269
diag(fndrA) <- diag(fndrA) + 1e-5
  is.positive.definite(fndrA)   #116
diag(GPsA) <- diag(GPsA) + 1e-5    #269
  is.positive.definite(GPsA)
diag(biphasichapccMat) <- diag(biphasichapccMat) + 1e-5   #1215
  is.positive.definite(biphasichapccMat)

# ### RM the extra fndr individuals that's not in the pedigree list
# RMextraFndr<-which(!rownames(fndrA)%in%rownames(biphasichapccMat))
# fndrA2<-fndrA[-RMextraFndr,-RMextraFndr]
#   dim(fndrA2)
#   fndrA2[29:30,29:30]
  
### check and see Is GPsA2 individuals all in the biphasichapccMat?
  which(!rownames(GPsA)%in%rownames(biphasichapccMat))  
GPsA<-GPsA[rownames(GPsA)%in%rownames(aMat),colnames(GPsA)%in%rownames(aMat)]  # 269 x 269
  which(!rownames(fndrsA)%in%rownames(biphasichapccMat))
is.positive.definite(fndrA)  

save(fndrA,GPsA,biphasichapccMat,file="Newly_saved_3_As_for_CovComb_0116_2021.Rdata")
save(hapOnePointer,dipTwoPointers,file="Pointers_from_biphasichapccMat_0116_2021.Rdata")


##{{Done in terminal}}
library(CovCombR)
load(here("TraitAnalyses201003/ReorderPedigree","Newly_saved_3_As_for_CovComb_0116_2021.Rdata"))
load(here("TraitAnalyses201003/ReorderPedigree","Pointers_from_biphasichapccMat_0116_2021.Rdata"))

CovList<-NULL
CovList[[1]]<-fndrA ## fndrsA

CovList[[2]]<-GPsA  
CovList[[3]]<- biphasichapccMat  

weights<-c(2,2,1)
outCovComb4<-CovComb(Klist=CovList,nu=2000,w=weights,Kinit=biphasichapccMat) 
####### HAS TO BE RE-ORDERED !!!
outCovComb4_HapOrder<-outCovComb4[match(rownames(biphasichapccMat),rownames(outCovComb4)),match(colnames(biphasichapccMat),colnames(outCovComb4))]

source("calcCCmatrixBiphasic.R")
outCovComb4_Hapconden<-condenseHapCCmat(outCovComb4_HapOrder,hapOnePointer,dipTwoPointers)

save(outCovComb4_HapOrder,outCovComb4_Hapconden,file="outCovComb4_hap_Conden_0116_2021.Rdata")


###{{}}


##### FOR JL: Ignore these part for now
#####
### compare dip and hap outCovComb
rm(list=ls())
# load(here("TraitAnalyses201003/ReorderPedigree","outCovComb_dip_0116_2021.Rdata"))
load(here("TraitAnalyses201003/data","outCovComb_dip_0116_2021.Rdata"))
  outCovComb4_dipOrder[1:4,1:4]
  outCovComb4_dipOrder[103:106,103:106]
hMat_dip<-outCovComb4_dipOrder

# load(here("TraitAnalyses201003/Making_haploid_CovComb","outCovComb4_hap_Conden_0116_2021.Rdata"))
# hMat_hap<-outCovComb4_HapOrder<-outCovComb4_Hapconden

load(here("TraitAnalyses201003/data","hMat_hap_0116_2021.Rdata")) # this outCovComb4_Hapconden=hMat_hap

#load(paste0("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/hMat_PedNH_CCmat_fndrMrkData_Both_PhotoScore23_withSGP_866.rdata"))

rownames(hMat_hap)<-rownames(hMat_dip)  # order of hMat_hap = biphasichapccMat=biphasicPedNH=hMat_dip
colnames(hMat_hap)<-colnames(hMat_dip)

save(hMat_hap,outCovComb4_Hapconden,file="/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Making_haploid_CovComb/hMat_hap_0116_2021.Rdata")

library(cultevo)
mantel.test(dist(as.matrix(hMat_dip)),dist(as.matrix(hMat_hap)),trials=99) # r =0.959  #47fndrs_in_fndrsA:0.923???
####

hMat_hap2<-hMat_hap
diag(hMat_hap2)<-NA

hMat_dip2<-hMat_dip
diag(hMat_dip2)<-NA

cor.test(c(hMat_dip),c(hMat_hap))  #0.977  #47fndrs_in_fndrsA: 0.964
cor(c(hMat_hap2),c(hMat_dip2),use="complete")  #0.972   #47fndrs_in_fndrsA: 0.955
cor.test(c(hMat_hap[lower.tri(hMat_hap)]),c(hMat_dip[lower.tri(hMat_dip)])) #0.972  #47fndrs_in_fndrsA: 0.955

cor(c(hMat_hap2[1:56,1:56]),c(hMat_dip2[1:56,1:56]),use="complete") #0.996  #47fndrs_in_fndrsA: 0.617
cor(c(hMat_hap[1:56,1:56]),c(hMat_dip[1:56,1:56]),use="complete") #0.997  #47fndrs_in_fndrsA: 0.943


# with diagonal  !!!
low_hapM<-hMat_hap
low_hapM[upper.tri(low_hapM)]<-NA
 
low_dipM<-hMat_dip
low_dipM[upper.tri(low_dipM)]<-NA


# without diagonal  !!!
low_hapM2<-low_hapM
diag(low_hapM2)<-NA

low_dipM2<-low_dipM
diag(low_dipM2)<-NA


corM<-function(n1,n2,matrx1,matrx2){
  if(n1==n2){
    cor(c(matrx1[n1,]),c(matrx2[n1,]),use="complete")
    }else{
      cor(c(matrx1[n1:n2,]),c(matrx2[n1:n2,]),use="complete")
      }
}

# there is a total of 11 different categories
#  rownames(cormatrix)<-c("fndrGP_geno","fndrGeno_GPNo","fndrNo_GPGeno","fndrGP_No","GP1","GP2","GP3","GP5","SP1","SP2","SP_GP","All")
# GP1, GP3 ---GPs genotyped
# GP2, GP4, GP5 ---GPs not genotyped
# SP1, SP2 --- SPs

n1<-c(1,57,60,94,105,264,332,443,544,788,789,1)
n2<-c(56,59,93,104,263,331,442,543,787,788,866,866)
  length(n1)==length(n2)
cormatrix<-matrix(nrow=length(n1),ncol=3)

# each category has its own staring and ending row number in the pedigree matrix
for (i in 1:nrow(cormatrix)){
  cormatrix[i,1]<-n2[i]-n1[i]+1
  cormatrix[i,2]<-corM(n1[i],n2[i],low_hapM,low_dipM)
  cormatrix[i,3]<-corM(n1[i],n2[i],low_hapM2,low_dipM2)
}

cormatrix
  colnames(cormatrix)<-c("samplesize","WithDiagonal","noDiagonal")
  rownames(cormatrix)<-c("fndrGP_geno","fndrGeno_GPNo","fndrNo_GPGeno","fndrGP_No","GP1","GP2","GP3","GP5","SP1","SP2","SP_GP","All")

  GPs_genotyped_rows<-c(105:263,332:442)
  GPs_NOTgenotyped_rows<-c(264:331,443:543,789:866)
  SPs_rows<-c(544:788)
  SPs_GPs_rows<-c(789:866)
 
cor2<-function(rows,data1=low_hapM,data2=low_dipM){
  cor2<-cor(c(data1[rows,rows]),c(data2[rows,rows]),use="complete")
  print(dim(data1[rows,rows]))
  print(cor2)
} 
cor.GPs_genotyped_rows<-cor2(rows=GPs_genotyped_rows)  
cor.GPs_NOTgenotyped_rows<-cor2(rows=GPs_NOTgenotyped_rows) # Included the SPs_GPs_row  
cor.SPs_rows<-cor2(SPs_rows)  
cor_SPs_GPs_rows<-cor2(SPs_GPs_rows)

  
BLUP_dip<-read.csv(here("allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_0116_2021_dip.csv"),sep=",",header=T)
BLUP_hap<-read.csv(here("allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_0116_2021_hap.csv"),sep=",",header=T)

for (i in 1:(ncol(BLUP_dip)-1)){
  print(cor(BLUP_dip[,(i+1)],BLUP_hap[,(i+1)]))
}

# #[1] 0.9921583
# [1] 0.9884357
# [1] 0.9908348
# [1] 0.9858998
# [1] 0.9899603

#58 in fndrsA 
  # samplesize WithDiagonal noDiagonal
  # fndrGP_geno           56    0.9988876  0.9984221
  # fndrGeno_GPNo          3    0.9980901  0.9920499
  # fndrNo_GPGeno         34    0.9932212  0.9866290
  # fndrGP_No             11    1.0000000 -0.1306208   ## Off-diag all 0s
  # GP1                  159    0.9811114  0.9568204
  # GP2                   68    0.9952201  0.9907836
  # GP3                  111    0.9768215  0.9667586
  # GP5                  101    0.9932301  0.9909059
  # SP1                  244    0.9779482  0.9751890
  # SP2                    1    0.9233533  0.9163822
  # SP_GP                 78    0.9651917  0.9543390
  # All                  866    0.9804502  0.9724217
  
  
#47 in fndrsA individuals
  # WithDiagonal noDiagonal
  # fndr1    0.9683710  0.6169740
  # fndr2    0.9139768  0.5249622
  # fndr3    0.9831948  0.9348435
  # fndr5    1.0000000 -0.1117056
  # GP1      0.9748336  0.9378237
  # GP2      0.9610352  0.9222910
  # GP3      0.9731881  0.9588681
  # GP5      0.9794020  0.9685276
  # SP1      0.9651113  0.9601893
  # SP2      0.9308626  0.9243578
  # SP_GP    0.9661335  0.9550919
  # All      0.9698842  0.9552952

# ##hMat_hap, Went to Step3_Cleaned_V2 for h2 estimation
