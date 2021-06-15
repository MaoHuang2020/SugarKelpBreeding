### Make EVD data
rm(list=ls())
library(here)

load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP","dataNHpi_withChk_3_sets_PhotoScore23.rdata")   ## Plot
load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP/dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi

load(here("CVData1920/data","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","malePar","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-Y1[,colKeep]
  head(Y2)
nrowchk<-nrow(Y2[Y2$crossID=="Check",]) #33
  
  #write.csv(Y2,"Y2.csv")
### parameters.R Script
#setwd('/local/workdir/mh865/GCA_SCA/CVData1920/output')

#phenotype.file <- here("CVData1920/data","CrossMerge1920.csv")          # path to pehnotype file
mm.file        <- here("CVData1920/data","A.csv")          # path to covariates file
weight.file     <- NULL    # path to weight file if NULL, no weights are used
ctr 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
std 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
weighting 	<- !is.null(weight.file) 	# Valid only if mm.file is not NULL

#P1 ###!!!!!
colIDy<-3
colnames(Y2)[colIDy]<-"P1"

#P2 ###!!!!!
colIDy <- 4  # !!!!! 3, and 4 CHANGE between GP1, and GP2
colnames(Y2)[colIDy]<-"P2"

#the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
# or NULL if ID is not used
# 4 #(basicallythe P2 column ---- 1126Notes)


### Making the EVD for P1 and P2
### Checks parameters
# cat('=1===> Reading and Checking Parameters');cat('\n')
# source('../input/parameters.R') 
# if(is.null(mm.file)){ stop('Please provide path and filename of molecular marker file') }
# if(is.null(mm.file)){ cat('Are you working with a grouping factor?. If not Please provide path and filename of covariate matrix file') }


### Reads data and perform consistency checks

cat('=2===> Reading Data ');cat('\n')
# if(!is.null(colIDy)){
  #Y <- read.table(phenotype.file,sep=',', header=TRUE, stringsAsFactors=FALSE )
#???!!!!!##   Y<-droplevels(Y2[!Y2$crossID=="Check",])  # Y: Only the experimental Crosses, Y2 contains both experimental and checks
  Y<-Y2
  IDs<-Y[,colIDy]  #All IDs need to be in the "G" matrix ???? This is all IDs, including those for checks
#   if(is.null(mm.file))Y=Y else rm(Y)
# }
  head(IDs)
G0 <- read.table(mm.file,sep=',',header=TRUE,stringsAsFactors=FALSE)
  dim(G0)  # 866x867 -> P1: 107 x 107 (106, experimental+ 1 check); P2: 121 x 121
  G0[1:4,1:5]
name <- G0[,1]
G <- G0[,-1]
colnames(G) <- name
rownames(G) <- name
G <- 2*as.matrix(G) ###### ?????? Diagonal becomes 4!
n<- dim(G)[1]

Index<-rownames(G)%in%IDs  ##P1:107, Adding index to ensure G has the list of individuals in the IDs (P1 list)
G<-G[Index,Index]  # subset the G to match up list of P1 names
  dim(G)  
# Set up the paramter for P1 and P2 seperately and the main code too
  
#SP plots+1 check  
dataNHpi_RMchk<-droplevels(Y[!Y$crossID=="Check",]) # Subset without chk levels
  dim(dataNHpi_RMchk)
phenoNamesFact<-factor(IDs,levels=rownames(G))   #Which levels should I use???? colnames are sorted alphabetically for the outCovComb 
msZ0<-model.matrix(~-1 +phenoNamesFact,data=dataNHpi_RMchk)  
    dim(msZ0)
    msZ0[1:4,1:5]
# Add checks
nrowchk<-which(is.na(phenoNamesFact))  # These are the checks that were not having any P1 nor P2
chkSpMat<-matrix(0, length(nrowchk), nrow(G))
rownames(chkSpMat)<-nrowchk
  dim(chkSpMat)   # 27 * 207 or 27 *121
#  
msZ<-rbind(msZ0,chkSpMat)
  dim(msZ)  # The rownumbers are the corresponding rows 

  
#### Should I be making the order of the Z rows the same as Y????????????
  
#if(!is.null(colIDy)){ stopifnot(all(IDs%in%rownames(G))) }
# if(!is.null(colIDy)){
  # Z<-rbind(msZ0,msZchk)  # expand the G matrix into Gc_P1 based
  #IDs<-factor(IDs,levels=rownames(G))
    str(IDs)
  #Z<-as.matrix(model.matrix(~IDs-1,data=Y))
    Z<-msZ[order(as.numeric(rownames(msZ))),]
    dim(Z)
    # order of Z is the same as that of G
    colnames(Z)<-stringr::str_remove(colnames(Z),"phenoNamesFact")
    Z[,which(colnames(Z)=="SL18-NL-7-FG3")]
    identical(colnames(Z),rownames(G))
    
  GCA<-tcrossprod(tcrossprod(Z,G),Z)  # ZGZ'
  dim(GCA)
# }

GCA1<-GCA
save(GCA,file=here("CVData1920/GP1","G.rda"))  ###!!! GP1/

GCA2<-GCA
save(GCA,file=here("CVData1920/GP2","G.rda"))  ###!!! GP2/

EVD<-eigen(GCA)
rownames(EVD$vectors)<-rownames(GCA)
## If you provide G, the RHKS will internally compute its eigen values
## provide to save time
save(EVD,file=here("CVData1920/GP1","EVD.rda")) ###!!! GP1/

save(EVD,file=here("CVData1920/GP2","EVD.rda")) ###!!! GP2/



######## Generating the EVD for P1xP2

G1.file <- here("CVData1920/GP1","G.rda")        # path to matrix file 1
G2.file <- here("CVData1920/GP2","G.rda")           # path to matrix file 2

## Parameters
# source('../input/parameters.R')

## Checks parameters
# cat('=1===> Reading and Checking Parameters');cat('\n')
# source('../input/parameters.R') 
#   if(is.null(G1.file) | is.null(G2.file)){ cat('Please provide path and filename G matrices 1 and 2') }
#   cat('The Eigen value decomposition will be perfomed with components \n')


load(G1.file)
G1 <- GCA
load(G2.file)
G2 <- GCA

GI <- G1*G2 
EVD <- eigen(GI)

save(GI,file=here("CVData1920/GP1P2","G.rda"))
save(EVD,file=here("CVData1920/GP1P2","EVD.rda"))



### CV scheme parameters.R
setwd(here("CVData1920/output"))

#phenotype.file <- here("CVData1920/data","CrossMerge1920.csv") 
#Y <- read.csv(phenotype.file,sep=',', header=TRUE, stringsAsFactors=FALSE)
  Y<-Y2
  head(Y)
  dim(Y)
#rownames(Y)<-Y$GID

colENV  <- NULL  # column in phenotype file that gives the id of the environment 
colVAR  <- which(colnames(Y)=="Crosses")  # column in phenotype file that gives the id of the variety 
colPhen <- which(colnames(Y)=="dryWgtPerM")  # phenotypes column
colCV <- ncol(Y)+1  # CV column
CV0 <- FALSE 
ESC <- FALSE 
r <- 1
set.seed(1)

nIter  <- 50000
burnIn <- 40000
#### CV scheme main code

#Load the BGLR library
library(BGLR)

	
#load(files)

y   = Y[,colPhen]
gid = Y[,colVAR]

if(ESC) { y=scale(y,center=TRUE,scale=TRUE) }


load(here("CVData1920/GP1","EVD.rda"))  # path to pehnotype file 
EVD1<-EVD
  rm(EVD)
load(here("CVData1920/GP2","EVD.rda"))       # path to pehnotype file 
EVD2<-EVD
  rm(EVD)
load(here("CVData1920/GP1P2","EVD.rda"))       # path to pehnotype file 
EVD3<-EVD
  rm(EVD)

# This is to define what model do you want to use, in BGLR, setting up the ETA respectively
###########
#### This is to allow each AB a separate ETA?
library(BGLR)

ETA<-list(list(~factor(popChk)+factor(Year),data=Y,model="FIXED"),
          list(~factor(line)+factor(block),data=Y,model="BRR"),
          list(V=EVD1$vectors,d=EVD1$values,model="RKHS"),
          list(V=EVD2$vectors,d=EVD2$values,model="RKHS"),
          list(V=EVD3$vectors,d=EVD3$values,model="RKHS")
)



### One Time prediction
########

fm=BGLR(y=y,
        ETA=ETA,
        nIter=nIter,
        burnIn=burnIn,
        saveAt="OneTime_With_Check",
        verbose=TRUE)

yHat=fm$yHat

cor(yHat,y,use="complete")

varE<-scan("OneTime_With_CheckvarE.dat")
varB<-scan("OneTime_With_CheckETA_2_varB.dat")

varU3<-scan("OneTime_With_CheckETA_3_varU.dat")
varU4<-scan("OneTime_With_CheckETA_4_varU.dat")
varU5<-scan("OneTime_With_CheckETA_5_varU.dat")

mean(varE)
mean(varU3)
mean(varU4)
mean(varU5)

plot(varU3,type='o',col=2,cex=.5)

par(mfrow=c(2,2))
plot(varE,type='o',col=2,cex=.5)
plot(varU3,type='o',col=1,cex=.5)
plot(varU4,type='o',col=1,cex=.5)
plot(varU5,type='o',col=1,cex=.5)

#
mean(varE)
0.03884398
mean(varB)
[1] 0.009886058
mean(varU3)
[1] 0.005178434
> mean(varU4)
[1] 0.001864376
> mean(varU5)
[1] 0.003799387


### Setting the CV
########



#### make more random samples
folds   <- 1:10    ### If no CV, then set a -999 for "folds" in parameters

sampleCV<-matrix(nrow=nrow(Y),ncol=500)
for (n in 1:500){
sets<-rep(1:10,26)[-c(1:2)]  # nrow(Y)=nrow(CrossMerge1920)=258 lines !!!!!!!!Add 2 fake Checks
sampleCV[,n]<-sets[order(runif(nrow(Y)))]
}

save(sampleCV,file="sampleCV_0215_2021.Rdata")

cycles<-2 # !!!
ntraits<-1  # !!!
cor<-matrix(nrow=cycles,ncol=ntraits)

for (i in 1:cycles){

	tmp<-NULL

for(fold in folds){

    yNA<-y
   # print(fold)


    colCV<-i
    testing=which(sampleCV[,colCV]==fold)
    
   
 #   if(CV0)
 #   {         
 #     testing <- which(gid %in% gid[testing])   #CV0 is DJ paper, random remove one at a time
 #   }  
   
  

    yNA[testing]=NA
    
    fm=BGLR(y=yNA,
	ETA=ETA,
	nIter=nIter,
	burnIn=burnIn,
	saveAt=paste0("CVData1920_10fold_Cycle",i),
	verbose=TRUE)

    	fm$y=y     ### !!!! Not clear about this step
    
    predictions=data.frame(testing,Individual=gid[testing], y=y[testing], yHat=fm$yHat[testing])
	tmp<-rbind(tmp,predictions)
}  # End folds
	cor[i,]<-cor(tmp$y,tmp$yHat,use="complete")

### Save the predictions from full 10 folds
   
dir.create(paste0('10folds_Cycle',i))
setwd(paste0('10folds_Cycle',i))   ### Creating the fold_# folder
write.table(tmp,file=paste("predictions_cycle",i,".csv",sep=""),row.names=FALSE,sep=",") 
             
  # print(str(fm))
  rm(fm) 
  
  unlink("*.dat")
  
  setwd('..')

} # End cycles

cor

## 
setwd("/local/workdir/mh865/GCA_SCA/CVData1920/output/GP1_P2_P1P2/10folds_Cycle1")
cor<-matrix(nrow=cycles,ncol=ntraits)

for (i in 1:cycles){
setwd(paste0("../10folds_Cycle",i))
predictionSaved<-read.csv(paste("predictions_cycle",i,".csv",sep=""),sep=",",header=TRUE)
cor[i,]<-cor(predictionSaved$y,predictionSaved$yHat,use="complete")
}

write.csv(cor,paste0("GP1_P2_P1P2_cor_cycles_0119_2020",cycles,".csv"))
# GP1+GP2+GP1P2 Model
colMeans(cor)
# 0.489798

######## Scan the variance
varUETA3<-scan("CVData1920_10fold_Cycle50ETA_3_varU.dat")
varUETA2<-scan("CVData1920_10fold_Cycle50ETA_2_varU.dat")
varUETA1<-scan("CVData1920_10fold_Cycle50ETA_1_varU.dat")
varE<-scan("CVData1920_10fold_Cycle50varE.dat")
par(mfrow=c(2,2))
plot(varE,type='o',col=2,cex=.5)
plot(varUETA3,type='o',col=1,cex=.5)
plot(varUETA2,type='o',col=1,cex=.5)
plot(varUETA1,type='o',col=1,cex=.5)

mean(varUETA3)
mean(varUETA2)
mean(varUETA1)

### Estimate all individuals together to get the var.component
### Using only GP1, then GP2 to do CV

####### GP1+GP2 only Model
# colMeans(cor)
# 0.481274
### CV scheme parameters.R

rm(list=ls())
setwd('/local/workdir/mh865/GCA_SCA/CVData1920/output')

folds   <- 1:10    ### If no CV, then set a -999 for "folds" in parameters
nIter  <- 50000
burnIn <- 40000
phenotype.file <- '../data/CrossMerge1920.csv' 
AB <- list() 
AB[[1]] <- '../GP1/EVD.rda'       # path to pehnotype file 
AB[[2]] <- '../GP2/EVD.rda'       # path to pehnotype file 

type <- c('RKHS','RKHS','RKHS') #!!!
colENV  <- NULL  # column in phenotype file that gives the id of the environment 
colVAR  <- 2  # column in phenotype file that gives the id of the variety 
colPhen <- 5  # phenotypes column
colCV <- 7   # CV column
CV0 <- FALSE 
ESC <- FALSE 
r <- 1
set.seed(1)


#### CV scheme main code

#Load the BGLR library
library(BGLR)

Y <- read.csv(phenotype.file,sep=',', header=TRUE, stringsAsFactors=FALSE )
##rownames(Y)<-Y$GID
	
#load(files)

y   = Y[,colPhen]
gid = Y[,colVAR]

if(ESC) { y=scale(y,center=TRUE,scale=TRUE) }


## This is to define what model do you want to use, in BGLR, setting up the ETA respectively
###########
#### WHY??? This is to allow each AB a separate ETA?

nk <-length(AB)  # what is the nk here? AB has GP1, GP2, GP1P2
ETA<-list(nk)

for(i in 1:nk){
  
  if(type[i]=='BRR')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=Z,model='BRR')
    rm(Z)
  }
  
  if(type[i]=='FIXED')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=Z,model='FIXED')
    rm(Z)
  }
  
  if(type[i]=='RKHS')
  {
    load(AB[[i]])
    ETA[[i]] <- list(V=EVD$vectors,d=EVD$values,model='RKHS')
    rm(EVD)
  }
  
  if(type[i]=='BayesA')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BayesA')
    rm(X)
  }
  
  if(type[i]=='BayesB')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BayesB')
    rm(X)
  }
  
  if(type[i]=='BayesC')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BayesC')
    rm(X)
  }
  
  if(type[i]=='BL')
  {
    load(AB[[i]])
    ETA[[i]] <- list(X=X,model='BL')
    rm(X)
  }  
  print(i)
}
########

### Setting the CV
########



#### make more random samples
sampleCV<-matrix(nrow=nrow(Y),ncol=500)
for (n in 1:500){
sets<-rep(1:10,26)[-c(1:4)]  # nrow(Y)=nrow(CrossMerge1920)=256 lines !!!!!!!!
sampleCV[,n]<-sets[order(runif(nrow(Y)))]
}

save(sampleCV,file="sampleCV.Rdata")


load("sampleCV.Rdata")

cycles<-50 # !!!
ntraits<-1  # !!!
cor<-matrix(nrow=cycles,ncol=ntraits)

for (i in 1:cycles){

	tmp<-NULL

for(fold in folds){

    yNA<-y
   # print(fold)


    colCV<-i
    testing=which(sampleCV[,colCV]==fold)
    
   
 #   if(CV0)
 #   {         
 #     testing <- which(gid %in% gid[testing])   #CV0 is DJ paper, random remove one at a time
 #   }  
   
  

    yNA[testing]=NA
    
    fm=BGLR(y=yNA,
	ETA=ETA,
	nIter=nIter,
	burnIn=burnIn,
	saveAt=paste0("CVData1920_10fold_Cycle",i),
	verbose=TRUE)

    	fm$y=y     ### !!!! Not clear about this step
    
    predictions=data.frame(testing,Individual=gid[testing], y=y[testing], yHat=fm$yHat[testing])
	tmp<-rbind(tmp,predictions)
}  # End folds
	cor[i,]<-cor(tmp$y,tmp$yHat,use="complete")

### Save the predictions from full 10 folds
   
dir.create(paste0('10folds_Cycle',i))
setwd(paste0('10folds_Cycle',i))   ### Creating the fold_# folder
write.table(tmp,file=paste("predictions_cycle",i,".csv",sep=""),row.names=FALSE,sep=",") 
             
  # print(str(fm))
  rm(fm) 
  
  unlink("*.dat")
  
  setwd('..')

} # End cycles

colMeans(cor)
## GP1+GP2 Model
# colMeans(cor)
# 0.481274

