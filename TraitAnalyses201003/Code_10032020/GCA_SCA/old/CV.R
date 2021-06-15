### Make EVD data
### Parameters
setwd<-'/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTimePrediction/output'
phenotype.file <- '../data/Y.csv'          # path to pehnotype file
mm.file        <- '../data/A.csv'          # path to covariates file
weight.file     <- NULL    # path to weight file if NULL no weights are used


ctr 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
std 		<- TRUE # TRUE/FALSE		# Valid only if mm.file is not NULL
weighting 	<- !is.null(weight.file) 	# Valid only if mm.file is not NULL



colIDy <- 3  # the P1 column. column in phenotype file that gives the IDs that link observations to covariates or grouping factor
# or NULL if ID is not used
# 4 #(basicallythe P2 column ---- 1126Notes)

## Checks parameters
cat('=1===> Reading and Checking Parameters');cat('\n')
source('../input/parameters.R') 
#if(is.null(mm.file)){ stop('Please provide path and filename of molecular marker file') }
if(is.null(mm.file)){ cat('Are you working with a grouping factor?. If not Please provide path and filename of covariate matrix file') }

### Reads data and perform consistency checks
cat('=2===> Reading Data ');cat('\n')
if(!is.null(colIDy)){
  Y <- read.table(phenotype.file,sep=',', header=TRUE, stringsAsFactors=FALSE )
  IDs<-Y[,colIDy] # All IDs need to be in the G matrix
  if(is.null(mm.file))Y=Y else rm(Y)
}

G <- read.table(mm.file,sep=',',header=TRUE,stringsAsFactors=FALSE)
name <- G[,1]
G <- G[,-1]
colnames(G) <- name
rownames(G) <- name
G <- 2*as.matrix(G) 
n<- dim(G)[1]

Index<-rownames(G)%in%IDs  ## Adding index to ensure G has the list of individuals in the IDs (P1 list)
G<-G[Index,Index]  # subset the G to match up list of P1 names

# Set up the paramter for P1 and P2 seperately and the main code too

if(!is.null(colIDy)){ stopifnot(all(IDs%in%rownames(G))) }

if(!is.null(colIDy)){
  IDs<-factor(IDs,levels=rownames(G))
  Z<-as.matrix(model.matrix(~IDs-1))  # expand the G matrix into Gc_P1 based
  G<-tcrossprod(tcrossprod(Z,G),Z)  # ZGZ'
}


save(G,file='../GP1/G.rda')

EVD<-eigen(G)
rownames(EVD$vectors)<-rownames(G)
## If you provide G, the RHKS will internally compute its eigen values
## provide to save time

save(EVD,file='../GP2/EVD.rda')







### CV scheme
folds   <- 1:5
nIter  <- 500
burnIn <- 400
phenotype.file <- '../../../../../../root.Data/Z2_1.csv' 
AB <- list() 
AB[[1]] <- '../../../../../../G/GcI.P1/output/EVD.rda'       # path to pehnotype file 
AB[[2]] <- '../../../../../../G/GcI.P2/output/EVD.rda'       # path to pehnotype file 
AB[[3]] <- '../../../../../../I/GcP1P2/output/EVD.rda'       # path to pehnotype file
AB[[4]] <- '../../../../../../G/E/output/EVD.rda'       # path to pehnotype file

type <- c('RKHS') #!!!
colENV  <- NULL  # column in phenotype file that gives the id of the environment 
colVAR  <- 3  # column in phenotype file that gives the id of the variety 
colPhen <- 17
colCV <- 14
CV0 <- FALSE 
ESC <- FALSE 
r <- 1
set.seed(1)
#Load the BGLR library
library(BGLR)

## Parameters
source('../input/parameters.R')

Y <- read.csv(phenotype.file,sep=',', header=TRUE, stringsAsFactors=FALSE )


#load(files)

y   = Y[,colPhen]
gid = Y[,colVAR]


if(ESC) { y=scale(y,center=TRUE,scale=TRUE) }





nk <-length(AB)
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


for(fold in folds)
{
  
  yNA<-y
  print(fold)
  
  if(fold != -999)
  {
    
    dir.create(paste('fold_',fold,sep=''))
    setwd(paste('fold_',fold,sep=''))
    
    testing=which(Y[,colCV]==fold)
    
    if(CV0)
    {         
      testing <- which(gid %in% gid[testing])   
    }  
    
    
    yNA=y
    yNA[testing]=NA
    
    
    fm=BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose=TRUE)
    fm$y=y
    
    predictions=data.frame(testing,Individual=gid[testing], y=y[testing], yHat=fm$yHat[testing])
    
    write.table(predictions,file=paste("predictions_",fold,".csv",sep=""),row.names=FALSE,sep=",") # Change to a unique name?
    
  }else{
    
    dir.create('fullData')
    setwd('fullData')
    
    fm=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose=TRUE)
    save(fm,file='fm_full.RData')
    
  }          
  
  print(str(fm))
  rm(fm)
  
  
  unlink("*.dat")
  
  setwd('..')
  
}

