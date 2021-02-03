# Further evaluate Ash content

###########  Further evaluating Ash!!!
# AshOnly: With_58fndrs_in_OutCovComb: Both, 0.110; 2019, 0.409; 2020, 0.234  
Ashdata<-droplevels(dataNHpiBoth_C[which(!is.na(dataNHpiBoth_C$Ash)),])
# plot out the Ash data for each of the same crosses, by year
#library(tidyverse)
df <- Ashdata[,colnames(Ashdata)%in%c("crossID","Year","plotNo","Region","plotNo_Numeric","femaPar","femaParLoc","malePar","maleParLoc","Crosses","CrossLoc","line","position","block","PhotoScore","Ash")]
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

plot(df$Year,df$Ash, pch=16)
pdf("Ash_2020vs2019_Boxplot.pdf", height=8, width=8)

plot(df$Ash,df$Year, pch=16)

library(ggplot2) 

plot<-ggplot(data=df,aes(Ash,Year))+
  geom_point(aes(color=Year))+
  geom_line(aes(group=as.factor(LocPair)))
print(plot)

dim(df)
head(df)

data<-df
# subset genomic relationship matrix
Ash_grm<-outCovComb4_dipOrder[rownames(outCovComb4_dipOrder)%in%df$Crosses,colnames(outCovComb4_dipOrder)%in%df$Crosses]
save(data,Ash_grm,file="Trait_Ash.Rdata")

#### Run in terminal
#R
library(asreml)
asreml.license.activate()
#enter this code CDEA-HECC-CDAH-FIED

mod1 <- asreml(Ash ~ Year,
               random= ~ us(Year):vm(Crosses,Ash_grm),
               data = data, maxiter=100, trace=TRUE)
# genetic cor: -0.45
mod3 <- asreml(Ash ~ Year+line+block,
               random= ~ us(Year):vm(Crosses,Ash_grm),
               data = data, maxiter=100, trace=TRUE)
# genetic cor: -0.49

dataNHpi$Year<-as.factor(dataNHpi$Year)
dataNHpi$popChk<-as.factor(dataNHpi$popChk)
dataNHpi$line<-as.factor(dataNHpi$line)
dataNHpi$block<-as.factor(dataNHpi$block)

All_grm<-outCovComb4_dipOrder
Trait_grm<-outCovComb4_dipOrder[rownames(outCovComb4_dipOrder)%in%as.character(dataNHpi$Crosses),colnames(outCovComb4_dipOrder)%in%as.character(dataNHpi$Crosses)]

save(dataNHpi,All_grm,Trait_grm,file="Asreml_OtherTraits.Rdata")

mod4<-asreml(AshFreedryWgtPerM~Year+line+block+popChk,
             random= ~ us(Year):vm(Crosses,Trait_grm),
             data = dataNHpi, maxiter=100, trace=TRUE)

summary(mod4)$varcomp
# r for AshFDwPM: 0.4455884 using Trait_grm
# r for AshFDwPM: 0.4390511 using All_grm
# warnings(): Crosses has levels in the data that are missing in All_grm ## Checks??  


summary(mod1)$varcomp
summary(mod3)$varcomp

snpRelMat<-Ash_grm
Gsnp=solve(snpRelMat+diag(.0005,length(snpRelMat[,1]),length(snpRelMat[,1])))
matrixattr(Gsnp,"INVERSE")=TRUE  ## did not work??

mod2<-asreml(Ash ~ Year,
             random= ~ us(Year):vm(Crosses,Gsnp),
             data = data, maxiter=100, trace=TRUE)
summary.asreml(mod2)$varcomp


#genetic cor: 0.24
#genetic cor= covarianceA_B/sqrt(varianceA*varianceB)
######### Done investigating Ash

