# Make hist on the plots
# Record the ones has the sorus tissue

# Find the ones that are the same cross between both years based on the Crosses INFO?????

#allBLUPs<-read.csv("allBLUPsDF_Ordred_Both_Plots.csv",sep=",",header=T,row.names=1)

progenySP <- which(nchar(rownames(allBLUPs)) > 15)
SP_allBLUPs<-as.data.frame(allBLUPs[progenySP,])
SP_allBLUPs<-SP_allBLUPs[order(-SP_allBLUPs$AshFDwPM),]
  dim(SP_allBLUPs)

plotBLUPs<-SP_allBLUPs
  head(plotBLUPs)
library(ggplot2)  

### Had GPs biomass, being used in making Crosses this year
Plot_Cross_Link<- read.csv("Plot_num_Cross_Link.csv",sep=",",header=T) 
  head(Plot_Cross_Link)    
  
library(expss)  
plotBLUPs$Isolate_Sorus<-vlookup(rownames(plotBLUPs),dict=Plot_Cross_Link,result_column = "Isolate_Sorus_2020_Need_Update",lookup_column = "Crosses")
  head(plotBLUPs)
  
plotBLUPs$Crossed2020<-vlookup(rownames(plotBLUPs),dict=Plot_Cross_Link,result_column = "Used_in2020Cross",lookup_column ="Crosses" )  

histo<-ggplot(plotBLUPs,aes(x=AshFDwPM))+
  geom_histogram(binwidth=0.01)

segment_data = data.frame(
  x = c(plotBLUPs$AshFDwPM[plotBLUPs$Isolate_Sorus=="Yes"]),
  xend = c(plotBLUPs$AshFDwPM[plotBLUPs$Isolate_Sorus=="Yes"]), 
  y = rep(0,length(plotBLUPs$AshFDwPM[plotBLUPs$Isolate_Sorus=="Yes"])),
  yend = rep(5, length(plotBLUPs$AshFDwPM[plotBLUPs$Isolate_Sorus=="Yes"]))
)

segment_data2 = data.frame(
  x = c(plotBLUPs$AshFDwPM[plotBLUPs$Crossed2020=="y"]),
  xend = c(plotBLUPs$AshFDwPM[plotBLUPs$Crossed2020=="y"]), 
  y = rep(0,length(plotBLUPs$AshFDwPM[plotBLUPs$Crossed2020=="y"])),
  yend = rep(4, length(plotBLUPs$AshFDwPM[plotBLUPs$Crossed2020=="y"]))
)

 histo +
  geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),linetype="dashed",color="green")+
   geom_segment(data=segment_data2,aes(x = x, y = y, xend = xend, yend = yend),linetype="dashed",color="red")+
   labs(x="AshFDwPM",y="Plot Count")
  
 
 
 
#### PDW  
 histo<-ggplot(plotBLUPs,aes(x=PDW))+
   geom_histogram(binwidth=0.01)
 
 segment_data = data.frame(
   x = c(plotBLUPs$PDW[plotBLUPs$Isolate_Sorus=="Yes"]),
   xend = c(plotBLUPs$PDW[plotBLUPs$Isolate_Sorus=="Yes"]), 
   y = rep(0,length(plotBLUPs$PDW[plotBLUPs$Isolate_Sorus=="Yes"])),
   yend = rep(10, length(plotBLUPs$PDW[plotBLUPs$Isolate_Sorus=="Yes"]))
 )
 
 histo +
   geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),linetype="dashed",color="red")+
   labs(x="PDW",y="Plot Count")

 #### DWpM  
 histo<-ggplot(plotBLUPs,aes(x=DWpM))+
   geom_histogram(binwidth=0.01)
 
 segment_data = data.frame(
   x = c(plotBLUPs$DWpM[plotBLUPs$Isolate_Sorus=="Yes"]),
   xend = c(plotBLUPs$DWpM[plotBLUPs$Isolate_Sorus=="Yes"]), 
   y = rep(0,length(plotBLUPs$DWpM[plotBLUPs$Isolate_Sorus=="Yes"])),
   yend = rep(5, length(plotBLUPs$DWpM[plotBLUPs$Isolate_Sorus=="Yes"]))
 )
 
 histo +
   geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend),linetype="dashed",color="red")+
   labs(x="DWpM",y="Plot Count")
 
 
# geom_vline(xintercept=c(plotBLUPs$index[plotBLUPs$Isolate_Sorus=="Yes"]),linetype="dashed",color="red")+
# geom_abline(xintercept=plotBLUPs$index[plotBLUPs$Isolate_Sorus=="Yes"],col="red")
#hist(plotBLUPs$index)
#abline(v=c(plotBLUPs$index[plotBLUPs$Isolate_Sorus=="Yes"]),col="red")

 ###h2PlIn<-read.csv("Heritabilities_Plot_Indi.csv",sep=",",header=T)
 
 #### Manually Input all the heritabilities
 # h2_19<-read.csv("heritability_all_2019.csv",sep=",",header=T)
 # h2_20<-read.csv("heritability_all_2020.csv",sep=",",header=T)
 # h2_Both<-read.csv("heritability_all_Both.csv",sep=",",header=T)
 # 
 
h2<-read.csv("Heritabilities_2019_2020_Both_with_0116_2021_outCovComb.csv",sep=",",header=T)

 h2_Both$h_2019<-vlookup(h2_Both$X,dict=h2_19,result_column = "x",lookup_column = "X")
 h2_Both$h_2020<-vlookup(h2_Both$X,dict=h2_20,result_column = "x",lookup_column = "X")
  head(h2_Both)
  head(h2_20)
  head(h2_19)
rownames(h2_Both)<-h2_Both$X
h2_Both<-h2_Both[,-1]

h2PlIn<-as.data.frame(t(h2_Both))
h2PlIn$Data<-c("Both","2019","2020")
  head(h2PlIn)
  str(h2PlIn) 
  
write.csv(h2PlIn,"heritabilities_h2PlIn.csv")  
###### Made the h2PlIn the format it needs
  
 #wide to long format
 library(tidyr)
# Gather different traits,what is their value called, from col1:colN 
h2PI<-gather(h2PlIn,Trait,Herit,AshFDwPM:stipeDiameter,factor_key = T)
 head(h2PI)

levels(h2PI$Trait) <-c("AshFDwPM","Ash","WWP","DWpM","pDW","BD","BL","BmWid","BTh","SL","SDia")
  head(h2PI)

ggplot(h2PI,aes(fill=Data,y=Herit,x=Trait))+
      geom_bar(position=position_dodge(),stat="identity")
      

# 
h2PL<-read.csv("Heritabilities_Plot_Scenarios.csv",sep=",",header=T)
  head(h2PL) 

h2PL.long<-gather(h2PL,Trait,Herit,densityBlades:percDryWgt,factor_key=T)   
  head(h2PL.long)

##compared Amat vs Hmat (both using Model2, where densityBlades are conditioned)
h2AH<-h2PL.long[h2PL.long$Model==2,]
h2AH<-h2AH[!is.na(h2AH$Herit),]
h2AH

ggplot(h2AH,aes(fill=Matrix,y=Herit,x=Trait))+
  geom_bar(position=position_dodge(),stat="identity")+
    facet_grid(rows=vars(Data))
## compared model2 (condition on blade density) and model3 (condition on development stage)

h2M23<-h2PL.long[h2PL.long$Scenario=="c" |h2PL.long$Scenario=="d",]
h2M23<-h2M23[!h2M23$Trait=="densityBlades",]
  h2M23
h2M23$Scenario<-factor(h2M23$Scenario)
  str(h2M23)
levels(h2M23$Scenario)<-c("BladeDensity covariate","Development covariate")

ggplot(h2M23,aes(fill=Scenario,y=Herit,x=Trait))+
  geom_bar(position=position_dodge(),stat="identity")
