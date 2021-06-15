## Evaluating the distribution of data
library(here)
here()
#here::i_am(here("TraitAnalyses201003/Code_10032020","Histogram_on_GEBVs_and_Ones_used_in_MiniCross.R"))

# Read in the BLUPs file
# Read in the link file

library(stringr)
find_par<-function(GPstring){
  SPs<-strsplit(as.character(GPstring), split="-", fixed=T)
  Par<-sapply(SPs, function(vec) paste(vec[1:3], collapse="-"))
  return(Par)
}

MiniCrossList<-read.csv(here("TraitAnalyses201003","MiniCross_status_02052021.csv"))
  str(MiniCrossList)
  
MiniCrossList$FemalePar<-find_par(MiniCrossList$Female.GP.ID)
MiniCrossList$MalePar<-find_par(MiniCrossList$Cross.To.Male_GP.ID)
  head(MiniCrossList)


MiniCross2019SPs<-unique(c(MiniCrossList$FemalePar,MiniCrossList$MalePar))
Rename<-as.data.frame(str_split_fixed(MiniCross2019SPs,"-",3))
Rename$V4<-paste0("S",Rename[,3])

MiniCross2019SPs<-paste(Rename[,1],Rename[,2],Rename[,4],sep ="-")

#### Use the updated file on the actual GOs being used for GOM_2020 Crosses.
SPCross2020<-read.csv(here("TraitAnalyses201003","2020_Actual_GOM_Crosses.csv"),sep=",",header=T)
GPs_Used<-unique(c(as.character(SPCross2020$FG),as.character(SPCross2020$MG)))
SP2019_GPs<-GPs_Used[stringr::str_detect(GPs_Used,"SL18-UCONN-",negate=FALSE)]
SP2019_plots<-unique(find_par(SP2019_GPs))  #11


#paste0(substring(MiniCross2019SPs, first=1, last=11), "S", substring(MiniCross2019SPs, first=12))

Plot_Cross_Link<- read.csv(here("TraitAnalyses201003","Plot_num_Cross_Link_02052021.csv"),sep=",",header=T) 
  head(Plot_Cross_Link)   
Plot_Cross_Link$Used_in_2020SPCross<-ifelse(Plot_Cross_Link$crossID%in%SP2019_plots,"y","no")  # An update
  
Plot_Cross_Link$Used_in_2021miniCross<-ifelse(Plot_Cross_Link$crossID%in%MiniCross2019SPs,"y","no")

##### Find the BLUPs file,

datafdr<-here(here(),"TraitAnalyses201003/data/")
#plotBLUPs<-read.csv(here("allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_0116_2021_hap.csv"),sep=",",header=T,row.names=1)
plotBLUPs<-read.csv(paste0(datafdr,"allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_0315_2021_hap_NonUpateBladeDensity.csv"),sep=",",header=T,row.names=1)

  head(plotBLUPs)
##### Then plot out
plotBLUPs$Isolate_Sorus<-expss::vlookup(rownames(plotBLUPs),dict=Plot_Cross_Link,result_column = "Isolate_Sorus_2020_Need_Update",lookup_column = "Crosses")
  head(plotBLUPs)

  unique(plotBLUPs$Crossed2020_Update)
plotBLUPs$Crossed2020<-expss::vlookup(rownames(plotBLUPs),dict=Plot_Cross_Link,result_column = "Used_in2020Cross",lookup_column ="Crosses" )  
plotBLUPs$Crossed2021<-expss::vlookup(rownames(plotBLUPs),dict=Plot_Cross_Link,result_column = "Used_in_2021miniCross",lookup_column ="Crosses" )
plotBLUPs$Crossed2020_Update<-expss::vlookup(rownames(plotBLUPs),dict=Plot_Cross_Link,result_column = "Used_in_2020SPCross",lookup_column ="Crosses" )

## Change trait here !!!!!!!
forPlot<-data.frame(Indiv=rownames(plotBLUPs),Trait=plotBLUPs$DWpM,Isolate_Sorus=plotBLUPs$Isolate_Sorus,Crossed2020=plotBLUPs$Crossed2020,Crossed2021=plotBLUPs$Crossed2021,Crossed2020_Update=plotBLUPs$Crossed2020_Update)

## Detect if it has "x"= SP plots
forPlot<-forPlot[stringr::str_detect(forPlot$Indiv,"x",negate=FALSE),]

library(ggplot2)
histo<-ggplot(forPlot,aes(x=Trait))+
  geom_histogram(binwidth=0.01)+
  labs(x="Trait GEBV",y="Plot Count")+
  theme_bw()+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"))

  

# segment_data = data.frame(
#   x = c(forPlot$Trait[forPlot$Crossed2020=="y"]),
#   xend = c(forPlot$Trait[forPlot$Crossed2020=="y"]), 
#   y = rep(0,length(forPlot$Trait[forPlot$Crossed2020=="y"])),
#   yend = rep(4, length(forPlot$Trait[forPlot$Crossed2020=="y"]))
# )
# its 2019_SP plot level GPs used for Plot Crossing in 2020 
segment_data = data.frame(
  x = c(forPlot$Trait[forPlot$Crossed2020_Update=="y"]),
  xend = c(forPlot$Trait[forPlot$Crossed2020_Update=="y"]), 
  y = rep(0,length(forPlot$Trait[forPlot$Crossed2020_Update=="y"])),
  yend = rep(4, length(forPlot$Trait[forPlot$Crossed2020_Update=="y"]))
)


# 2020 harvesting had sorus tissue
segment_data2 = data.frame(
  x = c(forPlot$Trait[forPlot$Isolate_Sorus=="Yes"]),
  xend = c(forPlot$Trait[forPlot$Isolate_Sorus=="Yes"]),
  y = rep(0,length(forPlot$Trait[forPlot$Isolate_Sorus=="Yes"])),
  yend = rep(4.5, length(forPlot$Trait[forPlot$Isolate_Sorus=="Yes"]))
)


# its GPs used for Mini Cross in 2021
segment_data3 = data.frame(
  x = c(forPlot$Trait[forPlot$Crossed2021=="y"]),
  xend = c(forPlot$Trait[forPlot$Crossed2021=="y"]), 
  y = rep(0,length(forPlot$Trait[forPlot$Crossed2021=="y"])),
  yend = rep(4, length(forPlot$Trait[forPlot$Crossed2021=="y"]))
)

sum(forPlot$Crossed2020=="y")  # data  13
sum(forPlot$Isolate_Sorus=="Yes") # data2  35
sum(forPlot$Crossed2021=="y")  # data 3 25
sum(forPlot$Crossed2020_Update=="y") # data  11

colors <- c("Produced sorus tissue in Yr2019" = "green", "Produced sorus tissue in Yr2020" = "red")

histo +
geom_segment(data=segment_data,aes(x = x, y = y, xend = xend, yend = yend,color="green"),linetype="dashed",show.legend=TRUE)+
geom_segment(data=segment_data3,aes(x = x, y = y, xend = xend, yend = yend,color="red"),linetype="dashed",show.legend=TRUE)+
   scale_color_manual(name = "",
                        values = c( "green", "red"),
                        labels = c("Selected in Yr2019", "Selected in Yr2020"))
  
#!!! The color=() must be inside the aes() statement, otherwise, legend cannot be viewed

#geom_segment(data = segment_data2, aes(x = x, y = y, xend = xend, yend = yend),linetype="dashed",color="steelblue")+
