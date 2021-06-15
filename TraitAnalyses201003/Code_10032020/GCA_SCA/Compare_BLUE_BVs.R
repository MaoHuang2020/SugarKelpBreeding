# Compare BLUEs vs BLUPs
# Compared files were saved in the "/Users/maohuang/Desktop/Kelp/GCA_SCA/OneTime1920"
### Y????
####
WD<-"/Users/maohuang/Desktop/Kelp/GCA_SCA/" # local
load(paste0(WD,"OneTime1920/data/","dataNHpi_withChk_3_sets_PhotoScore23.rdata")) 
Y1<-dataNHpiBoth_C
colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM","AshFreedryWgtPerM")
Y2<-droplevels(Y1[,colKeep])
head(Y2)
Y<-Y2 # Both Years


load(file=paste0(WD,"OneTime1920/data/","BLUE_DwPM_2vs1Year_Update03082021.rdata"))
ls()
  head(CrossBLUE)

wd<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/"
allBVs<-read.csv(paste0(wd,"allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_0116_2021_dip.csv"),sep=",",header=T)

  head(allBVs)
  dim(CrossBLUE)
  dim(allBVs)

rownames(CrossBLUE)<-CrossBLUE$CrossName  

CrossBLUE$BV_DWpM<-expss::vlookup(CrossBLUE$CrossName,dict=allBVs,lookup_column = "X",result_column = "DWpM")
write.csv(CrossBLUE,"Crosses_BLUE_and_BV_both_Years_Data.csv")

head(CrossBLUE)
cor(CrossBLUE$BV_DWpM,CrossBLUE$BLUE_DwPM_2Yrs,use="complete")  # 0.67

Y$BLUE_Trait_WithinYr<-ifelse(Y$Year==2019,
                     expss::vlookup(Y$Crosses,dict=CrossBLUE,result_column = "BLUE_DwPM_2019",lookup_column = "CrossName"),
                     expss::vlookup(Y$Crosses,dict=CrossBLUE,result_column = "BLUE_DwPM_2020",lookup_column = "CrossName"))
Y$BLUE_Trait_BothYr<-expss::vlookup(Y$Crosses,dict=CrossBLUE,result_column = "BLUE_DwPM_2Yrs",lookup_column = "CrossName")
### This is the DwPM
  head(Y)

library(ggplot2)
plot<-ggplot(data=Y,aes(BLUE_Trait_WithinYr,Year))+
  geom_point(aes(color=as.factor(Year)))+ 
  geom_line(aes(group=as.factor(Crosses)))
print(plot)
