### To JL on haploid and diploid data files

library(here)
load(here("TraitAnalyses201003/data","outCovComb_dip_0116_2021.Rdata"))
### The outCovComb4_dipOrder is the diploid relationship matrix

hMat_dip<-outCovComb4_dipOrder

load(here("TraitAnalyses201003/data","outCovComb4_hap_Conden_0116_2021.Rdata")) 
#"hMat_hap_0116_2021.Rdata" saved the hMat_hap=outCovComb4_Hapconden
hMat-hap<-outCovComb4_Hapconden   

rownames(hMat_hap)<-rownames(hMat_dip)  # order of hMat_hap = biphasichapccMat=biphasicPedNH=hMat_dip
colnames(hMat_hap)<-colnames(hMat_dip)

cultevo::mantel.test(dist(as.matrix(hMat_dip)),dist(as.matrix(hMat_hap)),trials=99) # r =0.959  

# row1<-1:56   #fndr+its GP both genotyped
# row2<-57:59  # fndr genotyped, its GP Not
# row3<-60:93  # fndr Not genotyped, its GP yes
# row4<-94:104 # fndr and GP Both Not genotyped
# row5<-c(105:263,332:442)  # GPs genotyped
# row6<-c(264:331,443:543,789:866)  # GPs Not genotyped
# row7<-c(544:788)  # SP rows
# row8<-c(789:866)  # GPs_from_Farm_SP rows



###### Below is how I calculated cor for each catogery
###### Some categories (GP1-5, SP1-2) were old definitions I have in excel for explaination 
corM<-function(n1,n2,matrx1,matrx2){
  if(n1==n2){
    cor(c(matrx1[n1,]),c(matrx2[n1,]),use="complete")
  }else{
    cor(c(matrx1[n1:n2,]),c(matrx2[n1:n2,]),use="complete")
  }
}

# with diagonal  !!!
low_hapM<-hMat_hap;low_hapM[upper.tri(low_hapM)]<-NA
low_dipM<-hMat_dip;low_dipM[upper.tri(low_dipM)]<-NA

# without diagonal  !!!
low_hapM2<-low_hapM;diag(low_hapM2)<-NA
low_dipM2<-low_dipM;diag(low_dipM2)<-NA

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

colnames(cormatrix)<-c("samplesize","WithDiagonal","noDiagonal")
rownames(cormatrix)<-c("fndrGP_geno","fndrGeno_GPNo","fndrNo_GPGeno","fndrGP_No","GP1","GP2","GP3","GP5","SP1","SP2","SP_GP","All")
  cormatrix


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

