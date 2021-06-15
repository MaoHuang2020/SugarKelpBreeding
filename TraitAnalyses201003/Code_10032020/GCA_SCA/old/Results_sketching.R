#Results comparison
#5 Reps, cor(BLUEs and GCA1+GCA2+SCA)

#1. All data predict itself
0.936

#5REPs varE   varB_ETA3    varGCA1     varGCA2     varSCA 
# 0.052979270 0.011043716 0.002925032 0.001710788 0.003170026

### (With and without femaLoc and maleLoc)


# 2. All data CV_10 fold


#3. BothGPs geno predict others
-0.070

#5REPs varE   varB_ETA3     varGCA1     varGCA2      varSCA 
# 0.070814898 0.010794574 0.001855974 0.002027203 0.001677588
 
#4. FG geno predict others (TP: 49!) 
#---- Q: can we only gentoype FGs? construct the A matrix/GCA1, differently, removing those Males genotyped "pretend no info"
print(mean(r))
0.076

#5REPs varE   varB_ETA3     varGCA1     varGCA2      varSCA 
#0.009133676 0.007052382 0.002248628 0.003859018 0.003592970

#5. yr19_20  (RM 5 common plots)
print(mean(r))
-0.116

#5REPs varE    varB_ETA3      varGCA1      varGCA2       varSCA 
#0.0145761019 0.0014965818 0.0002895085 0.0002806575 0.0001338398 

#6. yr20_19  (Include 5 common plots)
print(mean(r))
-0.028

#5REPs varE   varB_ETA3     varGCA1     varGCA2      varSCA 
#0.065364183 0.025325891 0.004340126 0.003323680 0.007853688 

## ASreml between Year genetic correlation, sqrt(10)*phenotypes of Year1 
0.491
#                                             component  std.error  z.ratio
# Year:vm(Crosses, Trait_grm3)!Year_2019:2019 0.04636942 0.02158158 2.148565
# Year:vm(Crosses, Trait_grm3)!Year_2020:2019 0.03428150 0.03408872 1.005655
# Year:vm(Crosses, Trait_grm3)!Year_2020:2020 0.10528296 0.03528587 2.983714
# units!R                                     0.10558092 0.01407421 7.501729


# 2. Between Loc prediction
##stdE
# CB         CC         JS         LD         NC         NL         
# 0.12432809 0.07807229 0.03218488 0.02579745 0.05501534 0.02505755  
# OI          SF 
# 0.07120138  0.10800176 
##Mean on BLUE
# CB          CC          JS          LD          NC          NL 
# 0.33343530  0.34441539 -0.08539057 -0.19884729 -0.08185030 -0.11030479 
# OI          SF 
# 0.03992838  0.09420910

#2. Predict that Loc
##Cor using GCA+SCA, 5 Reps
CB          CC          JS          LD          NC          NL
0.016     -0.088       0.266        0.233     0.093       0.214
OI          SF
0.067    -0.092

##Cor using yHat, 5 Reps
# CB          CC          JS          LD          NC          NL 
# 0.54001055 -0.09604023 -0.27680635  0.59474666  0.12345266  0.03154891 
# OI          SF 
# 0.37882775  0.51000054  




# All data itself cor: 0.962
#Yr19to20  cor:-0.1417018 / -0.0768
#Yr20to19  cor:0.1147057 /-0.1160252/ -0.106/ 0.05

#Loc to others cor: -0.08165424/ 0.082
#BothGP genotyped cor: 0.11

# Change G1, so the GCA1 will be different?


## The averages of 10K samples for Yr20 to pedict Yr19 
#  varE        varB       varU3       varU4       varU5 
#  0.0609   0.01738     0.00896     0.00376     0.00869

## The fm output  
# [1] 0.0582
# [1] 0.0173
# [1] 0.0089
# [1] 0.0038
# [1] 0.0090

## The averages of 10K samples for Yr19 to predict Yr20
# varE         varB        varU3        varU4        varU5 
# 0.00651      0.00115     0.00044    0.00044       0.00016

# The fm output
# varE
#  0.00634
# ETA[[2]]$varB
#  0.00116
# ETA[[3]]$varU
#  0.00042
# ETA[[4]]$varU
#  0.00045
# ETA[[5]]$varU
#  0.00018



### Alldata_output from 10K samples
# > mean(varE)
# [1] 0.04635725
# > mean(varB)
# [1] 0.01075573
# > mean(varU3)
# [1] 0.005894885
# > mean(varU4)
# [1] 0.001956617
# > mean(varU5)
# [1] 0.003598225

