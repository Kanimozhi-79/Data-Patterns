########################
##LOADING LIBRARY FILES
#######################

library(infotheo)
library(discretization)
library(arules)
library(microbenchmark)
library(caret)
library(e1071)
library(mlbench)
library(kernlab)
library(xlsx)
library(arulesCBA)
library(Matrix)
library(purrr)
library(dplyr)
library(qpcR) 

######################################
#     Reading Input data Partition
######################################
###       1.      "F:/Ph.D Work/R/Coding/GDS4275-PituitaryTumor/GDS4275_Pituitary Tumor_1_8310genes.csv"
###       2.      "F:/Ph.D Work/R/Coding/GDS4275-PituitaryTumor/GDS4275_Pituitary Tumor_2_8472genes.csv"
###       3.      "F:/Ph.D Work/R/Coding/GDS4275-PituitaryTumor/GDS4275_Pituitary Tumor_3_8497genes.csv"
###       4.      "F:/Ph.D Work/R/Coding/GDS4275-PituitaryTumor/GDS4275_Pituitary Tumor_4_8796genes.csv"
###       5.      "F:/Ph.D Work/R/Coding/GDS4275-PituitaryTumor/GDS4275_Pituitary Tumor_5_4854genes.csv"
###                Genes---38932
###               Samples---23
###               Tumor --  1:14
###               Normal--  15:23
start_time <- Sys.time()
a<-read.csv("F:/Ph.D Work/R/Coding/GDS4275-PituitaryTumor/GDS4275_Pituitary Tumor_1_8310genes.csv",header=TRUE)
b<-read.csv("F:/Ph.D Work/R/Coding/GDS4275-PituitaryTumor/GDS4275_Pituitary Tumor_2_8472genes.csv",header=TRUE)
c<-read.csv("F:/Ph.D Work/R/Coding/GDS4275-PituitaryTumor/GDS4275_Pituitary Tumor_3_8497genes.csv",header=TRUE)
d<-read.csv("F:/Ph.D Work/R/Coding/GDS4275-PituitaryTumor/GDS4275_Pituitary Tumor_4_8796genes.csv",header=TRUE)
e<-read.csv("F:/Ph.D Work/R/Coding/GDS4275-PituitaryTumor/GDS4275_Pituitary Tumor_5_4854genes.csv",header=TRUE)


inputdata1=data.frame(matrix(NA,nrow=23,ncol=0))
inputdata1<-cbind(inputdata1,a)
inputdata1<-cbind(inputdata1,b)
inputdata1<-cbind(inputdata1,c)
inputdata1<-cbind(inputdata1,d)
inputdata1<-cbind(inputdata1,e)

ncol(inputdata1)
ncol(a)
ncol(b)
ncol(c)
ncol(d)
ncol(e)

inputdata1=round(inputdata1,digits=1)
inputdata1<-as.data.frame(sapply(inputdata1,as.numeric))

rows=nrow(inputdata1)#Number of Rows
cols=ncol(inputdata1)#Number of Columns

#omitting Label Column ie., only features columns

##################################################
#    # Creating Column vector to store Gene Names
##################################################

Col_Name_EO<-c()

Col_Name_CE<-c()
Col_Name_CE_NT<-c() #Col_Name_PO_NT
Col_Name_CE_TN<-c()

Col_Name_TBE<-c()
Col_Name_TBE_NT<-c() #Col_Name_PO_NT
Col_Name_TBE_TN<-c()

Col_Name_BBE<-c()
Col_Name_BBE_NT<-c() #Col_Name_PO_NT
Col_Name_BBE_TN<-c()

New_DF_EO=data.frame(matrix(NA,nrow=23,ncol=0))

New_DF_CE=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_CE_NT=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_CE_TN=data.frame(matrix(NA,nrow=23,ncol=0))

New_DF_TBE=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_TBE_NT=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_TBE_TN=data.frame(matrix(NA,nrow=23,ncol=0))

New_DF_BBE=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_BBE_NT=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_BBE_TN=data.frame(matrix(NA,nrow=23,ncol=0))



#Partial Overlapping
#coMBINED NT AND TN PO
Col_Name_PO<-c()
#Col_Name_PO
#SEPARATE NT AND TN PO
Col_Name_PO_NT<-c()
#Col_Name_PO_NT
Col_Name_PO_TN<-c()
#Col_Name_PO_TN

#Top Bottom Boundary Overlap
#coMBINED NT AND TN TBBO
Col_Name_TBBO<-c()
#Col_Name_TBBO
#SEPARATE NT AND TN PO
Col_Name_TBBO_NT<-c()
#Col_Name_TBBO_NT
Col_Name_TBBO_TN<-c()
#Col_Name_TBBO_TN

#Complete Separation
#coMBINED NT AND TN CS
Col_Name_CS<-c()
Col_Name_CS
#SEPARATE NT AND TN PO
Col_Name_CS_NT<-c()
#Col_Name_CS_NT
Col_Name_CS_TN<-c()
#Col_Name_CS_TN
##################################################
#New_DF_EO=data.frame(matrix(NA,nrow=64,ncol=0))
# Creating new DataFrame to append partially overlapping classes and separation classes
#change nrow when Dataset changed
New_DF_PO=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_Names=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_PO_TN=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_PO_NT=data.frame(matrix(NA,nrow=23,ncol=0))

# Creating new DataFrame to append cOMPLETELY separation classes
#change nrow when Dataset changed
New_DF_TBBO=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_TBBO_Names=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_TBBO_TN=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_TBBO_NT=data.frame(matrix(NA,nrow=23,ncol=0))

# Creating new DataFrame to append cOMPLETELY separation classes
#change nrow when Dataset changed
New_DF_CS=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_CS_Names=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_CS_NT=data.frame(matrix(NA,nrow=23,ncol=0))
New_DF_CS_TN=data.frame(matrix(NA,nrow=23,ncol=0))

#Vectors for storing N_min, T_min, N_max & T_max
#I.EXACT Overlap NT
N_min_EO<-c()
N_max_EO<-c()
T_min_EO<-c()
T_max_EO<-c()
#II.Complete Embedding_NT
N_NT_min_CE<-c()
N_NT_max_CE<-c()
T_NT_min_CE<-c()
T_NT_max_CE<-c()

#III.Complete Embedding_TN
N_TN_min_CE<-c()
N_TN_max_CE<-c()
T_TN_min_CE<-c()
T_TN_max_CE<-c()

N_min_CE<-c()
N_max_CE<-c()
T_min_CE<-c()
T_max_CE<-c()

#IV.LMO-UME_NT
N_NT_min_TBE<-c()
N_NT_max_TBE<-c()
T_NT_min_TBE<-c()
T_NT_max_TBE<-c()

#V. LMO-UME_TN
N_TN_min_TBE<-c()
N_TN_max_TBE<-c()
T_TN_min_TBE<-c()
T_TN_max_TBE<-c()

N_min_TBE<-c()
N_max_TBE<-c()
T_min_TBE<-c()
T_max_TBE<-c()

#VI.UMO-LME_NT
N_NT_min_BBE<-c()
N_NT_max_BBE<-c()
T_NT_min_BBE<-c()
T_NT_max_BBE<-c()

#VII. UMO-LME_TN
N_TN_min_BBE<-c()
N_TN_max_BBE<-c()
T_TN_min_BBE<-c()
T_TN_max_BBE<-c()

N_min_BBE<-c()
N_max_BBE<-c()
T_min_BBE<-c()
T_max_BBE<-c()

#VIII.Partial Overlap NT
N_NT_min_PO<-c()
N_NT_max_PO<-c()
T_NT_min_PO<-c()
T_NT_max_PO<-c()

#IX. Partial Overlap TN
N_TN_min_PO<-c()
N_TN_max_PO<-c()
T_TN_min_PO<-c()
T_TN_max_PO<-c()

N_min_PO<-c()
N_max_PO<-c()
T_min_PO<-c()
T_max_PO<-c()

#X. TBBO NT
N_NT_min_TBBO<-c()
N_NT_max_TBBO<-c()
T_NT_min_TBBO<-c()
T_NT_max_TBBO<-c()

#XI.TBBO NT
N_TN_min_TBBO<-c()
N_TN_max_TBBO<-c()
T_TN_min_TBBO<-c()
T_TN_max_TBBO<-c()

N_min_TBBO<-c()
N_max_TBBO<-c()
T_min_TBBO<-c()
T_max_TBBO<-c()

#XII.CS NT
N_NT_min_CS<-c()
N_NT_max_CS<-c()
T_NT_min_CS<-c()
T_NT_max_CS<-c()

#XIII.CS TN
N_TN_min_CS<-c()
N_TN_max_CS<-c()
T_TN_min_CS<-c()
T_TN_max_CS<-c()

N_min_CS<-c()
N_max_CS<-c()
T_min_CS<-c()
T_max_CS<-c()

#Storing Column names 
Gene_ID<-colnames(inputdata1)

i=j=k=l=m=n=o=p=q=r=s=1
partial1=partial2=complete1=complete2=1
inr=1

EO_inr=1
CE_inr=1
TBE_inr=1
BBE_inr=1
PO_inr=1
CS_inr=1
TBBO_inr=1

CE_NT_inr=1
CE_TN_inr=1
TBE_NT_inr=1
TBE_TN_inr=1
BBE_NT_inr=1
BBE_TN_inr=1

PO_NT_inr=1
PO_TN_inr=1
TBBO_NT_inr=1
TBBO_TN_inr=1
CS_NT_inr=1
CS_TN_inr=1


start_time1 <- Sys.time()
for(i in 1:cols){ #n Number of features
  Tumor_min<-round(min(inputdata1[1:14,i]),digits=1) #No of Tumor Samples
  Tumor_max<-round(max(inputdata1[1:14,i]),digits=1)
  Tumor_min=as.numeric(Tumor_min)
  Tumor_max=as.numeric(Tumor_max)
  
  Normal_min<-round(min(inputdata1[15:23,i]),digits=1) #No of Normal Samples
  Normal_max<-round(max(inputdata1[15:23,i]),digits=1)
  Normal_min=as.numeric(Normal_min)
  Normal_max=as.numeric(Normal_max)
  # 
  # Tmin1<-Tumor_min
  # Tmax1<-Tumor_max
  # Nmin1<-Normal_min
  # Nmax1<-Normal_max
  # 
  # 
  #   
  
  # Case 1)Exact Overlap
  if(Tumor_min==Normal_min && Tumor_max==Normal_max){ 
    N_min_EO[EO_inr]=Normal_min
    N_max_EO[EO_inr]=Normal_max
    T_min_EO[EO_inr]=Tumor_min
    T_max_EO[EO_inr]=Tumor_max
    New_DF_EO=cbind(New_DF_EO,inputdata1[,i])
    Col_Name_EO[EO_inr]<-names(inputdata1[i])
    
    j=j+1
    inr=inr+1
    EO_inr=EO_inr+1
  }   
  #Case 2a)Complete embedding Tumor inside Normal
  else if(Normal_min<Tumor_min && Tumor_max>Normal_min && Normal_max>Tumor_min && Normal_max>Tumor_max){
    N_NT_min_CE[CE_NT_inr]=Normal_min
    N_NT_max_CE[CE_NT_inr]=Normal_max
    T_NT_min_CE[CE_NT_inr]=Tumor_min
    T_NT_max_CE[CE_NT_inr]=Tumor_max
    New_DF_CE_NT=cbind(New_DF_CE_NT,inputdata1[,i])
    Col_Name_CE_NT[CE_NT_inr]<-names(inputdata1[i])
   
    N_min_CE[CE_inr]=Normal_min
    N_max_CE[CE_inr]=Normal_max
    T_min_CE[CE_inr]=Tumor_min
    T_max_CE[CE_inr]=Tumor_max
    New_DF_CE=cbind(New_DF_CE,inputdata1[,i])
    Col_Name_CE[CE_inr]<-names(inputdata1[i])
   
     k=k+1
    inr=inr+1
    CE_inr=CE_inr+1
    CE_NT_inr=CE_NT_inr+1
  }
  #Case 2b)Complete embedding  Normal inside Tumor
  else if(Normal_min>Tumor_min && Normal_min<Tumor_max && Normal_max>Tumor_min && Normal_max<Tumor_max){
    N_TN_min_CE[CE_TN_inr]=Normal_min
    N_TN_max_CE[CE_TN_inr]=Normal_max
    T_TN_min_CE[CE_TN_inr]=Tumor_min
    T_TN_max_CE[CE_TN_inr]=Tumor_max
    New_DF_CE_TN=cbind(New_DF_CE_TN,inputdata1[,i])
    Col_Name_CE_TN[CE_TN_inr]<-names(inputdata1[i])
    
    N_min_CE[CE_inr]=Normal_min
    N_max_CE[CE_inr]=Normal_max
    T_min_CE[CE_inr]=Tumor_min
    T_max_CE[CE_inr]=Tumor_max
    New_DF_CE=cbind(New_DF_CE,inputdata1[,i])
    Col_Name_CE[CE_inr]<-names(inputdata1[i])
    
     l=l+1
     inr=inr+1
     CE_inr=CE_inr+1
     CE_TN_inr=CE_TN_inr+1
    }      
  #Case 3a)Exact Top Boundary and Complete Embedding Tumor inside Normal
  else if(Tumor_min==Normal_min && Tumor_max<Normal_max){
    N_NT_min_TBE[TBE_NT_inr]=Normal_min
    N_NT_max_TBE[TBE_NT_inr]=Normal_max
    T_NT_min_TBE[TBE_NT_inr]=Tumor_min
    T_NT_max_TBE[TBE_NT_inr]=Tumor_max
    New_DF_TBE_NT=cbind(New_DF_TBE_NT,inputdata1[,i])
    Col_Name_TBE_NT[TBE_NT_inr]<-names(inputdata1[i])
    
    N_min_TBE[TBE_inr]=Normal_min
    N_max_TBE[TBE_inr]=Normal_max
    T_min_TBE[TBE_inr]=Tumor_min
    T_max_TBE[TBE_inr]=Tumor_max
    New_DF_TBE=cbind(New_DF_TBE,inputdata1[,i])
    Col_Name_TBE[TBE_inr]<-names(inputdata1[i])
    
    m=m+1
    inr=inr+1
    TBE_inr=TBE_inr+1
    TBE_NT_inr=TBE_NT_inr+1
    
  }       
  #Case 3b)Exact Top Boundary and Complete Embedding Normal inside Tumor
  else if(Tumor_min==Normal_min &&  Normal_max<Tumor_max){ 
    N_TN_min_TBE[TBE_TN_inr]=Normal_min
    N_TN_max_TBE[TBE_TN_inr]=Normal_max
    T_TN_min_TBE[TBE_TN_inr]=Tumor_min
    T_TN_max_TBE[TBE_TN_inr]=Tumor_max
    New_DF_TBE_TN=cbind(New_DF_TBE_TN,inputdata1[,i])
    Col_Name_TBE_TN[TBE_TN_inr]<-names(inputdata1[i])
    
    N_min_TBE[TBE_inr]=Normal_min
    N_max_TBE[TBE_inr]=Normal_max
    T_min_TBE[TBE_inr]=Tumor_min
    T_max_TBE[TBE_inr]=Tumor_max
    New_DF_TBE=cbind(New_DF_TBE,inputdata1[,i])
    Col_Name_TBE[TBE_inr]<-names(inputdata1[i])
    
    n=n+1
    inr=inr+1
    TBE_inr=TBE_inr+1
    TBE_TN_inr=TBE_TN_inr+1
   }
  #Case 4a)Exact Bottom Boundary & Complete Embedding Tumor inside Normal
  else if(Tumor_max==Normal_max && Tumor_min<Normal_min){
    N_NT_min_BBE[BBE_NT_inr]=Normal_min
    N_NT_max_BBE[BBE_NT_inr]=Normal_max
    T_NT_min_BBE[BBE_NT_inr]=Tumor_min
    T_NT_max_BBE[BBE_NT_inr]=Tumor_max
    New_DF_BBE_NT=cbind(New_DF_BBE_NT,inputdata1[,i])
    Col_Name_BBE_NT[BBE_NT_inr]<-names(inputdata1[i])
    
    N_min_BBE[BBE_inr]=Normal_min
    N_max_BBE[BBE_inr]=Normal_max
    T_min_BBE[BBE_inr]=Tumor_min
    T_max_BBE[BBE_inr]=Tumor_max
    New_DF_BBE=cbind(New_DF_BBE,inputdata1[,i])
    Col_Name_BBE[BBE_inr]<-names(inputdata1[i])
    
    o=o+1
    inr=inr+1
    BBE_inr=BBE_inr+1
    BBE_NT_inr=BBE_NT_inr+1
  }
  #Case 4b)Exact Bottom Boundary & Complete Embedding Normal inside Tumor 
  else if(Tumor_max==Normal_max && Normal_min<Tumor_min){
    N_TN_min_BBE[BBE_TN_inr]=Normal_min
    N_TN_max_BBE[BBE_TN_inr]=Normal_max
    T_TN_min_BBE[BBE_TN_inr]=Tumor_min
    T_TN_max_BBE[BBE_TN_inr]=Tumor_max
    New_DF_BBE_TN=cbind(New_DF_BBE_TN,inputdata1[,i])
    Col_Name_BBE_TN[BBE_TN_inr]<-names(inputdata1[i])
    
    N_min_BBE[BBE_inr]=Normal_min
    N_max_BBE[BBE_inr]=Normal_max
    T_min_BBE[BBE_inr]=Tumor_min
    T_max_BBE[BBE_inr]=Tumor_max
    New_DF_BBE=cbind(New_DF_BBE,inputdata1[,i])
    Col_Name_BBE[BBE_inr]<-names(inputdata1[i])
    
    p=p+1
    inr=inr+1
    BBE_inr=BBE_inr+1
    BBE_TN_inr=BBE_TN_inr+1
  }
  
  #Case 5a) Partial Overlap  Normal and Tumor range
  else if(Normal_min<Tumor_min && Normal_min<Tumor_max && Normal_max>Tumor_min && Normal_max<Tumor_max){
    #print(paste("i=",i,"@",Gene_ID[i]))
    N_NT_min_PO[PO_NT_inr]=Normal_min
    N_NT_max_PO[PO_NT_inr]=Normal_max
    T_NT_min_PO[PO_NT_inr]=Tumor_min
    T_NT_max_PO[PO_NT_inr]=Tumor_max
    New_DF_PO_NT=cbind(New_DF_PO_NT,inputdata1[,i])
    Col_Name_PO_NT[PO_NT_inr]<-names(inputdata1[i])
    
    N_min_PO[PO_inr]=Normal_min
    N_max_PO[PO_inr]=Normal_max
    T_min_PO[PO_inr]=Tumor_min
    T_max_PO[PO_inr]=Tumor_max
    New_DF_PO=cbind(New_DF_PO,inputdata1[,i])
    Col_Name_PO[PO_inr]<-names(inputdata1[i])
    
    partial1=partial1+1
    inr=inr+1
    PO_inr=PO_inr+1
    PO_NT_inr=PO_NT_inr+1
  }
  
  #Case 5b) Partial Overlap Tumor and Normal range
  else if(Normal_min>Tumor_min && Normal_min<Tumor_max && Normal_max>Tumor_min && Normal_max>Tumor_max){ 
    N_TN_min_PO[PO_TN_inr]=Normal_min
    N_TN_max_PO[PO_TN_inr]=Normal_max
    T_TN_min_PO[PO_TN_inr]=Tumor_min
    T_TN_max_PO[PO_TN_inr]=Tumor_max
    New_DF_PO_TN=cbind(New_DF_PO_TN,inputdata1[,i])
    Col_Name_PO_TN[PO_TN_inr]<-names(inputdata1[i])
    
    N_min_PO[PO_inr]=Normal_min
    N_max_PO[PO_inr]=Normal_max
    T_min_PO[PO_inr]=Tumor_min
    T_max_PO[PO_inr]=Tumor_max
    New_DF_PO=cbind(New_DF_PO,inputdata1[,i])
    Col_Name_PO[PO_inr]<-names(inputdata1[i])
    
    partial2=partial2+1
    inr=inr+1
    PO_inr=PO_inr+1
    PO_TN_inr=PO_TN_inr+1
  }
  #Case 6a)Boundary Overlap  Normal and Tumor range
  else if(Normal_max==Tumor_min && Tumor_max>Normal_max){
    N_NT_min_TBBO[TBBO_NT_inr]<-Normal_min
    N_NT_max_TBBO[TBBO_NT_inr]<-Normal_max
    T_NT_min_TBBO[TBBO_NT_inr]<-Tumor_min
    T_NT_max_TBBO[TBBO_NT_inr]<-Tumor_max
    New_DF_TBBO_NT=cbind(New_DF_TBBO_NT,inputdata1[,i])
    Col_Name_TBBO_NT[TBBO_NT_inr]<-names(inputdata1[i])
    
    N_min_TBBO[TBBO_inr]=Normal_min
    N_max_TBBO[TBBO_inr]=Normal_max
    T_min_TBBO[TBBO_inr]=Tumor_min
    T_max_TBBO[TBBO_inr]=Tumor_max
    New_DF_TBBO=cbind(New_DF_TBBO,inputdata1[,i])
    Col_Name_TBBO[TBBO_inr]<-names(inputdata1[i])
    
    q=q+1
    inr=inr+1
    TBBO_inr=TBBO_inr+1
    TBBO_NT_inr=TBBO_NT_inr+1
    
  }
  #Case 6b)Boundary Overlap Tumor  and Normal range
  else if(Tumor_max==Normal_min && Normal_max>Tumor_max){
    N_TN_min_TBBO[TBBO_TN_inr]<-Normal_min
    N_TN_max_TBBO[TBBO_TN_inr]<-Normal_max
    T_TN_min_TBBO[TBBO_TN_inr]<-Tumor_min
    T_TN_max_TBBO[TBBO_TN_inr]<-Tumor_max
    New_DF_TBBO_TN=cbind(New_DF_TBBO_TN,inputdata1[,i])
    Col_Name_TBBO_TN[TBBO_TN_inr]<-names(inputdata1[i])
    
    N_min_TBBO[TBBO_inr]=Normal_min
    N_max_TBBO[TBBO_inr]=Normal_max
    T_min_TBBO[TBBO_inr]=Tumor_min
    T_max_TBBO[TBBO_inr]=Tumor_max
    New_DF_TBBO=cbind(New_DF_TBBO,inputdata1[,i])
    Col_Name_TBBO[TBBO_inr]<-names(inputdata1[i])
    
    r=r+1
    inr=inr+1
    TBBO_inr=TBBO_inr+1
    TBBO_TN_inr=TBBO_TN_inr+1
  }
  
  #Case 7a)Complete Separation  Normal and Tumor range
  else if(Tumor_min>Normal_max && Tumor_max>Normal_max){   # N and T range
    N_NT_min_CS[CS_NT_inr]<-Normal_min
    N_NT_max_CS[CS_NT_inr]<-Normal_max
    T_NT_min_CS[CS_NT_inr]<-Tumor_min
    T_NT_max_CS[CS_NT_inr]<-Tumor_max
    New_DF_CS_NT=cbind(New_DF_CS_NT,inputdata1[,i])
    Col_Name_CS_NT[CS_NT_inr]<-names(inputdata1[i])
    
    N_min_CS[CS_inr]=Normal_min
    N_max_CS[CS_inr]=Normal_max
    T_min_CS[CS_inr]=Tumor_min
    T_max_CS[CS_inr]=Tumor_max
     New_DF_CS=cbind(New_DF_CS,inputdata1[,i])
    Col_Name_CS[CS_inr]<-names(inputdata1[i])
    
    complete1=complete1+1
    inr=inr+1
    CS_inr=CS_inr+1
    CS_NT_inr=CS_NT_inr+1
  } 
  #Case 7b)Complete Separation  Tumor and Normal range
  else if(Normal_min>Tumor_max && Normal_max>Tumor_max){    # T and N range
    N_TN_min_CS[CS_TN_inr]<-Normal_min
    N_TN_max_CS[CS_TN_inr]<-Normal_max
    T_TN_min_CS[CS_TN_inr]<-Tumor_min
    T_TN_max_CS[CS_TN_inr]<-Tumor_max
    New_DF_CS_TN=cbind(New_DF_CS_TN,inputdata1[,i])
    Col_Name_CS_TN[CS_TN_inr]<-names(inputdata1[i]) 
    
    N_min_CS[CS_inr]=Normal_min
    N_max_CS[CS_inr]=Normal_max
    T_min_CS[CS_inr]=Tumor_min
    T_max_CS[CS_inr]=Tumor_max
    New_DF_CS=cbind(New_DF_CS,inputdata1[,i])
    Col_Name_CS[CS_inr]<-names(inputdata1[i])
    
    complete2=complete2+1
    inr=inr+1
    CS_inr=CS_inr+1
    CS_TN_inr=CS_TN_inr+1
  }
  
  else{
    s=s+1
  }
}#End of for Dataset with n features
end_time1 <- Sys.time()
inner_time=end_time1-start_time1
inner_time


i
#Displaying the 13 categories of Features
print(paste("Length of New_DF_EO=",length(New_DF_EO)))

print(paste("Length of New_DF_CE=",length(New_DF_CE)))
print(paste("Length of New_DF_CE_NT=",length(New_DF_CE_NT)))
print(paste("Length of New_DF_CE_TN=",length(New_DF_CE_TN)))

print(paste("Length of New_DF_TBE=",length(New_DF_TBE)))
print(paste("Length of New_DF_TBE_NT=",length(New_DF_TBE_NT)))
print(paste("Length of New_DF_TBE_TN=",length(New_DF_TBE_TN)))

print(paste("Length of New_DF_BBE=",length(New_DF_BBE)))
print(paste("Length of New_DF_BBE_NT=",length(New_DF_BBE_NT)))
print(paste("Length of New_DF_BBE_TN=",length(New_DF_BBE_TN)))

print(paste("Length of New DF=",length(New_DF_PO)))
print(paste("Length of New_DF_PO_TN=",length(New_DF_PO_NT)))
print(paste("Length of New_DF_PO_TN=",length(New_DF_PO_TN)))

print(paste("Length of New_DF_TBBO=",length(New_DF_TBBO)))
print(paste("Length of New_DF_TBBO_TN=",length(New_DF_TBBO_NT)))
print(paste("Length of New_DF_TBBO_TN=",length(New_DF_TBBO_TN)))

print(paste("Length of New_DF_Separation=",length(New_DF_CS)))
print(paste("Length of New_DF_CS_NT=",length(New_DF_CS_NT)))
print(paste("Length of New_DF_CS_TN=",length(New_DF_CS_TN)))

print(paste("Eliminated genes with Exact Overlap =",j-1))
print(paste("Eliminated genes with Complete Embed N & T =",k-1))
print(paste("Eliminated genes with Complete Embed T & N =",l-1))
print(paste("Eliminated genes with Top Boundary and Embedding N & T=",m-1))
print(paste("Eliminated genes with Top Boundary and Embedding T & N=",n-1))
print(paste("Eliminated genes with Bottom Boundary and Embedding N & T =",o-1))
print(paste("Eliminated genes with Bottom Boundary and Embedding T & N =",p-1))

print(paste("Normal and Tumor Range partial Overlapping =",partial1-1))
print(paste("Tumor and Normal Range partial Overlapping =",partial2-1))
print(paste("Boundary Overlap Normal and Tumor Range =",q-1))
print(paste("Boundary Overlap Tumor and Normal Range =",r-1))
print(paste("Normal and Tumor Range complete separation =",complete1-1))
print(paste("Tumor and Normal Range complete separation =",complete2-1))
print(paste("Conditionless Else=",s-1))

Gene_Num_Summary<-data.frame(GENECOUNT=cols,EO=j-1,CE_NT=k-1,CE_TN=l-1,TOE_NT=m-1,TOE_TN=n-1,BOE_NT=o-1,BOE_TN=p-1,
                  PO_NT=partial1-1,PO_TN=partial2-1,TBBO_NT=q-1,TBBO_TN=r-1,CS_NT=complete1-1,CS_TN=complete2-1)
Gene_Num_Summary
write.csv(Gene_Num_Summary,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/GDS4275_Gene_Summary_DF.csv",row.names = F)

DF_Summary<-data.frame(New_DF_EO=length(New_DF_EO), New_DF_CE=length(New_DF_CE), New_DF_TBE=length(New_DF_TBE),
                       New_DF_BBE=length(New_DF_BBE), New_DF_PO=length(New_DF_PO),New_DF_TBBO=length(New_DF_TBBO),
                       New_DF_CS=length(New_DF_CS))
DF_Summary
write.csv(DF_Summary,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/GDS4275_DF_Summary.csv",row.names = F)

DF_Summary_Sub<-data.frame(New_DF_EO=length(New_DF_EO), New_DF_CE=length(New_DF_CE), New_DF_CE_NT=length(New_DF_CE_NT), New_DF_CE_TN=length(New_DF_CE_TN),
                       New_DF_TBE=length(New_DF_TBE), New_DF_TBE_NT=length(New_DF_TBE_NT), New_DF_TBE_TN=length(New_DF_TBE_TN),
                       New_DF_BBE=length(New_DF_BBE), New_DF_BBE_NT=length(New_DF_BBE_NT), New_DF_BBE_TN=length(New_DF_BBE_TN),
                       New_DF_PO=length(New_DF_PO), New_DF_PO_NT=length(New_DF_PO_NT), New_DF_PO_TN=length(New_DF_PO_TN),
                       New_DF_TBBO=length(New_DF_TBBO), New_DF_TBBO_NT=length(New_DF_TBBO_NT), New_DF_TBBO_TN=length(New_DF_TBBO_TN),
                       New_DF_CS=length(New_DF_CS), New_DF_CS_NT=length(New_DF_CS_NT), New_DF_CS_TN=length(New_DF_CS_TN))
DF_Summary_Sub
write.csv(DF_Summary_Sub,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/GDS4275_DF_Summary_Sub_Summary.csv",row.names = F)

inr
################################################
## CASE I)     PARTIAL OVERLAPPING
################################################
PO_DF_Summary<-data.frame(New_DF_PO=length(New_DF_PO),New_DF_PO_NT=length(New_DF_PO_NT),
                          New_DF_PO_TN=length(New_DF_PO_TN))
PO_DF_Summary
write.csv(PO_DF_Summary,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/PO/GDS4275_PO_DF_Summary.csv",row.names = F)

##########################################################
#             PO_MinMax_Vectors_NT and TN
##########################################################
po1<-data.frame()
po1<-rbind(po1,N_NT_min_PO)
po1<-rbind(po1,N_NT_max_PO)
po1<-rbind(po1,T_NT_min_PO)
po1<-rbind(po1,T_NT_max_PO)
names(po1)[] <- Col_Name_PO_NT
po1<-as.data.frame(sapply(po1,as.numeric))
write.csv(po1,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/PO/GDS4275_PO_MinMax_Vectors_NT.csv",row.names = F)
length(po1)

po2<-data.frame()
po2<-rbind(po2,N_TN_min_PO)
po2<-rbind(po2,N_TN_max_PO)
po2<-rbind(po2,T_TN_min_PO)
po2<-rbind(po2,T_TN_max_PO)
names(po2)[] <- Col_Name_PO_TN
po2<-as.data.frame(sapply(po2,as.numeric))
write.csv(po2,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/PO/GDS4275_PO_MinMax_Vectors_TN.csv",row.names = F)
length(po2)

po<-data.frame()
po<-rbind(po,N_min_PO)
po<-rbind(po,N_max_PO)
po<-rbind(po,T_min_PO)
po<-rbind(po,T_max_PO)
names(po)[] <- Col_Name_PO
po<-as.data.frame(sapply(po,as.numeric))
write.csv(po,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/PO/GDS4275_PO_MinMax_Vectors.csv",row.names = F)
length(po)

#Col_Name_Vector
DF_Col_PO_len<-length(Col_Name_PO)
DF_Col_PO_len
DF_Col_PO_NT_len<-length(Col_Name_PO_NT)
DF_Col_PO_NT_len
DF_Col_PO_TN_len<-length(Col_Name_PO_TN)
DF_Col_PO_TN_len

# Changing col names
New_DF_PO<-round(New_DF_PO,digits=1)
New_DF_PO<-as.data.frame(sapply(New_DF_PO,as.numeric))
#print(sapply(New_DF_PO, class))
is.numeric(New_DF_PO)
#names(df)[names(df) == "old.name"] <- "new.name"
names(New_DF_PO)[] <- Col_Name_PO

New_DF_PO_NT<-round(New_DF_PO_NT,digits=1)
New_DF_PO_NT<-as.data.frame(sapply(New_DF_PO_NT,as.numeric))
names(New_DF_PO_NT)[] <- Col_Name_PO_NT

New_DF_PO_TN<-round(New_DF_PO_TN,digits=1)
New_DF_PO_TN<-as.data.frame(sapply(New_DF_PO_TN,as.numeric))
names(New_DF_PO_TN)[] <- Col_Name_PO_TN


##########################################################
#             LUMO_MinMax_Vectors_NT and TN
##########################################################
tbbo1<-data.frame()
tbbo1<-rbind(tbbo1,N_NT_min_TBBO)
tbbo1<-rbind(tbbo1,N_NT_max_TBBO)
tbbo1<-rbind(tbbo1,T_NT_min_TBBO)
tbbo1<-rbind(tbbo1,T_NT_max_TBBO)
names(tbbo1)[] <- Col_Name_TBBO_NT
tbbo1<-as.data.frame(sapply(tbbo1,as.numeric))
write.csv(tbbo1,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/LUMO/GDS4275_LUMO_MinMax_Vectors_NT.csv",row.names = F)
length(tbbo1)

tbbo2<-data.frame()
tbbo2<-rbind(tbbo2,N_TN_min_TBBO)
tbbo2<-rbind(tbbo2,N_TN_max_TBBO)
tbbo2<-rbind(tbbo2,T_TN_min_TBBO)
tbbo2<-rbind(tbbo2,T_TN_max_TBBO)
names(tbbo2)[] <- Col_Name_TBBO_TN
tbbo2<-as.data.frame(sapply(tbbo2,as.numeric))
write.csv(tbbo2,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/LUMO/GDS4275_LUMO_MinMax_Vectors_TN.csv",row.names = F)
length(tbbo2)

tbbo<-data.frame()
tbbo<-rbind(tbbo,N_min_TBBO)
tbbo<-rbind(tbbo,N_max_TBBO)
tbbo<-rbind(tbbo,T_min_TBBO)
tbbo<-rbind(tbbo,T_max_TBBO)
names(tbbo)[] <- Col_Name_TBBO
tbbo<-as.data.frame(sapply(tbbo,as.numeric))
write.csv(tbbo,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/LUMO/GDS4275_LUMO_MinMax_Vectors.csv",row.names = F)
length(tbbo)

#Col_Name_Vector
DF_Col_TBBO_len<-length(Col_Name_TBBO)
DF_Col_TBBO_len
DF_Col_TBBO_NT_len<-length(Col_Name_TBBO_NT)
DF_Col_TBBO_NT_len
DF_Col_TBBO_TN_len<-length(Col_Name_TBBO_TN)
DF_Col_TBBO_TN_len

# Changing col names
New_DF_TBBO<-round(New_DF_TBBO,digits=1)
New_DF_TBBO<-as.data.frame(sapply(New_DF_TBBO,as.numeric))
#print(sapply(New_DF_TBBO, class))
is.numeric(New_DF_TBBO)
#names(df)[names(df) == "old.name"] <- "new.name"
names(New_DF_TBBO)[] <- Col_Name_TBBO

New_DF_TBBO_NT<-round(New_DF_TBBO_NT,digits=1)
New_DF_TBBO_NT<-as.data.frame(sapply(New_DF_TBBO_NT,as.numeric))
names(New_DF_TBBO_NT)[] <- Col_Name_TBBO_NT

New_DF_TBBO_TN<-round(New_DF_TBBO_TN,digits=1)
New_DF_TBBO_TN<-as.data.frame(sapply(New_DF_TBBO_TN,as.numeric))
names(New_DF_TBBO_TN)[] <- Col_Name_TBBO_TN


##########################################################
#             CS_MinMax_Vectors_NT and TN
##########################################################
cs1<-data.frame()
cs1<-rbind(cs1,N_NT_min_CS)
cs1<-rbind(cs1,N_NT_max_CS)
cs1<-rbind(cs1,T_NT_min_CS)
cs1<-rbind(cs1,T_NT_max_CS)
names(cs1)[] <- Col_Name_CS_NT
cs1<-as.data.frame(sapply(cs1,as.numeric))
write.csv(cs1,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/CS/GDS4275_CS_MinMax_Vectors_NT.csv",row.names = F)
length(cs1)

cs2<-data.frame()
cs2<-rbind(cs2,N_TN_min_CS)
cs2<-rbind(cs2,N_TN_max_CS)
cs2<-rbind(cs2,T_TN_min_CS)
cs2<-rbind(cs2,T_TN_max_CS)
names(cs2)[] <- Col_Name_CS_TN
cs2<-as.data.frame(sapply(cs2,as.numeric))
write.csv(cs2,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/CS/GDS4275_CS_MinMax_Vectors_TN.csv",row.names = F)
length(cs2)

cs<-data.frame()
cs<-rbind(cs,N_min_CS)
cs<-rbind(cs,N_max_CS)
cs<-rbind(cs,T_min_CS)
cs<-rbind(cs,T_max_CS)
names(cs)[] <- Col_Name_CS
cs<-as.data.frame(sapply(cs,as.numeric))
write.csv(cs,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/CS/GDS4275_CS_MinMax_Vectors.csv",row.names = F)
length(cs)

#Col_Name_Vector
DF_Col_CS_len<-length(Col_Name_CS)
DF_Col_CS_len
DF_Col_CS_NT_len<-length(Col_Name_CS_NT)
DF_Col_CS_NT_len
DF_Col_CS_TN_len<-length(Col_Name_CS_TN)
DF_Col_CS_TN_len

# Changing col names
New_DF_CS<-round(New_DF_CS,digits=1)
New_DF_CS<-as.data.frame(sapply(New_DF_CS,as.numeric))
#print(sapply(New_DF_CS, class))
is.numeric(New_DF_CS)
#names(df)[names(df) == "old.name"] <- "new.name"
names(New_DF_CS)[] <- Col_Name_CS

New_DF_CS_NT<-round(New_DF_CS_NT,digits=1)
New_DF_CS_NT<-as.data.frame(sapply(New_DF_CS_NT,as.numeric))
names(New_DF_CS_NT)[] <- Col_Name_CS_NT

New_DF_CS_TN<-round(New_DF_CS_TN,digits=1)
New_DF_CS_TN<-as.data.frame(sapply(New_DF_CS_TN,as.numeric))
names(New_DF_CS_TN)[] <- Col_Name_CS_TN


###########################################################
# NON INFORMATIVE PATTERNS
############################################################
##########################################################
#             EO_MinMax_Vectors
##########################################################
eo<-data.frame()
eo<-rbind(eo,N_min_EO)
eo<-rbind(eo,N_max_EO)
eo<-rbind(eo,T_min_EO)
eo<-rbind(eo,T_max_EO)
names(eo)[] <- Col_Name_EO
eo<-as.data.frame(sapply(eo,as.numeric))
write.csv(eo,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/EO/GDS4275_EO_MinMax_Vectors.csv",row.names = F)
length(eo)

#Col_Name_Vector
DF_Col_EO_len<-length(Col_Name_EO)
DF_Col_EO_len
# Changing col names
#EO
New_DF_EO<-round(New_DF_EO,digits=1)
New_DF_EO<-as.data.frame(sapply(New_DF_EO,as.numeric))
names(New_DF_EO)[] <- Col_Name_EO

##########################################################
#             CE_MinMax_Vectors
##########################################################
ce1<-data.frame()
ce1<-rbind(ce1,N_NT_min_CE)
ce1<-rbind(ce1,N_NT_max_CE)
ce1<-rbind(ce1,T_NT_min_CE)
ce1<-rbind(ce1,T_NT_max_CE)
names(ce1)[] <- Col_Name_CE_NT
ce1<-as.data.frame(sapply(ce1,as.numeric))
write.csv(ce1,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/CE/GDS4275_CE_MinMax_Vectors_NT.csv",row.names = F)
length(ce1)

ce2<-data.frame()
ce2<-rbind(ce2,N_TN_min_CE)
ce2<-rbind(ce2,N_TN_max_CE)
ce2<-rbind(ce2,T_TN_min_CE)
ce2<-rbind(ce2,T_TN_max_CE)
names(ce2)[] <- Col_Name_CE_TN
ce2<-as.data.frame(sapply(ce2,as.numeric))
write.csv(ce2,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/CE/GDS4275_CE_MinMax_Vectors_TN.csv",row.names = F)
length(ce2)

ce<-data.frame()
ce<-rbind(ce,N_min_CE)
ce<-rbind(ce,N_max_CE)
ce<-rbind(ce,T_min_CE)
ce<-rbind(ce,T_max_CE)
names(ce)[] <- Col_Name_CE
ce<-as.data.frame(sapply(ce,as.numeric))
write.csv(ce,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/CE/GDS4275_CE_MinMax_Vectors.csv",row.names = F)
length(ce)

#Col_Name_Vector
DF_Col_CE_len<-length(Col_Name_CE)
DF_Col_CE_len
DF_Col_CE_NT_len<-length(Col_Name_CE_NT)
DF_Col_CE_NT_len
DF_Col_CE_TN_len<-length(Col_Name_CE_TN)
DF_Col_CE_TN_len
New_DF_CE<-round(New_DF_CE,digits=1)
New_DF_CE<-as.data.frame(sapply(New_DF_CE,as.numeric))
names(New_DF_CE)[] <- Col_Name_CE

New_DF_CE_NT<-round(New_DF_CE_NT,digits=1)
New_DF_CE_NT<-as.data.frame(sapply(New_DF_CE_NT,as.numeric))
names(New_DF_CE_NT)[] <- Col_Name_CE_NT

New_DF_CE_TN<-round(New_DF_CE_TN,digits=1)
New_DF_CE_TN<-as.data.frame(sapply(New_DF_CE_TN,as.numeric))
names(New_DF_CE_TN)[] <- Col_Name_CE_TN
##########################################################
#             TBE_MinMax_Vectors
##########################################################
tbe1<-data.frame()
tbe1<-rbind(tbe1,N_NT_min_TBE)
tbe1<-rbind(tbe1,N_NT_max_TBE)
tbe1<-rbind(tbe1,T_NT_min_TBE)
tbe1<-rbind(tbe1,T_NT_max_TBE)
names(tbe1)[] <- Col_Name_TBE_NT
tbe1<-as.data.frame(sapply(tbe1,as.numeric))
write.csv(tbe1,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/LMO-UME/GDS4275_LMO-UME_MinMax_Vectors_NT.csv",row.names = F)
length(tbe1)

tbe2<-data.frame()
tbe2<-rbind(tbe2,N_TN_min_TBE)
tbe2<-rbind(tbe2,N_TN_max_TBE)
tbe2<-rbind(tbe2,T_TN_min_TBE)
tbe2<-rbind(tbe2,T_TN_max_TBE)
names(tbe2)[] <- Col_Name_TBE_TN
tbe2<-as.data.frame(sapply(tbe2,as.numeric))
write.csv(tbe2,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/LMO-UME/GDS4275_LMO-UME_MinMax_Vectors_TN.csv",row.names = F)
length(tbe2)

tbe<-data.frame()
tbe<-rbind(tbe,N_min_TBE)
tbe<-rbind(tbe,N_max_TBE)
tbe<-rbind(tbe,T_min_TBE)
tbe<-rbind(tbe,T_max_TBE)
names(tbe)[] <- Col_Name_TBE
tbe<-as.data.frame(sapply(tbe,as.numeric))
write.csv(tbe,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/LMO-UME/GDS4275_LMO-UME_MinMax_Vectors.csv",row.names = F)
length(tbe)

#Col_Name_Vector
DF_Col_TBE_len<-length(Col_Name_TBE)
DF_Col_TBE_len
DF_Col_TBE_NT_len<-length(Col_Name_TBE_NT)
DF_Col_TBE_NT_len
DF_Col_TBE_TN_len<-length(Col_Name_TBE_TN)
DF_Col_TBE_TN_len
#TBE
New_DF_TBE<-round(New_DF_TBE,digits=1)
New_DF_TBE<-as.data.frame(sapply(New_DF_TBE,as.numeric))
names(New_DF_TBE)[] <- Col_Name_TBE

New_DF_TBE_NT<-round(New_DF_TBE_NT,digits=1)
New_DF_TBE_NT<-as.data.frame(sapply(New_DF_TBE_NT,as.numeric))
names(New_DF_TBE_NT)[] <- Col_Name_TBE_NT

New_DF_TBE_TN<-round(New_DF_TBE_TN,digits=1)
New_DF_TBE_TN<-as.data.frame(sapply(New_DF_TBE_TN,as.numeric))
names(New_DF_TBE_TN)[] <- Col_Name_TBE_TN
##########################################################
#             BBE_MinMax_Vectors
##########################################################
#TBE
#N_NT_min_BBE
#N_NT_max_BBE
#T_NT_min_BBE
#T_NT_max_BBE
bbe1<-data.frame()
bbe1<-rbind(bbe1,N_NT_min_BBE)
bbe1<-rbind(bbe1,N_NT_max_BBE)
bbe1<-rbind(bbe1,T_NT_min_BBE)
bbe1<-rbind(bbe1,T_NT_max_BBE)
names(bbe1)[] <- Col_Name_BBE_NT
bbe1<-as.data.frame(sapply(bbe1,as.numeric))
write.csv(bbe1,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/UMO-LME/GDS4275_UMO-LME_MinMax_Vectors_NT.csv",row.names = F)
length(bbe1)

#N_TN_min_BBE
#N_TN_max_BBE
#T_TN_min_BBE
#T_TN_max_BBE
bbe2<-data.frame()
bbe2<-rbind(bbe2,N_TN_min_BBE)
bbe2<-rbind(bbe2,N_TN_max_BBE)
bbe2<-rbind(bbe2,T_TN_min_BBE)
bbe2<-rbind(bbe2,T_TN_max_BBE)
names(bbe2)[] <- Col_Name_BBE_TN
bbe2<-as.data.frame(sapply(bbe2,as.numeric))
write.csv(bbe2,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/UMO-LME/GDS4275_UMO-LME_MinMax_Vectors_TN.csv",row.names = F)
length(bbe2)

#N_min_BBE
#N_max_BBE
#T_min_BBE
#T_max_BBE
bbe<-data.frame()
bbe<-rbind(bbe,N_min_BBE)
bbe<-rbind(bbe,N_max_BBE)
bbe<-rbind(bbe,T_min_BBE)
bbe<-rbind(bbe,T_max_BBE)
names(bbe)[] <- Col_Name_BBE
bbe<-as.data.frame(sapply(bbe,as.numeric))
write.csv(bbe,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/UMO-LME/GDS4275_UMO-LME_MinMax_Vectors_TN.csv",row.names = F)
length(bbe)

#Col_Name_Vector
DF_Col_BBE_len<-length(Col_Name_BBE)
DF_Col_BBE_len
DF_Col_BBE_NT_len<-length(Col_Name_BBE_NT)
DF_Col_BBE_NT_len
DF_Col_BBE_TN_len<-length(Col_Name_BBE_TN)
DF_Col_BBE_TN_len
New_DF_BBE<-round(New_DF_BBE,digits=1)
New_DF_BBE<-as.data.frame(sapply(New_DF_BBE,as.numeric))
names(New_DF_BBE)[] <- Col_Name_BBE

New_DF_BBE_NT<-round(New_DF_BBE_NT,digits=1)
New_DF_BBE_NT<-as.data.frame(sapply(New_DF_BBE_NT,as.numeric))
names(New_DF_BBE_NT)[] <- Col_Name_BBE_NT

New_DF_BBE_TN<-round(New_DF_BBE_TN,digits=1)
New_DF_BBE_TN<-as.data.frame(sapply(New_DF_BBE_TN,as.numeric))
names(New_DF_BBE_TN)[] <- Col_Name_BBE_TN

write.csv(New_DF_EO,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/EO/GDS4275_New_DF_EO.csv",row.names = F)

write.csv(New_DF_CE,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/CE/GDS4275_New_DF_CE.csv",row.names = F)
write.csv(New_DF_CE_NT,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/CE/GDS4275_New_DF_CE_NT.csv",row.names = F)
write.csv(New_DF_CE_TN,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/CE/GDS4275_New_DF_CE_TN.csv",row.names = F)

write.csv(New_DF_TBE,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/LMO-UME/GDS4275_New_DF_TBE.csv",row.names = F)
write.csv(New_DF_TBE_NT,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/LMO-UME/GDS4275_New_DF_TBE_NT.csv",row.names = F)
write.csv(New_DF_TBE_TN,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/LMO-UME/GDS4275_New_DF_TBE_TN.csv",row.names = F)

write.csv(New_DF_BBE,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/UMO-LME/GDS4275_New_DF_BBE.csv",row.names = F)
write.csv(New_DF_BBE_NT,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/UMO-LME/GDS4275_New_DF_BBE_NT.csv",row.names = F)
write.csv(New_DF_BBE_TN,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/UMO-LME/GDS4275_New_DF_BBE_TN.csv",row.names = F)

############################################################
#Writing Partial Overlap DF to a CSV File                     
write.csv(New_DF_PO,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/PO/GDS4275_New_DF_PO.csv",row.names = F)
write.csv(New_DF_PO_NT,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/PO/GDS4275_New_DF_PO_NT.csv",row.names = F)
write.csv(New_DF_PO_TN,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/PO/GDS4275_New_DF_PO_TN.csv",row.names = F)
length(New_DF_PO)
length(New_DF_PO_NT)
length(New_DF_PO_TN)

#Writing LUMO to a CSV File                     
write.csv(New_DF_TBBO,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/LUMO/GDS4275_New_DF_LUMO.csv",row.names = F)
write.csv(New_DF_TBBO_NT,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/LUMO/GDS4275_New_DF_LUMO_NT.csv",row.names = F)
write.csv(New_DF_TBBO_TN,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/LUMO/GDS4275_New_DF_LUMO_TN.csv",row.names = F)
length(New_DF_TBBO)
length(New_DF_TBBO_NT)
length(New_DF_TBBO_TN)

#Writing CS DF to a CSV File                     
write.csv(New_DF_CS,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/CS/GDS4275_New_DF_CS.csv",row.names = F)
write.csv(New_DF_CS_NT,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/CS/GDS4275_New_DF_CS_NT.csv",row.names = F)
write.csv(New_DF_CS_TN,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/CS/GDS4275_New_DF_CS_TN.csv",row.names = F)
length(New_DF_CS)
length(New_DF_CS_NT)
length(New_DF_CS_TN)

#######
############ binding Dataframes and Min-Max Vectors
a<-data.frame()
a<-rbind(a,New_DF_EO)
a<-rbind(a,eo)

b<-data.frame()
b<-rbind(b,New_DF_CE)
b<-rbind(b,ce)

c<-data.frame()
c<-rbind(c,New_DF_CE_NT)
c<-rbind(c,ce1)

d<-data.frame()
d<-rbind(d,New_DF_CE_TN)
d<-rbind(d,ce2)                     

e<-data.frame()
e<-rbind(e,New_DF_TBE)
e<-rbind(e,tbe)

f<-data.frame()
f<-rbind(f,New_DF_TBE_NT)
f<-rbind(f,tbe1)

g<-data.frame()
g<-rbind(g,New_DF_TBE_TN)
g<-rbind(g,tbe2)
    
h<-data.frame()                  
h<-rbind(h,New_DF_BBE)
h<-rbind(h,bbe)

j<-data.frame()
j<-rbind(j,New_DF_BBE_NT)
j<-rbind(j,bbe1)

k<-data.frame()
k<-rbind(k,New_DF_BBE_TN)
k<-rbind(k,bbe2)   

l<-data.frame()               
l<-rbind(l,New_DF_PO)
l<-rbind(l,po)

m<-data.frame()
m<-rbind(m,New_DF_PO_NT)
m<-rbind(m,po1)

n<-data.frame()
n<-rbind(n,New_DF_PO_TN)
n<-rbind(n,po2)  

o<-data.frame()
o<-rbind(o,New_DF_TBBO)
o<-rbind(o,tbbo)

p<-data.frame()
p<-rbind(p,New_DF_TBBO_NT)
p<-rbind(p,tbbo1)

q<-data.frame()
q<-rbind(q,New_DF_TBBO_TN)
q<-rbind(q,tbbo2)

r<-data.frame()
r<-rbind(r,New_DF_CS)
r<-rbind(r,cs)

s<-data.frame()
s<-rbind(s,New_DF_CS_NT)
s<-rbind(s,cs1)

t<-data.frame()
t<-rbind(t,New_DF_CS_TN)
t<-rbind(t,cs2)

write.csv(a,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/EO/GDS4275_EO_DF_Min-Max.csv",row.names = F)

write.csv(b,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/CE/GDS4275_CE_DF_Min-Max.csv",row.names = F)
write.csv(c,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/CE/GDS4275_CE_NT_DF_Min-Max.csv",row.names = F)
write.csv(d,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/CE/GDS4275_CE_TN_DF_Min-Max.csv",row.names = F)

write.csv(e,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/LMO-UME/GDS4275_TBE_DF_Min-Max.csv",row.names = F)
write.csv(f,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/LMO-UME/GDS4275_TBE_NT_DF_Min-Max.csv",row.names = F)
write.csv(g,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/LMO-UME/GDS4275_TBE_TN_DF_Min-Max.csv",row.names = F)

write.csv(h,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/UMO-LME/GDS4275_BBE_DF_Min-Max.csv",row.names = F)
write.csv(j,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/UMO-LME/GDS4275_BBE_NT_DF_Min-Max.csv",row.names = F)
write.csv(k,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Insignificant-Patterns/UMO-LME/GDS4275_BBE_TN_DF_Min-Max.csv",row.names = F)

write.csv(l,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/PO/GDS4275_PO_DF_Min-Max.csv",row.names = F)
write.csv(m,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/PO/GDS4275_PO_NT_DF_Min-Max.csv",row.names = F)
write.csv(n,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/PO/GDS4275_PO_TN_DF_Min-Max.csv",row.names = F)

write.csv(o,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/LUMO/GDS4275_LUMO_DF_Min-Max.csv",row.names = F)
write.csv(p,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/LUMO/GDS4275_LUMO_NT_DF_Min-Max.csv",row.names = F)
write.csv(q,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/LUMO/GDS4275_LUMO_TN_DF_Min-Max.csv",row.names = F)

write.csv(r,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/CS/GDS4275_CS_DF_Min-Max.csv",row.names = F)
write.csv(s,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/CS/GDS4275_CS_NT_DF_Min-Max.csv",row.names = F)
write.csv(t,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/Significant-Patterns/CS/GDS4275_CS_TN_DF_Min-Max.csv",row.names = F)


end_time <- Sys.time()
time=end_time-start_time

inner_time
time
Full_Run_Time<-data.frame(Pattern_Categorization_Time=inner_time,Full_Run_Time=time)
write.csv(Full_Run_Time,"F:/Ph.D Work/R/Coding-Towards Perfection/GDS4275-PituitaryTumor/GDS4275-PituitaryTumor_38929genes_Full_Runtime.csv",row.names = F)

#####      STOP RUNNING HERE
