#########################################
# Import Packages
#########################################
library(dplyr)
library(infotheo)
library(FSelector)

#########################################
#   Read Input Files
#########################################
inputdata<-read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/GDS3057_New_DF_EO.csv",header=TRUE)
nrow(inputdata)
ncol(inputdata)
df<-inputdata

############################
##### Read CLASS LABEL
############################
Class_lab<-read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/GDS3057-Leukemia_Label_Numeric.csv",header=TRUE)
Fact_lab<-read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/GDS3057-Leukemia_Label_Factor.csv",header=TRUE)

nrow(df)
ncol(df)
#########################################
##### cOMPUTING Pearson Coefficient
#########################################
gname<-names(inputdata)

DF_Len<-length(df)
y1 = Class_lab$Cancer.Presence
Cancer.Presence<-Fact_lab$Cancer.Presence
cor_value<-c()

start_time <- Sys.time()
for(i in 1:(DF_Len)){
  x1=df[i]
  cor_value[i] = cor(x1, y1, method = "pearson") 
}
end_time <- Sys.time()
time=end_time-start_time

time_estimate<-data.frame()
time_estimate<-data.frame(Pattern="GDS3057_EO_Pearson",Runtime_846genes=time)
write.csv(time_estimate,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Runtime.csv",row.names=F)

cor_value
gname<-names(df)
#length(gname)

Pearson_Rank_DF<-data.frame()
Pearson_Rank_DF<-data.frame(GENE_NAME=gname,Pearson_value=cor_value)
write.csv(Pearson_Rank_DF,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Cor_Rank.csv",row.names =F)

#########################################
##### Sorting Rank Vector and Dataframe
#########################################

CORR_Vector_Len<-length(cor_value)
Gene_Index_inr=1
New_DF_Gene<-data.frame(matrix(NA,nrow=64,ncol=0))
New_DF_Gene_Name<-c()

Gene_Rank_Vec<-c()
Gene_Rank_Vec<-cor_value
Gene_Rank_Vec_Name<-c()
Gene_Rank_Vec_Name<-gname
i=1
j=1

for(i in 1:(CORR_Vector_Len)){
  max_index <- i
  for(j in (i+1):CORR_Vector_Len-1){
    if(Gene_Rank_Vec[j]>Gene_Rank_Vec[max_index]){
      max_index = j
    }# End of if
  }# End of j
  
  temp=Gene_Rank_Vec[i]
  Gene_Rank_Vec[i]=Gene_Rank_Vec[max_index]
  Gene_Rank_Vec[max_index]=temp
  
  temp1=Gene_Rank_Vec_Name[i]
  Gene_Rank_Vec_Name[i]=Gene_Rank_Vec_Name[max_index]
  Gene_Rank_Vec_Name[max_index]=temp1
  
  New_DF_Gene=cbind(New_DF_Gene,inputdata[max_index])
  New_DF_Gene_Name[i]<-names(inputdata[max_index])
}# End of i
i
j


#New_DF_Gene=cbind(New_DF_Gene,inputdata[j])
names(New_DF_Gene)[] <- New_DF_Gene_Name
#New_DF_Gene
length(New_DF_Gene)
length(New_DF_Gene_Name)
############################################################################
###   Writing Rank Vectors and Full Dataframe
############################################################################
write.csv(New_DF_Gene,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_DF_Sorted.csv",row.names = F)
Pearson_Gene_Rank<-data.frame(Gene_ID=Gene_Rank_Vec_Name,Pearson_Value=Gene_Rank_Vec)
#Pearson_Gene_Rank
#nrow(Pearson_Gene_Rank)
write.csv(Pearson_Gene_Rank,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Gene_RankList_Sorted.csv",row.names = F)

############################################################################
###   For Dataframe partition work from here
############################################################################


df_50<-New_DF_Gene[,1:50]
df_50<-cbind(df_50,Cancer.Presence)
write.csv(df_50,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Top50_genes.csv",row.names =F)

df_40<-New_DF_Gene[,1:40]
df_40<-cbind(df_40,Cancer.Presence)
write.csv(df_40,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Top40_genes.csv",row.names = F)

df_30<-New_DF_Gene[,1:30]
df_30<-cbind(df_30,Cancer.Presence)
write.csv(df_30,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Top30_genes.csv",row.names = F)

df_25<-New_DF_Gene[,1:25]
df_25<-cbind(df_25,Cancer.Presence)
write.csv(df_25,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Top25_genes.csv",row.names = F)

df_20<-New_DF_Gene[,1:20]
df_20<-cbind(df_20,Cancer.Presence)
write.csv(df_20,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Top20_genes.csv",row.names = F)

df_15<-New_DF_Gene[,1:15]
df_15<-cbind(df_15,Cancer.Presence)
write.csv(df_15,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Top15_genes.csv",row.names = F)

df_10<-New_DF_Gene[,1:10]
df_10<-cbind(df_10,Cancer.Presence)
write.csv(df_10,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Top10_genes.csv",row.names = F)

df_5<-New_DF_Gene[,1:5]
df_5<-cbind(df_5,Cancer.Presence)
write.csv(df_5,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/Pearson/GDS3057-EO_Pearson_Top5_genes.csv",row.names = F)

######################################################
###########   STOP HERE   ###########################
########################################################



