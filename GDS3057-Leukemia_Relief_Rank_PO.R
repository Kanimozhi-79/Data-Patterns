#########################################
# Import Packages
#########################################
library(FSinR)

#########################################
#   Read Input Files
#########################################
inputdata<-read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/GDS3057_New_DF_PO.csv",header=TRUE)

nrow(inputdata)

gname<-c()
gname<-names(inputdata)
#length(gname)

############################
##### Read CLASS LABEL
############################
#Label needs to be numeric for calculating relief values
label<-read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/GDS3057-Leukemia_Label_Numeric.csv",header=TRUE)
inputdata<-cbind(inputdata,label)
ncol(inputdata)

#Label needs to be Factor for calculating running SVM
Fact_lab<-read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/GDS3057-Leukemia_Label_Factor.csv",header=TRUE)
Cancer.Presence<-Fact_lab$Cancer.Presence
nrow(label)
nrow(Fact_lab)


#########################################
##### cOMPUTING Relief values
#########################################
relief_evaluator <- relief()
# Evaluate the features (parameters: dataset, target variable and features)

df_len<-ncol(inputdata)
df_len
index=1
Relief_values<-c()
Gene_Col_Vector<-c()
i=1

start_time1 <- Sys.time()
for(i in 1:df_len){
Gene_Col_Vector<-inputdata[i]  
Relief_values[index]<-relief_evaluator(inputdata,'Cancer.Presence',c(Gene_Col_Vector))
index=index+1
}
end_time1 <- Sys.time()
inner_time=end_time1-start_time1
time_estimate<-data.frame()
time_estimate<-data.frame(Pattern="GDS3057_PO_Relief",Runtime_7933genes=inner_time)
write.csv(time_estimate,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Runtime.csv",row.names=F)
#Relief_values
#index
#i
#Gene_Col_Vector
Relief_values<-Relief_values[-c(7934)]

Relief_Rank_DF<-data.frame()
Relief_Rank_DF<-data.frame(Gene_Name=gname,Rank_Relief=Relief_values)
write.csv(Relief_Rank_DF,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Cor_Rank.csv",row.names = F)
end_time <- Sys.time()


########################################################################################################

#####     Sorting Rank Vector and Dataframe ---sORTING WORKS ONLY FOR NON NA VALUES 
#####     IF NA VALUES COMES AS CORRELATION VALUE DO THE RANK SORTING MANUALLY
#####     ELSE SKIP THIS SORTING PART AND GO TO LINE 120 FOR Writing Rank Vectors and Full Dataframe

########################################################################################################

cor_value<-Relief_values
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
write.csv(New_DF_Gene,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_DF_Sorted.csv",row.names = F)
Relief_Gene_Rank<-data.frame(Gene_ID=Gene_Rank_Vec_Name,Pearson_Value=Gene_Rank_Vec)
#Relief_Gene_Rank
nrow(Relief_Gene_Rank)
write.csv(Relief_Gene_Rank,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Gene_RankList_Sorted.csv",row.names = F)


############################################################################
###   For Dataframe partition work from here
############################################################################


New_DF_Gene<-read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_DF_Sorted.csv",header=TRUE)
length(New_DF_Gene)
nrow(New_DF_Gene)


df_50<-New_DF_Gene[,1:50]
df_50<-cbind(df_50,Cancer.Presence)
write.csv(df_50,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Top50_genes.csv",row.names = F)

df_40<-New_DF_Gene[,1:40]
df_40<-cbind(df_40,Cancer.Presence)
write.csv(df_40,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Top40_genes.csv",row.names = F)

df_30<-New_DF_Gene[,1:30]
df_30<-cbind(df_30,Cancer.Presence)
write.csv(df_30,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Top30_genes.csv",row.names = F)

df_25<-New_DF_Gene[,1:25]
df_25<-cbind(df_25,Cancer.Presence)
write.csv(df_25,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Top25_genes.csv",row.names = F)

df_20<-New_DF_Gene[,1:20]
df_20<-cbind(df_20,Cancer.Presence)
write.csv(df_20,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Top20_genes.csv",row.names = F)

df_15<-New_DF_Gene[,1:15]
df_15<-cbind(df_15,Cancer.Presence)
write.csv(df_15,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Top15_genes.csv",row.names = F)

df_10<-New_DF_Gene[,1:10]
df_10<-cbind(df_10,Cancer.Presence)
write.csv(df_10,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Top10_genes.csv",row.names = F)

df_5<-New_DF_Gene[,1:5]
df_5<-cbind(df_5,Cancer.Presence)
write.csv(df_5,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Relief/GDS3057-PO_Relief_Top5_genes.csv",row.names = F)


####################################################################
#       PROGRAM ENDS HERE
####################################################################


