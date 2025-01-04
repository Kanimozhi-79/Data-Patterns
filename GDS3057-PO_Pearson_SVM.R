library(kernlab)
library(e1071)
library(caret)
library(caTools)
library(class)
# Run one DataFrame at a time
# cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_Top5_genes.csv",header = TRUE)
# cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_Top10_genes.csv",header = TRUE)
# cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_Top15_genes.csv",header = TRUE)
# cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_Top20_genes.csv",header = TRUE)
# cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_Top25_genes.csv",header = TRUE)

# I.5 genes
#cancerdata
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_Top5_genes.csv",header = TRUE)

Vec_SVM_5_DF<-c()
Sno=1


for(i in 1:25){
index<-1:nrow(cancerdata)
length(index)
testindex<-sample(index,trunc(length(index)/3))
#testindex
trainset<-cancerdata[-testindex,]
#trainset
testset<-cancerdata[testindex,]
#testset
fitControl<-trainControl(method="LOOCV",number=10,classProbs=TRUE)
#fitControl
svm_model<-train(Cancer.Presence ~.,data=trainset,method='svmRadial',trControl=fitControl,tuneLength=3)
svm_model
testset$svm_pred<-predict(svm_model,testset[,])
testset$svm_pred
confusionMatrix(as.factor(testset$Cancer.Presence),as.factor(testset$svm_pred))
testset$pred_svm_prob<-predict(object=svm_model,testset[,],type='prob')
summary(testset$svm_pred)

#plot(testset$svm_pred)

cm <- table(testset$Cancer.Presence, testset$svm_pred)
cm
TP=cm[1]
FP=cm[3]
FN=cm[2]
TN=cm[4]
Accuracy=round((TP+TN)/(TP+FP+FN+TN)*100,digits=2)
Precision=round(TP/(TP+FP),digits=2)
Recall=round(TP/(TP+FN),digits=2)
Sensitivity=Recall
TPR=Recall
Fmeasure=round((2*Precision*Recall)/(Precision+Recall),digits=2)
Specificity=round(TN/(TN+FP),digits=2)
FPR=round(FP/(FP+TN),digits=2)
ROC=round((1+TPR-FPR)/2,digits=2)
Err_Rate=round((FP+FN)/(TP+FP+FN+TN),digits=2)


svm_5_DF<-data.frame()
svm_5_DF<-data.frame(Sno=Sno, TP=TP, FP=FP, FN=FN, TN=TN, Accuracy=Accuracy, Precision=Precision, Sensitivity_Recall_TPR=TPR,
                     Fmeasure=Fmeasure, Specificity=Specificity, FPR=FPR,ROC_AUC=ROC,Misclassification_Error_Rate=Err_Rate)
Vec_SVM_5_DF<-rbind(Vec_SVM_5_DF,svm_5_DF)
Sno<-Sno+1
}


vec1<-c()
vec2<-c()
vec3<-c()
vec4<-c()
vec5<-c()
vec6<-c()
vec7<-c()
vec8<-c()
vec9<-c()
vec10<-c()
vec11<-c()
vec12<-c()



vec1<-Vec_SVM_5_DF$TP
vec2<-Vec_SVM_5_DF$FP
vec3<-Vec_SVM_5_DF$FN
vec4<-Vec_SVM_5_DF$TN
vec5<- Vec_SVM_5_DF$Accuracy
vec6<-Vec_SVM_5_DF$Precision
vec7<-Vec_SVM_5_DF$Sensitivity_Recall_TPR
vec8<- Vec_SVM_5_DF$Fmeasure
vec9<-Vec_SVM_5_DF$Specificity
vec10<-Vec_SVM_5_DF$FPR
vec11<-Vec_SVM_5_DF$ROC_AUC
vec12<-Vec_SVM_5_DF$Misclassification_Error_Rate


vec1<-unlist(vec1)
vec2<-unlist(vec2)
vec3<-unlist(vec3)
vec4<-unlist(vec4)
vec5<-unlist(vec5)
vec6<-unlist(vec6)
vec7<-unlist(vec7)
vec8<-unlist(vec8)
vec9<-unlist(vec9)
vec10<-unlist(vec10)
vec11<-unlist(vec11)
vec12<-unlist(vec12)


vec1<-as.numeric(vec1)
vec2<-as.numeric(vec2)
vec3<-as.numeric(vec3)
vec4<-as.numeric(vec4)
vec5<-as.numeric(vec5)
vec6<-as.numeric(vec6)
vec7<-as.numeric(vec7)
vec8<-as.numeric(vec8)
vec9<-as.numeric(vec9)
vec10<-as.numeric(vec10)
vec11<-as.numeric(vec11)
vec12<-as.numeric(vec12)



TP=round(mean(vec1),digits=0)
FP=round(mean(vec2),digits=0)
FN=round(mean(vec3),digits=0)
TN=round(mean(vec4),digits=0)
Accuracy=round(mean(vec5),digits=2)
Precision=round(mean(vec6),digits=2)
Recall=round(mean(vec7),digits=2)
Sensitivity=Recall
TPR=Recall
Fmeasure=round(mean(vec8),digits=2)
Specificity=round(mean(vec9),digits=2)
FPR=round(mean(vec10),digits=2)
ROC=round(mean(vec11),digits=2)
Err_Rate=round(mean(vec12),digits=2)



svm_5_DF<-data.frame(Sno="Mean", TP=TP, FP=FP, FN=FN, TN=TN, Accuracy=Accuracy, Precision=Precision, Sensitivity_Recall_TPR=TPR,
                     Fmeasure=Fmeasure, Specificity=Specificity, FPR=FPR,ROC_AUC=ROC,Misclassification_Error_Rate=Err_Rate)
Vec_SVM_5_DF<-rbind(Vec_SVM_5_DF,svm_5_DF)
SD_Accuracy<-round(sd(Vec_SVM_5_DF$Accuracy),digits=1)
svm_5_DF<-data.frame(Sno="", TP="", FP="", FN="", TN="", Accuracy=SD_Accuracy, Precision="", Sensitivity_Recall_TPR="",
                     Fmeasure="", Specificity="", FPR="",ROC_AUC="",Misclassification_Error_Rate="")
Vec_SVM_5_DF<-rbind(Vec_SVM_5_DF,svm_5_DF)
###############################
# Writing to a CSV file
###############################
write.csv(Vec_SVM_5_DF,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_SVM_5_Runs.csv",row.names = F)



####################################################
#II. 10 GENES
###################################################
#cancerdata
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_Top10_genes.csv",header = TRUE)

Vec_SVM_5_DF<-c()
Sno=1


for(i in 1:25){
index<-1:nrow(cancerdata)
length(index)
testindex<-sample(index,trunc(length(index)/3))
#testindex
trainset<-cancerdata[-testindex,]
#trainset
testset<-cancerdata[testindex,]
#testset
start_time <- Sys.time()
fitControl<-trainControl(method="LOOCV",number=10,classProbs=TRUE)
#fitControl
svm_model<-train(Cancer.Presence ~.,data=trainset,method='svmRadial',trControl=fitControl,tuneLength=3)
#svm_model
testset$svm_pred<-predict(svm_model,testset[,])
testset$svm_pred
end_time <- Sys.time()
time=end_time-start_time

# time_estimate<-data.frame()
# time_estimate<-data.frame(Pattern="GDS3057_PO_Pearson_SVM",Runtime_10genes=time)
# write.csv(time_estimate,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_SVM_10g_Runtime.csv",row.names=F)
# 


confusionMatrix(as.factor(testset$Cancer.Presence),as.factor(testset$svm_pred))
testset$pred_svm_prob<-predict(object=svm_model,testset[,],type='prob')
summary(testset$svm_pred)

#plot(testset$svm_pred)
cm <- table(testset$Cancer.Presence, testset$svm_pred)
cm
TP=cm[1]
FP=cm[3]
FN=cm[2]
TN=cm[4]
Accuracy=round((TP+TN)/(TP+FP+FN+TN)*100,digits=2)
Precision=round(TP/(TP+FP),digits=2)
Recall=round(TP/(TP+FN),digits=2)
Sensitivity=Recall
TPR=Recall
Fmeasure=round((2*Precision*Recall)/(Precision+Recall),digits=2)
Specificity=round(TN/(TN+FP),digits=2)
FPR=round(FP/(FP+TN),digits=2)
ROC=round((1+TPR-FPR)/2,digits=2)
Err_Rate=round((FP+FN)/(TP+FP+FN+TN),digits=2)


svm_5_DF<-data.frame()
svm_5_DF<-data.frame(Sno=Sno, TP=TP, FP=FP, FN=FN, TN=TN, Accuracy=Accuracy, Precision=Precision, Sensitivity_Recall_TPR=TPR,
                     Fmeasure=Fmeasure, Specificity=Specificity, FPR=FPR,ROC_AUC=ROC,Misclassification_Error_Rate=Err_Rate)
Vec_SVM_5_DF<-rbind(Vec_SVM_5_DF,svm_5_DF)
Sno<-Sno+1
}


vec1<-c()
vec2<-c()
vec3<-c()
vec4<-c()
vec5<-c()
vec6<-c()
vec7<-c()
vec8<-c()
vec9<-c()
vec10<-c()
vec11<-c()
vec12<-c()



vec1<-Vec_SVM_5_DF$TP
vec2<-Vec_SVM_5_DF$FP
vec3<-Vec_SVM_5_DF$FN
vec4<-Vec_SVM_5_DF$TN
vec5<- Vec_SVM_5_DF$Accuracy
vec6<-Vec_SVM_5_DF$Precision
vec7<-Vec_SVM_5_DF$Sensitivity_Recall_TPR
vec8<- Vec_SVM_5_DF$Fmeasure
vec9<-Vec_SVM_5_DF$Specificity
vec10<-Vec_SVM_5_DF$FPR
vec11<-Vec_SVM_5_DF$ROC_AUC
vec12<-Vec_SVM_5_DF$Misclassification_Error_Rate


vec1<-unlist(vec1)
vec2<-unlist(vec2)
vec3<-unlist(vec3)
vec4<-unlist(vec4)
vec5<-unlist(vec5)
vec6<-unlist(vec6)
vec7<-unlist(vec7)
vec8<-unlist(vec8)
vec9<-unlist(vec9)
vec10<-unlist(vec10)
vec11<-unlist(vec11)
vec12<-unlist(vec12)


vec1<-as.numeric(vec1)
vec2<-as.numeric(vec2)
vec3<-as.numeric(vec3)
vec4<-as.numeric(vec4)
vec5<-as.numeric(vec5)
vec6<-as.numeric(vec6)
vec7<-as.numeric(vec7)
vec8<-as.numeric(vec8)
vec9<-as.numeric(vec9)
vec10<-as.numeric(vec10)
vec11<-as.numeric(vec11)
vec12<-as.numeric(vec12)



TP=round(mean(vec1),digits=0)
FP=round(mean(vec2),digits=0)
FN=round(mean(vec3),digits=0)
TN=round(mean(vec4),digits=0)
Accuracy=round(mean(vec5),digits=2)
Precision=round(mean(vec6),digits=2)
Recall=round(mean(vec7),digits=2)
Sensitivity=Recall
TPR=Recall
Fmeasure=round(mean(vec8),digits=2)
Specificity=round(mean(vec9),digits=2)
FPR=round(mean(vec10),digits=2)
ROC=round(mean(vec11),digits=2)
Err_Rate=round(mean(vec12),digits=2)



svm_5_DF<-data.frame(Sno="Mean", TP=TP, FP=FP, FN=FN, TN=TN, Accuracy=Accuracy, Precision=Precision, Sensitivity_Recall_TPR=TPR,
                     Fmeasure=Fmeasure, Specificity=Specificity, FPR=FPR,ROC_AUC=ROC,Misclassification_Error_Rate=Err_Rate)
Vec_SVM_5_DF<-rbind(Vec_SVM_5_DF,svm_5_DF)
SD_Accuracy<-round(sd(Vec_SVM_5_DF$Accuracy),digits=1)
svm_5_DF<-data.frame(Sno="", TP="", FP="", FN="", TN="", Accuracy=SD_Accuracy, Precision="", Sensitivity_Recall_TPR="",
                     Fmeasure="", Specificity="", FPR="",ROC_AUC="",Misclassification_Error_Rate="")
Vec_SVM_5_DF<-rbind(Vec_SVM_5_DF,svm_5_DF)
###############################
# Writing to a CSV file
###############################
write.csv(Vec_SVM_5_DF,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_SVM_10_Runs.csv",row.names = F)



####################################################
#III. 15 GENES
###################################################
#cancerdata
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_Top15_genes.csv",header = TRUE)

Vec_SVM_5_DF<-c()
Sno=1

for(i in 1:25){
index<-1:nrow(cancerdata)
length(index)
testindex<-sample(index,trunc(length(index)/3))
#testindex
trainset<-cancerdata[-testindex,]
#trainset
testset<-cancerdata[testindex,]
#testset
start_time <- Sys.time()
fitControl<-trainControl(method="LOOCV",number=10,classProbs=TRUE)
#fitControl
svm_model<-train(Cancer.Presence ~.,data=trainset,method='svmRadial',trControl=fitControl,tuneLength=3)
#svm_model
testset$svm_pred<-predict(svm_model,testset[,])
testset$svm_pred
end_time <- Sys.time()
time=end_time-start_time

# time_estimate<-data.frame()
# time_estimate<-data.frame(Pattern="GDS3057_PO_Pearson_SVM",Runtime_15genes=time)
# write.csv(time_estimate,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_SVM_Runtime.csv",row.names=F)

confusionMatrix(as.factor(testset$Cancer.Presence),as.factor(testset$svm_pred))
testset$pred_svm_prob<-predict(object=svm_model,testset[,],type='prob')
summary(testset$svm_pred)

#plot(testset$svm_pred)
cm <- table(testset$Cancer.Presence, testset$svm_pred)
cm
TP=cm[1]
FP=cm[3]
FN=cm[2]
TN=cm[4]
Accuracy=round((TP+TN)/(TP+FP+FN+TN)*100,digits=2)
Precision=round(TP/(TP+FP),digits=2)
Recall=round(TP/(TP+FN),digits=2)
Sensitivity=Recall
TPR=Recall
Fmeasure=round((2*Precision*Recall)/(Precision+Recall),digits=2)
Specificity=round(TN/(TN+FP),digits=2)
FPR=round(FP/(FP+TN),digits=2)
ROC=round((1+TPR-FPR)/2,digits=2)
Err_Rate=round((FP+FN)/(TP+FP+FN+TN),digits=2)


svm_5_DF<-data.frame()
svm_5_DF<-data.frame(Sno=Sno, TP=TP, FP=FP, FN=FN, TN=TN, Accuracy=Accuracy, Precision=Precision, Sensitivity_Recall_TPR=TPR,
                     Fmeasure=Fmeasure, Specificity=Specificity, FPR=FPR,ROC_AUC=ROC,Misclassification_Error_Rate=Err_Rate)
Vec_SVM_5_DF<-rbind(Vec_SVM_5_DF,svm_5_DF)
Sno<-Sno+1
}


vec1<-c()
vec2<-c()
vec3<-c()
vec4<-c()
vec5<-c()
vec6<-c()
vec7<-c()
vec8<-c()
vec9<-c()
vec10<-c()
vec11<-c()
vec12<-c()



vec1<-Vec_SVM_5_DF$TP
vec2<-Vec_SVM_5_DF$FP
vec3<-Vec_SVM_5_DF$FN
vec4<-Vec_SVM_5_DF$TN
vec5<- Vec_SVM_5_DF$Accuracy
vec6<-Vec_SVM_5_DF$Precision
vec7<-Vec_SVM_5_DF$Sensitivity_Recall_TPR
vec8<- Vec_SVM_5_DF$Fmeasure
vec9<-Vec_SVM_5_DF$Specificity
vec10<-Vec_SVM_5_DF$FPR
vec11<-Vec_SVM_5_DF$ROC_AUC
vec12<-Vec_SVM_5_DF$Misclassification_Error_Rate


vec1<-unlist(vec1)
vec2<-unlist(vec2)
vec3<-unlist(vec3)
vec4<-unlist(vec4)
vec5<-unlist(vec5)
vec6<-unlist(vec6)
vec7<-unlist(vec7)
vec8<-unlist(vec8)
vec9<-unlist(vec9)
vec10<-unlist(vec10)
vec11<-unlist(vec11)
vec12<-unlist(vec12)


vec1<-as.numeric(vec1)
vec2<-as.numeric(vec2)
vec3<-as.numeric(vec3)
vec4<-as.numeric(vec4)
vec5<-as.numeric(vec5)
vec6<-as.numeric(vec6)
vec7<-as.numeric(vec7)
vec8<-as.numeric(vec8)
vec9<-as.numeric(vec9)
vec10<-as.numeric(vec10)
vec11<-as.numeric(vec11)
vec12<-as.numeric(vec12)



TP=round(mean(vec1),digits=0)
FP=round(mean(vec2),digits=0)
FN=round(mean(vec3),digits=0)
TN=round(mean(vec4),digits=0)
Accuracy=round(mean(vec5),digits=2)
Precision=round(mean(vec6),digits=2)
Recall=round(mean(vec7),digits=2)
Sensitivity=Recall
TPR=Recall
Fmeasure=round(mean(vec8),digits=2)
Specificity=round(mean(vec9),digits=2)
FPR=round(mean(vec10),digits=2)
ROC=round(mean(vec11),digits=2)
Err_Rate=round(mean(vec12),digits=2)



svm_5_DF<-data.frame(Sno="Mean", TP=TP, FP=FP, FN=FN, TN=TN, Accuracy=Accuracy, Precision=Precision, Sensitivity_Recall_TPR=TPR,
                     Fmeasure=Fmeasure, Specificity=Specificity, FPR=FPR,ROC_AUC=ROC,Misclassification_Error_Rate=Err_Rate)
Vec_SVM_5_DF<-rbind(Vec_SVM_5_DF,svm_5_DF)
SD_Accuracy<-round(sd(Vec_SVM_5_DF$Accuracy),digits=1)
svm_5_DF<-data.frame(Sno="", TP="", FP="", FN="", TN="", Accuracy=SD_Accuracy, Precision="", Sensitivity_Recall_TPR="",
                     Fmeasure="", Specificity="", FPR="",ROC_AUC="",Misclassification_Error_Rate="")
Vec_SVM_5_DF<-rbind(Vec_SVM_5_DF,svm_5_DF)
###############################
# Writing to a CSV file
###############################
write.csv(Vec_SVM_5_DF,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/Significant-Patterns/PO/Pearson/GDS3057-PO_Pearson_SVM_15_Runs.csv",row.names = F)



###############################
# pROGRAM ENDS HERE
###############################

















####################################################
#IV. 20 GENES
#################################################### IV.20 genes
#cancerdata
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3257-LungCancer/InSignificant-Patterns/EO/IG/GDS3257-EO_IG_Top20_genes.csv",header = TRUE)
index<-1:nrow(cancerdata)
length(index)
testindex<-sample(index,trunc(length(index)/3))
#testindex
trainset<-cancerdata[-testindex,]
#trainset
testset<-cancerdata[testindex,]
#testset
fitControl<-trainControl(method="LOOCV",number=10,classProbs=TRUE)
#fitControl
svm_model<-train(Cancer.Presence ~.,data=trainset,method='svmRadial',trControl=fitControl,tuneLength=3)
#svm_model
testset$svm_pred<-predict(svm_model,testset[,])
testset$svm_pred
confusionMatrix(as.factor(testset$Cancer.Presence),as.factor(testset$svm_pred))
testset$pred_svm_prob<-predict(object=svm_model,testset[,],type='prob')
summary(testset$svm_pred)
plot(testset$svm_pred)


####################################################
#V. 25 GENES
###################################################
#cancerdata
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3257-LungCancer/InSignificant-Patterns/EO/IG/GDS3257-EO_IG_Top25_genes.csv",header = TRUE)
index<-1:nrow(cancerdata)
length(index)
testindex<-sample(index,trunc(length(index)/3))
#testindex
trainset<-cancerdata[-testindex,]
#trainset
testset<-cancerdata[testindex,]
#testset
fitControl<-trainControl(method="LOOCV",number=10,classProbs=TRUE)
#fitControl
svm_model<-train(Cancer.Presence ~.,data=trainset,method='svmRadial',trControl=fitControl,tuneLength=3)
#svm_model
testset$svm_pred<-predict(svm_model,testset[,])
testset$svm_pred
confusionMatrix(as.factor(testset$Cancer.Presence),as.factor(testset$svm_pred))
testset$pred_svm_prob<-predict(object=svm_model,testset[,],type='prob')
summary(testset$svm_pred)
plot(testset$svm_pred)