library(e1071)
library(caTools)
library(caret)




# I.5 genes
#cancerdata
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top5_genes.csv",header = TRUE)
Vec_SVM_5_DF<-c()
Sno=1


####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(10)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(20)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(30)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(40)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(50)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))

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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(60)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(70)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(80)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(90)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))

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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(100)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(110)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(120)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(130)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(140)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(150)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(160)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(170)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(180)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(190)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(200)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(210)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(220)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(230)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(240)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
set.seed(250)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################
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
Vec_SVM_5_DF
###############################
# Writing to a CSV file
###############################
write.csv(Vec_SVM_5_DF,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_NB_5_Runs.csv",row.names = F)






# II.10 genes
#cancerdata
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top10_genes.csv",header = TRUE)
Vec_SVM_5_DF<-c()
Sno=1


####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(10)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(20)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(30)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(40)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(50)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(60)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(70)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(80)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(90)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(100)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(110)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(120)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(130)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(140)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(150)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(160)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(170)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(180)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(190)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(200)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(210)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(220)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(230)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(240)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
set.seed(250)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################
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
Vec_SVM_5_DF
###############################
# Writing to a CSV file
###############################
write.csv(Vec_SVM_5_DF,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_NB_10_Runs.csv",row.names = F)








# III.15 genes
#cancerdata
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top15_genes.csv",header = TRUE)
Vec_SVM_5_DF<-c()
Sno=1


####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(10)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(20)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(30)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(40)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(50)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(60)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(70)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(80)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(90)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(100)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(110)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(120)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(130)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(140)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))

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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(150)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(160)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(170)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(180)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(190)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(200)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(210)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(220)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(230)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(240)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################




####################################################
#   single run
####################################################
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
set.seed(250)  # Setting Seed
classifier_cl <- naiveBayes(Cancer.Presence~ ., data = train_cl)
#classifier_cl
# Predicting on test data'
y_pred <- predict(classifier_cl, newdata = test_cl)
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, y_pred)
cm
# Model Evaluation
confusionMatrix(as.factor(test_cl$Cancer.Presence),as.factor(y_pred))


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
#Vec_SVM_5_DF
Sno<-Sno+1
####################################################
#   End of single run
####################################################
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
Vec_SVM_5_DF
###############################
# Writing to a CSV file
###############################
write.csv(Vec_SVM_5_DF,"E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_NB_15_Runs.csv",row.names = F)

