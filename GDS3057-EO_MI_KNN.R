library(e1071)
library(caTools)
library(class)
library(kernlab)

# cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top5_genes.csv",header = TRUE)
# cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top10_genes.csv",header = TRUE)
# cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top15_genes.csv",header = TRUE)
# cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top20_genes.csv",header = TRUE)
# cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top25_genes.csv",header = TRUE)




####################################################
#I. 5 GENES
###################################################
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top5_genes.csv",header = TRUE)
# Splitting data into train
# and test data
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:5])
test_scale <- scale(test_cl[, 1:5])
# Fitting KNN Model 
# to training dataset
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 1)
classifier_knn
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# Model Evaluation - Choosing K Calculate out of Sample error
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 3
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 3)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 5
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 5)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 7
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 7)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 9
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 15)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 11
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 19)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm


####################################################
#II. 10 GENES
###################################################
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top10_genes.csv",header = TRUE)
# Splitting data into train
# and test data
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:10])
test_scale <- scale(test_cl[, 1:10])
# Fitting KNN Model 
# to training dataset
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 1)
classifier_knn
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# Model Evaluation - Choosing K Calculate out of Sample error
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 3
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 3)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 5
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 5)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 7
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 7)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 9
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 9)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 11
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 11)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm



####################################################
#III. 15 GENES
###################################################
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top15_genes.csv",header = TRUE)
# Splitting data into train
# and test data
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:15])
test_scale <- scale(test_cl[, 1:15])
# Fitting KNN Model 
# to training dataset
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 1)
classifier_knn
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# Model Evaluation - Choosing K Calculate out of Sample error
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 3
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 3)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 5
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 5)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 7
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 7)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 9
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 9)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# K = 11
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 11)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm


####################################################
#IV. 20 GENES
###################################################
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top20_genes.csv",header = TRUE)
# Splitting data into train
# and test data
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:20])
test_scale <- scale(test_cl[, 1:20])
# Fitting KNN Model 
# to training dataset
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 1)
classifier_knn
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# Model Evaluation - Choosing K Calculate out of Sample error
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 3
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 3)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 5
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 5)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 7
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 7)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 9
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 9)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 11
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 11)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))





####################################################
#V. 25 GENES
###################################################
cancerdata <- read.csv("E:/Ph.D Work/R/Coding-Towards Perfection/GDS3057-Leukemia/InSignificant-Patterns/EO/MI/GDS3057-EO_MI_Top25_genes.csv",header = TRUE)
# Splitting data into train
# and test data
split <- sample.split(cancerdata, SplitRatio = 0.7)
train_cl <- subset(cancerdata, split == "TRUE")
test_cl <- subset(cancerdata, split == "FALSE")
train_scale <- scale(train_cl[, 1:25])
test_scale <- scale(test_cl[, 1:25])
# Fitting KNN Model 
# to training dataset
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 1)
classifier_knn
# Confusion Matrix
cm <- table(test_cl$Cancer.Presence, classifier_knn)
cm
# Model Evaluation - Choosing K Calculate out of Sample error
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 3
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 3)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 5
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 5)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 7
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 7)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 9
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 9)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
# K = 11
classifier_knn <- knn(train = train_scale,
                      test = test_scale,
                      cl = train_cl$Cancer.Presence,
                      k = 11)
misClassError <- mean(classifier_knn != test_cl$Cancer.Presence)
print(paste('Accuracy =', 1-misClassError))
