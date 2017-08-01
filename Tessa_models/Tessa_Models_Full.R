# Change Based on User's Preference 
setwd("~/Desktop")

#load in needed libraries
library(phyloseq)
library(nnet)
library(car)
library(caret)
library(e1071)
library(randomForest)

#import data
std_dta  <- import_biom("f120_r1k.biom") #whole data
m_dta <- data.frame(sample_data(std_dta)) #meta data
OTU_dta <- data.frame(otu_table(std_dta)) #OTU data

# Data Processing - Meta & OTU
m_dta$BMI <- as.numeric(m_dta$BMI)
m_dta$Age <- as.numeric(m_dta$Age)
m_dta$Estimated_PMI <- factor(m_dta$Estimated_PMI,levels=c("12","<24",">24",">48",">72"), ordered=T)
m_dta$Race <- factor(m_dta$Race)
m_dta$Manner.of.Death <- factor(m_dta$Manner.of.Death,levels=c("Natural","Accident","Suicide","Homicide"),ordered=T)
m_dta$Season <- factor(m_dta$Season,levels=c("Spring","Summer","Autumn","Winter"),ordered=T)
m_dta$Sex <- factor(m_dta$Sex)
m_dta$Weight_Status <- factor(m_dta$Weight_Status,levels=c("Underweight","Normal Weight","Overweight","Obese","Severe Obesity","Morbid Obesity","Super Obese"))
m_dta$Event_Location <- factor(m_dta$Event_Location)
m_dta$Age_Group <- cut(m_dta$Age,c(0,25,50,75,100))
m_dta$BMI_Group <- cut(m_dta$BMI,c(0,20,30,40,50,60))
mt_dta<-m_dta[!duplicated(m_dta$Pack_ID),] #patients
row.names(mt_dta)<-mt_dta$Pack_ID
mt_dta$Estimated_PMI[mt_dta$Estimated_PMI == "12"] <- "<24"
mt_dta$Estimated_PMI <- factor(as.character(mt_dta$Estimated_PMI))
meta_dta <-subset(mt_dta,select=c("Estimated_PMI","Race","Manner.of.Death","Season","Sex","Weight_Status","Event_Location","BMI_Group","Age_Group"))


ave_otu_dta<-matrix(NA,nrow=dim(OTU_dta)[1],ncol=length(unique(m_dta$Pack_ID)))
row.names(ave_otu_dta)<-rownames(OTU_dta)
colnames(ave_otu_dta)<-unique(m_dta$Pack_ID)
for (i in unique(m_dta$Pack_ID)){
  temp_id<-rownames(m_dta[m_dta$Pack_ID%in%i,])
  ave_otu_dta[,i]<-rowMeans(OTU_dta[,temp_id])
}
otu_dta <-t(round(data.frame(ave_otu_dta)))

# Place to store analysis numbers
feat_num <- NULL
stand_dev <- NULL
CI <- matrix(NA, 4,2)
accuracy <- NULL

#Testing & Training Sets
set.seed(101)
train<-sample(1:120,80)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

#Cross Validation: Feature Selection with Filter Method (P-Values)
# We'll use 8-fold CV since 80 divides nicely into 8 folds

# Divide Samples into 8 folds
new_rand_order <-sample(1:80,80)
fold_1 <- new_rand_order[1:10]
fold_2 <- new_rand_order[11:20]
fold_3 <- new_rand_order[21:30]
fold_4 <- new_rand_order[31:40]
fold_5 <- new_rand_order[41:50]
fold_6 <- new_rand_order[51:60]
fold_7 <- new_rand_order[61:70]
fold_8 <- new_rand_order[71:80]
fold_samples = list(fold_1, fold_2, fold_3, fold_4, fold_5, fold_6, fold_7, fold_8)

# Cross Validation
f_sizes = c(seq(10, 934, by = 10), 934)
acc_f_b = matrix(NA, nrow = length(f_sizes), ncol = 8) # Accuracy 
for(b in 1:8) {
  cat(".")
  MT_test = train_mt[fold_samples[[b]],]
  OTU_test = train_otu[fold_samples[[b]],] 
  MT_train = train_mt[-fold_samples[[b]],]
  OTU_train = train_otu[-fold_samples[[b]],] 
  
  ## Feature selection.
  p_val <- NULL
  dta <- data.frame(MT_train, OTU_train)
  for(i in 2:935){
    selected <- dta[,c(1,i)]
    pmi_fit <- multinom(Estimated_PMI ~ ., selected, model = TRUE)
    test_fit <- Anova(pmi_fit)
    p_val[i] <- test_fit$`Pr(>Chisq)`
  }
  oo = order(p_val)
  for(f in 1:length(f_sizes)) {
    feature_set = c(1,oo[1:f_sizes[f]])
    dta_train = dta[,feature_set]
    dta_test = data.frame(MT_test,OTU_test)
    dta_test = dta_test[,feature_set]
    pmi_mult <- multinom(Estimated_PMI ~ ., data = dta_train, MaxNWts=15000)
    pred_pmi <- predict(pmi_mult, newdata = dta_test)
    c_matr <- confusionMatrix(pred_pmi, dta_test$Estimated_PMI)
    acc_f_b[f, b] = c_matr$overall[1]
  }
}
acc_f_m_1 = rowMeans(acc_f_b)
f <- f_sizes[which(acc_f_m_1 == max(acc_f_m_1))]
feat_num[1] <- f
plot(f_sizes, acc_f_m_1, type = "l")

# 470 features w/ 28.95% accuracy
p_val_otu <- NULL
training <- data.frame(train_mt,train_otu)
for(i in 2:935){
  selected <- training[,c(1,i)]
  pmi_fit <- multinom(Estimated_PMI ~ ., selected, model = TRUE)
  test_fit <- Anova(pmi_fit)
  p_val[i] <- test_fit[9,3]
}
oo = order(p_val)
select_training <- training[,oo[1:f]]
testing <- data.frame(test_mt,test_otu)

pmi_mult <- multinom(Estimated_PMI ~ ., data = select_training, MaxNWts=20000)
pred_pmi <- predict(pmi_mult, newdata = testing)
(c_matr <- confusionMatrix(pred_pmi, testing$Estimated_PMI))
accuracy[1] <- c_matr$overall[1]

# Bootstrap CI for Accuracy: (0.53, 0.84)
# Stand. Dev. 0.08128782
B <- 1000
dta <- data.frame(meta_dta, otu_dta)
acc_b <- NULL
for(b in 1:B){
  train <- sample(1:120,80, replace = TRUE)
  train <- dta[train,]
  test <- sample(1:120,40,replace = TRUE)
  test <- dta[test,]
  pmi_mult <- multinom(Estimated_PMI ~ ., data = train[,oo[1:f]], MaxNWts=20000)
  pred_pmi <- predict(pmi_mult, newdata = test)
  c_matr <- confusionMatrix(pred_pmi, test$Estimated_PMI)
  acc_b[b] <- c_matr$overall[1]
}
hist(acc_b, xlab = "Accuracy", main = "Multinomial w/ Filter Method Accuracies")
CI[1, ] <- quantile(acc_b, c(0.025, 0.975))
stand_dev[1] <- sd(acc_b)

# Multinomial Model CV w/ Forward Stepwise Method Feature Selection (8-fold)
# Be sure to have run lines 54-64 before running the loop
# Note: This takes a long time to run, recommend changing f_sizes smaller value
f_sizes = c(1:20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 934) #or 1:15, 1:20, 1:50, etc.
acc_f_b <- matrix(NA, length(f_sizes), 8) # Accuracy
for(b in 1:8) {
  MT_test = train_mt[fold_samples[[b]],]
  OTU_test = train_otu[fold_samples[[b]],] 
  MT_train = train_mt[-fold_samples[[b]],]
  OTU_train = train_otu[-fold_samples[[b]],] 
  
  ## Feature selection.
  train_set <- data.frame(MT_train,OTU_train)
  test_set <- data.frame(MT_test,OTU_test)
  decr = TRUE
  old = 0
  select <- 1
  decr = TRUE
  acc_j = NULL
  for(j in 1:length(f_sizes)){
    acc <- NULL
    for(i in (1:935)[-select]){
      pmi_mult <- multinom(Estimated_PMI ~ ., data = train_set[,c(select,i)], MaxNWts=20000)
      pred_pmi <- predict(pmi_mult, newdata = test_set)
      c_matr <- confusionMatrix(pred_pmi, test_set$Estimated_PMI)
      acc[i] <- c_matr$overall[1]
    }
    new_select <- which(acc == max(acc, na.rm = TRUE))[1]
    new <- max(acc, na.rm = TRUE)
    dif <- new - old
    acc_f_b[j, b] <- new
    select <- c(select, new_select)
    old <- new
  }
}
acc_f_m_2 = rowMeans(acc_f_b)
num_f <- 1:f_sizes(which(acc_f_m_2 == max(acc_f_m_2)))[1]
feat_num[2] <- f_sizes(which(acc_f_m_2 == max(acc_f_m_2)))[1]
lines(f_sizes, acc_f_m_2)

# Building final model w/ selected features
dta_train = data.frame(train_mt, train_otu)
dta_test = data.frame(test_mt,test_otu)
old = 0
select <- 1
decr = TRUE
for(j in num_f){
  acc <- NULL
  for(i in (1:935)[-select]){
    pmi_mult <- multinom(Estimated_PMI ~ ., data = dta_train[,c(select,i)], MaxNWts=20000)
    pred_pmi <- predict(pmi_mult, newdata = dta_test)
    c_matr <- confusionMatrix(pred_pmi, dta_test$Estimated_PMI)
    acc[i] <- c_matr$overall[1]
  }
  new_select <- which(acc == max(acc, na.rm = TRUE))[1]
  new <- max(acc, na.rm = TRUE)
  dif <- new - old
  select <- c(select, new_select)
  old <- new
}
pmi_mult <- multinom(Estimated_PMI ~ ., data = dta_train[,c(select)], MaxNWts=20000)
pred_pmi <- predict(pmi_mult, newdata = dta_test)
c_matr <- confusionMatrix(pred_pmi, dta_test$Estimated_PMI)
accuracy[2] <- c_matr$overall[1]
# 17 features - 97.5% accuracy

# 95% CI for accuracy is (0.40, 0.85)
# Stand. Dev. 0.1133545
B <- 1000
dta <- data.frame(meta_dta, otu_dta)
acc_b <- NULL
for(b in 1:B){
  train <- sample(1:120,80, replace = TRUE)
  train <- dta[train,]
  test <- sample(1:120,40,replace = TRUE)
  test <- dta[test,]
  pmi_mult <- multinom(Estimated_PMI ~ ., data = train[,c(select)], MaxNWts=20000)
  pred_pmi <- predict(pmi_mult, newdata = test)
  c_matr <- confusionMatrix(pred_pmi, test$Estimated_PMI)
  acc_b[b] <- c_matr$overall[1]
}
hist(acc_b, xlab = "Accuracy", main = "Multinomial w/ Wrapper Method Accuracies")
CI[2, ] <- quantile(acc_b, c(0.025, 0.975))
stand_dev[2] <- sd(acc_b)


# Naive Bayes Model 
# Testing & Training Sets
set.seed(101)
train<-sample(1:120,80)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

# Naives Bayes CV w/ Wrapper Method Feature Selection (Forward Stepwise) - 75% Accurate
# We'll use 8-fold CV since 80 divides nicely into 8 folds

# Divide Samples into 8 folds
new_rand_order <-sample(1:80,80)
fold_1 <- new_rand_order[1:10]
fold_2 <- new_rand_order[11:20]
fold_3 <- new_rand_order[21:30]
fold_4 <- new_rand_order[31:40]
fold_5 <- new_rand_order[41:50]
fold_6 <- new_rand_order[51:60]
fold_7 <- new_rand_order[61:70]
fold_8 <- new_rand_order[71:80]
fold_samples = list(fold_1, fold_2, fold_3, fold_4, fold_5, fold_6, fold_7, fold_8)

# Cross Validation
f_sizes = c(1:20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 934)
acc_j_b <- matrix(NA, length(f_sizes), 8) # Accuracy
selected <- NULL
for(b in 1:8) {
  MT_test = train_mt[fold_samples[[b]],]
  OTU_test = train_otu[fold_samples[[b]],] 
  MT_train = train_mt[-fold_samples[[b]],]
  OTU_train = train_otu[-fold_samples[[b]],] 
  
  ## Feature selection.
  train_set <- data.frame(MT_train,OTU_train)
  test_set <- data.frame(MT_test,OTU_test)
  old = 0
  select <- 1
  decr = TRUE
  for(j in 1:20){
    acc <- NULL
    for(i in (1:935)[-select]){
      pmi_bayes <- naiveBayes(Estimated_PMI ~ ., data = train_set[,c(select,i)])
      pred_pmi <- predict(pmi_bayes, newdata = test_set)
      c_matr <- confusionMatrix(pred_pmi, test_set$Estimated_PMI)
      acc[i] <- c_matr$overall[1]
    }
    new_select <- which(acc == max(acc, na.rm = TRUE))[1]
    new <- max(acc, na.rm = TRUE)
    dif <- new - old
    acc_j_b[j, b] <- new
    select <- c(select, new_select)
    old <- new
  }
}
acc_f_b = rowMeans(acc_j_b)
num_f <- 1:which(acc_f_b == max(acc_f_b))[1]
feat_num[3] <- which(acc_f_b == max(acc_f_b))[1]
lines(f_sizes, acc_f_b)

dta_train = data.frame(train_mt, train_otu)
dta_test = data.frame(test_mt,test_otu)
old = 0
select <- 1
decr = TRUE
for(j in num_f){
  acc <- NULL
  for(i in (1:935)[-select]){
    pmi_bayes <- naiveBayes(Estimated_PMI ~ ., data = dta_train[,c(select,i)])
    pred_pmi <- predict(pmi_bayes, newdata = dta_test)
    c_matr <- confusionMatrix(pred_pmi, dta_test$Estimated_PMI)
    acc[i] <- c_matr$overall[1]
  }
  new_select <- which(acc == max(acc, na.rm = TRUE))[1]
  new <- max(acc, na.rm = TRUE)
  dif <- new - old
  select <- c(select, new_select)
  old <- new
}
pmi_bayes <- naiveBayes(Estimated_PMI ~ ., data = dta_train[,select])
pred_pmi <- predict(pmi_bayes, newdata = dta_test)
c_matr <- confusionMatrix(pred_pmi, dta_test$Estimated_PMI)
accuracy[3] <- c_matr$overall[1]
# 15 features - 75% accuracy

#  95% CI for accuracy is (0.099375, 0.525)
# Stand. Dev. 0.1154541
B <- 1000
dta <- data.frame(meta_dta, otu_dta)
acc_b <- NULL
for(b in 1:B){
  train <- sample(1:120,80, replace = TRUE)
  train <- dta[train,]
  test <- sample(1:120,40,replace = TRUE)
  test <- dta[test,]
  pmi_bayes <- naiveBayes(Estimated_PMI ~ ., data = train[,select])
  pred_pmi <- predict(pmi_bayes, newdata = test)
  c_matr <- confusionMatrix(pred_pmi, test$Estimated_PMI)
  acc_b[b] <- c_matr$overall[1]
}
hist(acc_b, xlab = "Accuracy", main = "Naive Bayes w/ Wrapper Method Accuracies")
CI[3, ] <- quantile(acc_b, c(0.025, 0.975))
stand_dev[3] <- sd(acc_b)

# Random Forests
# Remove NAs
otu_dta<- otu_dta[!is.na(meta_dta$BMI),]
meta_dta<-meta_dta[!is.na(meta_dta$BMI),]

#Testing & Training Sets
set.seed(101)
train <- sample(1:114,75)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

# Random Forest CV w/ Filter Method Feature Selection (Using Importance function as filter)
# We'll use 5-fold CV since 75 divides nicely into 5 folds

# Divide Samples into 5 folds
new_rand_order <-sample(1:75,75)
fold_1 <- new_rand_order[1:15]
fold_2 <- new_rand_order[16:30]
fold_3 <- new_rand_order[31:45]
fold_4 <- new_rand_order[46:60]
fold_5 <- new_rand_order[61:75]
fold_samples = list(fold_1, fold_2, fold_3, fold_4, fold_5)

# Cross Validation
f_sizes = seq(2, 934, by = 2)
acc_f_b = matrix(NA, nrow = length(f_sizes), ncol = 5) # Accuracy
for(b in 1:5) {
  cat(".")
  MT_test = train_mt[fold_samples[[b]],]
  OTU_test = train_otu[fold_samples[[b]],] 
  MT_train = train_mt[-fold_samples[[b]],]
  OTU_train = train_otu[-fold_samples[[b]],] 
  
  dta_train = data.frame(MT_train, OTU_train)
  dta_test = data.frame(MT_test,OTU_test)
  ## Feature selection.
  pmi_rf <- randomForest(Estimated_PMI ~ ., data = dta_train, MaxNWts=2000)
  oo <- order(importance(pmi_rf), decreasing=TRUE)
  for(f in 1:length(f_sizes)) {
    feature_set = c(1,oo[1:f_sizes[f]])
    data_train <- dta_train[,feature_set]
    data_test <- dta_test[,feature_set]
    pmi_rf <- randomForest(Estimated_PMI ~ ., data = data_train, MaxNWts=2000)
    pred_pmi <- predict(pmi_rf, newdata = data_test)
    c_matr <- confusionMatrix(pred_pmi, dta_test$Estimated_PMI)
    acc_f_b[f, b] = c_matr$overall[1]
    
  }
}
acc_f_r = rowMeans(acc_f_b)
f <- f_sizes[which(acc_f_r == max(acc_f_r))]
feat_num[4] <- f
lines(f_sizes, acc_f_r)

dta_train = data.frame(train_mt, train_otu)
dta_test = data.frame(test_mt,test_otu)
pmi_rf <- randomForest(Estimated_PMI ~ ., data = dta_train, MaxNWts=2000)
oo <- order(importance(pmi_rf), decreasing=TRUE)
feature_set = c(1,oo[1:f])
data_train <- dta_train[,feature_set]
data_test <- dta_test[,feature_set]
pmi_rf <- randomForest(Estimated_PMI ~ ., data = data_train, MaxNWts=2000)
pred_pmi <- predict(pmi_rf, newdata = data_test)
c_matr <- confusionMatrix(pred_pmi, dta_test$Estimated_PMI)
accuracy[4] <- c_matr$overall[1]
# 656 Features w/ 41.03% Accuracy

# 95% CI for accuracy is (0.5641026, 0.8684211)
# Stand. Dev. 0.07909884
B <- 1000
dta <- data.frame(meta_dta, otu_dta)
acc_b <- NULL
for(b in 1:B){
  valid <- FALSE
  train <- sample(1:120,80, replace = TRUE)
  train <- dta[train,]
  test <- sample(1:120,40,replace = TRUE)
  test <- dta[test,]
  while(valid == FALSE){
    for(i in table(train$Estimated_PMI) > 0){
      if(i == FALSE){
        train <- sample(1:120,80, replace = TRUE)
        train <- dta[train,]
      }
      else {valid <- TRUE}
    }
  }
  data_train <- train[,feature_set]
  data_test <- test[,feature_set]
  pmi_rf <- randomForest(Estimated_PMI ~ ., data = data_train, MaxNWts=2000, na.action = na.omit)
  pred_pmi <- predict(pmi_rf, newdata = data_test)
  c_matr <- confusionMatrix(pred_pmi, data_test$Estimated_PMI)
  acc_b[b] <- c_matr$overall[1]
}
hist(acc_b, xlab = "Accuracy", main = "Random Forest w/ Filter Method Accuracies")
CI[4, ] <- quantile(acc_b, c(0.025, 0.975))
stand_dev[4] <- sd(acc_b)

#print num of features
feat_num
accuracy
CI
stand_dev
