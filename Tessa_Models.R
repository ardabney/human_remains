# Change Based on User's Preference 
setwd("~/Desktop")

#load in needed libraries
library(phyloseq)
library(nnet)
library(car)
library(caret)
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

# Feature Selection (Note: This section takes a bit of time to run.)
p_val_otu <- NULL
for(i in 1:926){
  selected <- data.frame(meta_dta,otu_dta[,i])
  pmi_fit <- multinom(Estimated_PMI ~ ., selected, model = TRUE)
  test_fit <- Anova(pmi_fit)
  p_val_otu[i] <- test_fit[9,3]
}
select <- which(p_val_otu < 0.05)
select_otu <- otu_dta[,select]

#Creating Test and Training Data
set.seed(101)
train<-sample(1:120,80)
train_mt <-meta_dta[train,]
train_otu  <-select_otu[train,]
test_mt <- meta_dta[-train,]
test_otu <- select_otu[-train,]

# Model creation and Testing
training <- data.frame(train_mt,train_otu)
pmi_mult <- multinom(Estimated_PMI ~ ., data = training, MaxNWts=1500)
testing <- data.frame(test_mt,test_otu)
pred_pmi <- predict(pmi_mult, newdata = testing)
(c_matr <- confusionMatrix(pred_pmi, testing$Estimated_PMI))

#Cross Validation: Feature Selection with Filter Method

#Testing & Training Sets
set.seed(101)
train<-sample(1:120,80)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

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
f_sizes = c(seq(10, 926, by = 10), 926)
acc_f_b = matrix(NA, nrow = length(f_sizes), ncol = 8) # Misclassification Rate
for(b in 1:8) {
  cat(".")
  MT_test = train_mt[fold_samples[[b]],]
  OTU_test = train_otu[fold_samples[[b]],] 
  MT_train = train_mt[-fold_samples[[b]],]
  OTU_train = train_otu[-fold_samples[[b]],] 
  
  ## Feature selection.
p_val_otu <- NULL
for(i in 1:926){
    selected <- data.frame(MT_train,OTU_train[,i])
    pmi_fit <- multinom(Estimated_PMI ~ ., selected, model = TRUE)
    test_fit <- Anova(pmi_fit)
    p_val_otu[i] <- test_fit[9,3]
}
oo = order(p_val_otu)
for(f in 1:length(f_sizes)) {
    feature_set = oo[1:f_sizes[f]]
    otu_train = OTU_train[,feature_set]
    otu_test = OTU_test[,feature_set]
    dta_test = data.frame(MT_test,otu_test)
    dta_train = data.frame(MT_train,otu_train)
    pmi_mult <- multinom(Estimated_PMI ~ ., data = dta_train, MaxNWts=15000)
    pred_pmi <- predict(pmi_mult, newdata = dta_test)
    c_matr <- confusionMatrix(pred_pmi, dta_test$Estimated_PMI)
    acc_f_b[f, b] = 1 - c_matr$overall[1]
    
  }
}
acc_f = rowMeans(acc_f_b)

# 420 features appears to be the optimal number of classifiers - not sure, if done correctly
p_val_otu <- NULL
for(i in 1:926){
  selected <- data.frame(train_mt,train_otu[,i])
  pmi_fit <- multinom(Estimated_PMI ~ ., selected, model = TRUE)
  test_fit <- Anova(pmi_fit)
  p_val_otu[i] <- test_fit[9,3]
}
oo = order(p_val_otu)
select_otu <- train_otu[,oo[1:420]]

training <- data.frame(train_mt,train_otu)
pmi_mult <- multinom(Estimated_PMI ~ ., data = training, MaxNWts=20000)
testing <- data.frame(test_mt,test_otu)
pred_pmi <- predict(pmi_mult, newdata = testing)
table(pred_pmi, testing$Estimated_PMI) 
(c_matr <- confusionMatrix(pred_pmi, testing$Estimated_PMI))

# Multinomial Model w/ Forward Stepwise Method Feature Selection (note: although it has a extremely high - 95% accuracy - this model doesn't use any of the meta variables, so I'm not sure it is the best to use)
# This method stops when the accuracy stops increasing

#Testing & Training Sets
set.seed(101)
train<-sample(1:120,80)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

train_set <- data.frame(train_mt,train_otu)
test_set <- data.frame(test_mt,test_otu)
decr = TRUE
old = 0
select <- 1
decr = TRUE
while(decr == TRUE){
  acc <- NULL
  for(i in (1:935)[-select]){
    pmi_mult <- multinom(Estimated_PMI ~ ., data = train_set[,c(select,i)], MaxNWts=2000)
    pred_pmi <- predict(pmi_mult, newdata = test_set)
    c_matr <- confusionMatrix(pred_pmi, test_set$Estimated_PMI)
    acc[i] <- c_matr$overall[1]
  }
 new_select <- which(acc == max(acc, na.rm = TRUE))[1]
 new <- max(acc, na.rm = TRUE)
 dif <- new - old
 if(dif <= 0){
   decr = FALSE
 }
 if(dif > 0){
 select <- c(select, new_select)
 old <- new
 }
}

pmi_mult <- multinom(Estimated_PMI ~ ., data = train_set[,c(select)], MaxNWts=2000)
pred_pmi <- predict(pmi_mult, newdata = test_set)
confusionMatrix(pred_pmi, test_set$Estimated_PMI)

#Note: Based on resampling, it would appear this model isn't robust with accuracy ranging from .35 to .95 depending on the seed set 
acc <- matrix(NA, 1000,1)
for(i in 1:1000){
  set.seed(i)
  train<-sample(1:120,80)
  train_mt <- meta_dta[train,]
  train_otu <-  otu_dta[train,]
  test_mt <- meta_dta[-train,]
  test_otu <- otu_dta[-train,]
  train_set <- data.frame(train_mt,train_otu)
  test_set <- data.frame(test_mt,test_otu)
  pmi_mult <- multinom(Estimated_PMI ~ ., data = train_set[,c(select)], MaxNWts=2000)
  pred_pmi <- predict(pmi_mult, newdata = test_set)
  c_matr <- confusionMatrix(pred_pmi, test_set$Estimated_PMI)
  acc[i] <- c_matr$overall[1]
}
hist(acc)
table(acc)
boxplot(acc)

# Multinomial Model w/ Forward Stepwise Method Feature Selection (note: although it has a extremely high - 95% accuracy - this model doesn't use any of the meta variables, so I'm not sure it is the best to use)
# This method calculates all accuracies and picks highest

#Testing & Training Sets
set.seed(101)
train<-sample(1:120,80)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

train_set <- data.frame(train_mt,train_otu)
test_set <- data.frame(test_mt,test_otu)
decr = TRUE
old = 0
select <- 1
decr = TRUE
acc_j = NULL
for(j in 1:934){
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
  acc_j[j] <- new
  select <- c(select, new_select)
  old <- new
}

#Note: Manually Stopped after 5 hours, only 156 variables were processed, but this was plenty to see a general pattern

# 97.5% accurate 
pmi_mult <- multinom(Estimated_PMI ~ ., data = train_set[,c(select[1:14])], MaxNWts=2000)
pred_pmi <- predict(pmi_mult, newdata = test_set)
confusionMatrix(pred_pmi, test_set$Estimated_PMI)

#Note: Based on resampling, it would appear this model isn't robust with accuracy ranging from .375 to .975 depending on the seed set 
acc <- matrix(NA, 1000,1)
for(i in 1:1000){
  set.seed(i)
  train<-sample(1:120,80)
  train_mt <- meta_dta[train,]
  train_otu <-  otu_dta[train,]
  test_mt <- meta_dta[-train,]
  test_otu <- otu_dta[-train,]
  train_set <- data.frame(train_mt,train_otu)
  test_set <- data.frame(test_mt,test_otu)
  pmi_mult <- multinom(Estimated_PMI ~ ., data = train_set[,c(select[1:14])], MaxNWts=2000)
  pred_pmi <- predict(pmi_mult, newdata = test_set)
  c_matr <- confusionMatrix(pred_pmi, test_set$Estimated_PMI)
  acc[i] <- c_matr$overall[1]
}
hist(acc)
table(acc)
boxplot(acc)


# Backward Stepwise Feature Selection w/ Multinomial model
# Note: Haven't run yet, may need to debug
set.seed(101)
train<-sample(1:120,80)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

train_set <- data.frame(train_mt,train_otu)
test_set <- data.frame(test_mt,test_otu)
decr = TRUE
old = 0
pmi <- 1
decr = TRUE
remv <- 0
while(decr == TRUE){
  acc <- NULL
  for(i in (1:935)[-c(pmi,remv)]){
    pmi_mult <- multinom(Estimated_PMI ~ ., data = train_set[,-c(remv,i)], MaxNWts=20000)
    pred_pmi <- predict(pmi_mult, newdata = test_set)
    c_matr <- confusionMatrix(pred_pmi, test_set$Estimated_PMI)
    acc[i] <- c_matr$overall[1]
  }
  new_remove <- which(acc == max(acc, na.rm = TRUE))[1]
  new <- max(acc, na.rm = TRUE)
  dif <- new - old
  if(dif <= 0){
    decr = FALSE
  }
  if(dif > 0){
    remv <- c(remv, new_remove)
    old <- new
  }
}

pmi_mult <- multinom(Estimated_PMI ~ ., data = train_set[,c(select)], MaxNWts=2000)
pred_pmi <- predict(pmi_mult, newdata = test_set)
confusionMatrix(pred_pmi, test_set$Estimated_PMI)

# KNN - Error thrown regarding NAs, but NAs were removed so I'm unsure of how to proceed 
library(class)
otu_dta<- otu_dta[!is.na(meta_dta$BMI),]
meta_dta<-meta_dta[!is.na(meta_dta$BMI),]
pred_KNN <- rep(NA, 935)
set.seed(101)
train<-sample(1:114,75)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]
train_set <- data.frame(train_mt,train_otu)
test_set <- data.frame(test_mt,test_otu)
for(k in 1:935){
  pred_KNN_k <- knn(train_set, test_set, train_set$Estimated_PMI, k=k)
  pred_KNN[k] <- mean(pred_KNN_k == test_set$Estimated_PMI)
}
plot(1:935, pred_KNN)


# Forward Stepwise Method with random forest - Infinite Loop possibly, appears to not work
otu_dta<- otu_dta[!is.na(meta_dta$BMI),]
meta_dta<-meta_dta[!is.na(meta_dta$BMI),]
set.seed(101)
train<-sample(1:114,75)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

train_set <- data.frame(train_mt,train_otu)
test_set <- data.frame(test_mt,test_otu)
decr = TRUE
old = 0
select <- 1
decr = TRUE
while(decr == TRUE){
  acc <- NULL
  for(i in (1:935)[-select]){
    pmi_rf <- randomForest(Estimated_PMI ~ ., data = train_set[,c(select,i)], MaxNWts=2000)
    pred_rf <- predict(pmi_rf, newdata = test_set)
    c_matr <- confusionMatrix(pred_rf, test_set$Estimated_PMI)
    acc[i] <- c_matr$overall[1]
  }
  new_select <- which(acc == max(acc, na.rm = TRUE))[1]
  new <- max(acc, na.rm = TRUE)
  dif <- new - old
  if(dif <= 0){
    decr = FALSE
  }
  if(dif > 0){
    select <- c(select, new_select)
    old <- new
  }
}

pmi_mult <- multinom(Estimated_PMI ~ ., data = train_set[,c(select)], MaxNWts=2000)
pred_pmi <- predict(pmi_mult, newdata = test_set)
confusionMatrix(pred_pmi, test_set$Estimated_PMI)

# Naive Bayes model - Forward Stepwise Feature Selection
library(e1071)
set.seed(101)
train<-sample(1:120,80)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

train_set <- data.frame(train_mt,train_otu)
test_set <- data.frame(test_mt,test_otu)
decr = TRUE
old = 0
select <- 1
decr = TRUE
while(decr == TRUE){
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
  if(dif <= 0){
    decr = FALSE
  }
  if(dif > 0){
    select <- c(select, new_select)
    old <- new
  }
}
pmi_bayes <- naiveBayes(Estimated_PMI ~ ., data = train_set[,c(select)])
pred_pmi <- predict(pmi_bayes, newdata = test_set)
confusionMatrix(pred_pmi, test_set$Estimated_PMI)

#Robustness of Naive Bayes Model: Also not robust, model has lower accuracy overall from the multinomial model
acc <- matrix(NA, 1000,1)
for(i in 1:1000){
  set.seed(i)
  train<-sample(1:120,80)
  train_mt <- meta_dta[train,]
  train_otu <-  otu_dta[train,]
  test_mt <- meta_dta[-train,]
  test_otu <- otu_dta[-train,]
  train_set <- data.frame(train_mt,train_otu)
  test_set <- data.frame(test_mt,test_otu)
  pmi_bayes <- naiveBayes(Estimated_PMI ~ ., data = train_set[,c(select)])
  pred_pmi <- predict(pmi_bayes, newdata = test_set)
  c_matr <- confusionMatrix(pred_pmi, test_set$Estimated_PMI)
  acc[i] <- c_matr$overall[1]
}
hist(acc)
table(acc)
boxplot(acc)

# Random Forest w/Boruta Feature Selection (note: I'm not sure how to do CV with this)
install.packages("Boruta")
library(Boruta)
mta_dta<-meta_dta[!is.na(meta_dta$BMI),]
otu_dta<- otu_dta[!is.na(meta_dta$BMI),]

data <- data.frame(mta_dta,otu_dta)
train<-sample(1:114,75)
dta_train <- data[train,]
dta_test <- data[-train,]
boruta.train <- Boruta(Estimated_PMI ~ ., data = dta_train, doTrace = 2)
print(boruta.train)
final.boruta <- TentativeRoughFix(boruta.train)
f <- getConfirmedFormula(final.boruta)
randomForest(f, data=dta_train)

# Multinomial Model CV w/ Forward Stepwise Method Feature Selection

#Testing & Training Sets
set.seed(101)
train<-sample(1:120,80)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

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
f_sizes = 1:934
acc_f_b <- matrix(NA, length(f_sizes), 8) # Accuracy
selected <- NULL
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
  for(j in 1:934){
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
    acc_f_b[j] <- new
    select <- c(select, new_select)
    old <- new
  }
}
acc_f = rowMeans(acc_f_b)

#Code for Random Forest CV w/ Feature Selection based on Importance Function
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

# We'll use 5-fold CV since 75 divides nicely into 5 folds

# Divide Samples into 8 folds
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
  oo <- order(importance(pmi_rf))
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
acc_f = rowMeans(acc_f_b)

