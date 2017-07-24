# Change Based on User's Preference 
setwd("~/Desktop")

#load in needed libraries
library(phyloseq)
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
# **************************** #

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
acc_f = rowMeans(acc_f_b)
f <- f_sizes[which(acc_f == max(acc_f))]

# 656 Features w/ 41.03% Accuracy
dta_train = data.frame(train_mt, train_otu)
dta_test = data.frame(test_mt,test_otu)
pmi_rf <- randomForest(Estimated_PMI ~ ., data = dta_train, MaxNWts=2000)
oo <- order(importance(pmi_rf), decreasing=TRUE)
feature_set = c(1,oo[1:f])
data_train <- dta_train[,feature_set]
data_test <- dta_test[,feature_set]
pmi_rf <- randomForest(Estimated_PMI ~ ., data = data_train, MaxNWts=2000)
pred_pmi <- predict(pmi_rf, newdata = data_test)
confusionMatrix(pred_pmi, dta_test$Estimated_PMI)

# 95% CI for accuracy is (0.5, 0.9)
B <- 1000
dta <- data.frame(meta_dta, otu_dta)
acc_b <- NULL
for(b in 1:B){
  valid <- FALSE
  train <- sample(1:120,80, replace = TRUE)
  train <- dta[train,]
  test <- sample(1:120,20,replace = TRUE)
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
hist(acc_b)
quantile(acc_b, c(0.025, 0.975))

# Random Forest CV w/ Wrapper Method Feature Selection (Forward Stepwise)
f_sizes = 1:20
acc_j_b <- matrix(NA, length(f_sizes), 5) # Accuracy
selected <- NULL
for(b in 1:5) {
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
      pmi_rf <- randomForest(Estimated_PMI ~ ., data = train_set[,c(select,i)])
      pred_pmi <- predict(pmi_rf, newdata = test_set)
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
acc_j_b
