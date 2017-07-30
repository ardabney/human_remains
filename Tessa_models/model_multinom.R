# Change Based on User's Preference 
setwd("~/Desktop")

#load in needed libraries
library(phyloseq)
library(nnet)
library(car)
library(caret)

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

#Testing & Training Sets
set.seed(101)
train<-sample(1:120,80)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

# Cross Validation: Feature Selection with Filter Method (P-Values)
# We'll use 8-fold CV since 80 divides nicely into 8 folds
# ********************** #

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
f_sizes = c(seq(10, 935, by = 10), 935)
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
plot(f_sizes, acc_f_m_1, type = "both")

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

# Bootstrap CI For Accuracy: 
# Stand. Dev. 0.1152175 (approx. 12%)
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
quantile(acc_b, c(0.025, 0.975))
sd(acc_b)

# Multinomial Model CV w/ Forward Stepwise Method Feature Selection (8-fold)
# Note: This takes a long time to run, recommend changing f_sizes smaller value
# *************** #
f_sizes = c(1:20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 934) #or 1:15, 1:20, 1:50, etc.
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
  for(j in f_sizes){
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
num_f <- 1:which(acc_f_m_2 == max(acc_f_m_2))[1]
points(f_sizes, acc_f_m_2)
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
confusionMatrix(pred_pmi, dta_test$Estimated_PMI)
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
quantile(acc_b, c(0.025, 0.975))
sd(acc_b)
