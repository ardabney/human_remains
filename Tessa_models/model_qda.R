# Change Based on User's Preference 
setwd("~/Desktop")

#load in needed libraries (Note: Change/Add depending on your model type)
library(phyloseq)
library(MASS)
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

#LDA CV w/ Forward Stepwise Wrapper Method
f_sizes = 1:20 #or 1:15, 1:20, 1:50, etc.
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
  for(j in f_sizes){
    acc <- NULL
    for(i in (1:935)[-select]){
      pmi_lda <- qda(Estimated_PMI ~ ., data = train_set[,c(select,i)], prior = c(1,1,1,1)/4)
      pred_pmi <- predict(pmi_lda, newdata = test_set)
      c_matr <- confusionMatrix(pred_pmi$class, test_set$Estimated_PMI)
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
acc_f = rowMeans(acc_f_b)

#Throws Error: Error in qda.default(x, grouping, ...) : rank deficiency in group >48
# Unsure of how to handle