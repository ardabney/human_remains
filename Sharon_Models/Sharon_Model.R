# Change Based on User's Preference 
setwd("C:\\Users\\blueb\\Desktop\\STAT 485")

#load in needed libraries
library(phyloseq)
library(nnet)
library(car)
library(caret)
library(leaps)

#import data
std_dta  <- import_biom("f120_r1k.biom") #whole data
m_dta <- data.frame(sample_data(std_dta)) #meta data
OTU_dta <- data.frame(otu_table(std_dta)) #OTU data

# Adjust variables into correct formats 
m_dta$BMI <- as.numeric(m_dta$BMI)
m_dta$Age <- as.numeric(m_dta$Age)
m_dta$Estimated_PMI <- factor(m_dta$Estimated_PMI,levels=c("12","<24",">24",">48",">72"), ordered=T)
m_dta$Race <- factor(m_dta$Race)
m_dta$MOD <- factor(m_dta$Manner.of.Death,levels=c("Natural","Accident","Suicide","Homicide"),ordered=T)
m_dta$Season <- factor(m_dta$Season,levels=c("Spring","Summer","Autumn","Winter"),ordered=T)
m_dta$Sex <- factor(m_dta$Sex)
m_dta$Weight_Status <- factor(m_dta$Weight_Status,levels=c("Underweight","Normal Weight","Overweight","Obese","Severe Obesity","Morbid Obesity","Super Obese"))
m_dta$Event_Location <- factor(m_dta$Event_Location)

# Adjust Age & BMI into groups/factors for easier analysis 
m_dta$Age_Group <- cut(m_dta$Age,c(0,25,50,75,100))
m_dta$BMI_Group <- cut(m_dta$BMI,c(0,20,30,40,50,60))

# Data Processing - Meta
# Create data Set - one with all the patients' samples combined so we can explore Race, Sex, etc
mt_dta<-m_dta[!duplicated(m_dta$Pack_ID),] #patients
row.names(mt_dta)<-mt_dta$Pack_ID
mt_dta$Estimated_PMI[mt_dta$Estimated_PMI == "12"] <- "<24"
mt_dta$Estimated_PMI <- factor(as.character(mt_dta$Estimated_PMI))
meta_dta <-subset(mt_dta,select=c("Estimated_PMI","Race","MOD","Season","Sex","Weight_Status","Event_Location","BMI_Group","Age_Group"))

#Data Processing - OTU
ave_otu_dta<-matrix(NA,nrow=dim(OTU_dta)[1],ncol=length(unique(m_dta$Pack_ID)))
row.names(ave_otu_dta)<-rownames(OTU_dta)
colnames(ave_otu_dta)<-unique(m_dta$Pack_ID)
for (i in unique(m_dta$Pack_ID)){
  temp_id<-rownames(m_dta[m_dta$Pack_ID%in%i,])
  ave_otu_dta[,i]<-rowMeans(OTU_dta[,temp_id])
}
otu_dta <-t(data.frame(ave_otu_dta))


#seperate training and testing
set.seed(120)
train<-sample(1:120,80)
train_mt <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_mt <- meta_dta[-train,]
test_otu <- otu_dta[-train,]


#data frame selected data
select1 <- cbind(train_mt, train_otu)
select2 <- na.omit(select1) #76 obs

#select predictors
null <- multinom(MOD~ 1,select2,model = TRUE)
full <- multinom(MOD~.,select2, model = TRUE, MaxNWts=4000)
step.fwd <- step(null, scope=list(lower=null, upper=full), direction="forward")  

#backward selection is very slow
#step.bwd <- step(full, scope=list(lower=null, upper=full), direction="backward")

#see the selected predictors
step.fwd

# Multinomial Model creation and Testing
#porblem: code not robusted, I manually inputed the x variables seen in last step, could be improved 
mod_mult <- multinom(MOD ~ Race + denovo43942 + denovo211190 + denovo160622 + 
                       denovo113977 + denovo79567 + denovo159516 + denovo26819 + 
                       denovo13288 + denovo145648, data = select2, MaxNWts=4000)

testing <- data.frame(test_mt,test_otu)
pred_mod <- predict(mod_mult, newdata = testing)
confusionMatrix(pred_mod, testing$MOD)

# ----------------------------------------------------------------------------
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
    mod_fit <- multinom(MOD ~ ., selected, model = TRUE)
    test_fit <- Anova(mod_fit)
    p_val_otu[i] <- test_fit[9,3]
  }
  oo = order(p_val_otu)
  for(f in 1:length(f_sizes)) {
    feature_set = oo[1:f_sizes[f]]
    otu_train = OTU_train[,feature_set]
    otu_test = OTU_test[,feature_set]
    dta_test = data.frame(MT_test,otu_test)
    dta_train = data.frame(MT_train,otu_train)
    mod_mult <- multinom(MOD ~ ., data = dta_train, MaxNWts=15000)
    pred_mod <- predict(mod_mult, newdata = dta_test)
    c_matr <- confusionMatrix(pred_mod, dta_test$MOD)
    acc_f_b[f, b] = 1 - c_matr$overall[1]
    
  }
}
acc_f = rowMeans(acc_f_b)
max(acc_f)

#------------------------------------------------
# 31 features appears to be the optimal number of classifiers
p_val_otu <- NULL
for(i in 1:926){
  selected <- data.frame(train_mt,train_otu[,i])
  mod_fit <- multinom(MOD ~ ., selected, model = TRUE)
  test_fit <- Anova(mod_fit)
  p_val_otu[i] <- test_fit[9,3]
}
oo = order(p_val_otu)
select_otu <- train_otu[,oo[1:31]]

training <- data.frame(train_mt,train_otu)
mod_mult <- multinom(MOD ~ ., data = training, MaxNWts=20000)
testing <- data.frame(test_mt,test_otu)
pred_mod <- predict(mod_mult, newdata = testing)
table(pred_mod, testing$MOD) 
(c_matr <- confusionMatrix(pred_mod, testing$MOD))
# about 30% accuracy 
