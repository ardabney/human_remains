# Change Based on User's Preference 
setwd("C:\\Users\\blueb\\Desktop\\STAT 485")

#load in needed libraries
library(phyloseq)
library(nnet)
library(car)
library(caret)
library(leaps)
library(randomForest)
##-------------------------------process initial data-------------------------------------
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


#omit missing values in variable data set  
na.omit(otu_dta) 
na.omit(meta_dta)

#seperate training(80) and testing sets for the whole data set
set.seed(120)
train <- sample(1:120,80)
train_meta <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_meta <- meta_dta[-train,]
test_otu <- otu_dta[-train,]


## ==========Multinomial Filter Method========
##--------------------------------perform Cross Validation--------------------------------------------

#Cross Validation: Feature Selection with Filter Method (P-Values)

# Divide training set into 8 folds
fold_size <-sample(1:80,80)
fold_1 <- fold_size[1:10]
fold_2 <- fold_size[11:20]
fold_3 <- fold_size[21:30]
fold_4 <- fold_size[31:40]
fold_5 <- fold_size[41:50]
fold_6 <- fold_size[51:60]
fold_7 <- fold_size[61:70]
fold_8 <- fold_size[71:80]
fold_samples = list(fold_1, fold_2, fold_3, fold_4, fold_5, fold_6, fold_7, fold_8)


## Now do CV. For each CV fold, we'll build and assess the accuracy of feature set sizes 
f_sizes = c(seq(10, 926, by = 10), 926) #yields 93 numbers
acc_f_b = matrix(NA, nrow = length(f_sizes), ncol = 8) #93*8


for(b in 1:8) { #loop for each fold
  cat(".")
  # meta & otu training and testing set for each fold
  fold_test_meta = train_meta[fold_samples[[b]],]
  fold_test_otu = train_otu[fold_samples[[b]],] 
  fold_train_meta = train_meta[-fold_samples[[b]],]
  fold_train_otu = train_otu[-fold_samples[[b]],] 
  
  #combine meta and otu training and testing for CV
  cross_train = data.frame(fold_train_meta, fold_train_otu)
  cross_test = data.frame(fold_test_meta,fold_test_otu)

  
  #Feature selection
  p_val_otu <- NULL
  for(i in c(1:2,4:ncol(cross_train))){  # 3rd coloum of cross_train is response variable: MOD
    mod_fit <- multinom(cross_train[,3]~cross_train[,i],  model = TRUE)  # 3rd coloum(MOD) ~ other predictive variables in cross-training set
    aov_mod_fit <- Anova(mod_fit) 
    p_val_otu[i] <- aov_mod_fit$`Pr(>Chisq)`  #extract p-value from model's ANOVA
  }
  p_val_otu[3]=2  # assign the 3rd coloum(MOD) a big p-value so that it is not chosen
  oo = order(p_val_otu) #rank p-values of predictors
  
  for(f in 1:length(f_sizes)) {
    feature_set = oo[1:f_sizes[f]]
    
    #data frame of selected variables with top p-values
    select_train = data.frame(fold_train_meta,fold_train_otu)[,c(3,feature_set)] 
    select_test = data.frame(fold_test_meta,fold_test_otu)[,c(3,feature_set)]
    
    fit_mod <- multinom(MOD ~ ., data = select_train, MaxNWts=15000) #fit model using selected predictors 
    pred_mod <- predict(fit_mod, newdata = select_test)
    
    c_matr <- confusionMatrix(pred_mod, select_test$MOD) #one confusion matrix for one f
    
    #my code is as below, Tessa uses 1-c_matr$overall[1]
    acc_f_b[f, b] = c_matr$overall[1] #93*8 matrix
  }
}
acc_f = rowMeans(acc_f_b) #average accuracy for each f of all 8 folds 
f_select= f_sizes[which.max(acc_f)] #find maximum accuracy location
plot(f_sizes, acc_f, type = "both", xlab='feature size (f_sizes)', ylab='accuracy (acc_f)', 
     main='Feature size vs. Accuracy in Multinom (Filter) ')

##------------------------extract selected variables and fit model----------------------------------------

# f = 9 --> f_select = 90 features selected
# fit data to multinomial model
# accuracy = 0.4

#extract selected variables
p_val <- NULL
fit_train <- data.frame(train_meta, train_otu) #here use the train & test set for the whole data
fit_test <- data.frame(test_meta, test_otu)

for(i in c(1:2,4:ncol(fit_train))){
  multi_fit <- multinom(fit_train[,3]~fit_train[,i],  model = TRUE)
  aov_multi_fit <- Anova(multi_fit)
  p_val[i] <- aov_multi_fit$`Pr(>Chisq)`
}
p_val[3]=2
oo_fit = order(p_val)
select_var_fit = data.frame(fit_train[,oo_fit[1:f_select]]) #get top ranked predictors


mod_mult <- multinom(fit_train[,3]~., data = select_var_fit, MaxNWts=20000) #use selected variables to fit multinomial model
pred_mult <- predict(mod_mult, newdata = fit_test) #model predicts using test set
(c_matr <- confusionMatrix(pred_mult, fit_test$MOD))

##---------------------------bootstrap--------------------------------
# 95% CI = (0.200, 0.575)
# sd(accuracy) = 0.09810584
# mean(accuracy) = 0.3767375
s=10000
accuracy=numeric(s)
for(i in 1:s)
{
  boot_sample <- sample(1:nrow(fit_train),replace=T)
  subres=fit_train[boot_sample,3]
  subpre=select_var_fit[boot_sample,]
  
  boot_sample1 <- sample(1:nrow(fit_test),replace=T)
  subtest=fit_test[boot_sample1,]
  
  mod_mult <- multinom(subres~.,data=subpre, MaxNWts=20000)
  pred_mult <- predict(mod_mult, newdata = subtest)
  #accuracy[i]=sum(as.vector(fit_test$MOD)==as.vector(pred_mult))/length(pred_mult)
  c_matr <- confusionMatrix(pred_mult,subtest$MOD)
  (accuracy[i] = c_matr$overall[1])
}
quantile(accuracy,c(0.025,0.975))
sd(accuracy)
mean(accuracy)
median(accuracy)




##=========Multinomial Wrapper Method====================

#omit missing values in variable data set  
na.omit(otu_dta) 
na.omit(meta_dta)

#seperate training(80) and testing sets for the whole data set
set.seed(120)
train <- sample(1:120,80)
train_meta <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_meta <- meta_dta[-train,]
test_otu <- otu_dta[-train,]


#combine meta and otu training and testing for CV

cross_train = data.frame(train_meta, train_otu)
cross_test = data.frame(test_meta,test_otu)
dta=data.frame(meta_dta,otu_dta)

null <- multinom(MOD~ 1,data=na.omit(cross_train),model = TRUE)
full <- multinom(MOD~.,data=na.omit(cross_train), model = TRUE, MaxNWts=4000)
step.both <- step(null, scope=list(lower=null, upper=full), direction="both")  

pred_mult <- predict(step.both, newdata =cross_test) #model predicts using test set
(c_matr <- confusionMatrix(pred_mult, cross_test$MOD))
#accuracy = 0.5
#95% CI : (0.338, 0.662)

##---------------------------bootstrap--------------------

s=10000
accuracy=numeric(s)
for(i in 1:s)
{
  dta_train = data.frame(na.omit(cross_train))
  dta_test = data.frame((na.omit(cross_test)))
  
  boot_sample <- sample(1:nrow(dta_train),replace=T)
  #subres=dta_train[boot_sample,3] #responses
  subpre=dta_train[boot_sample,] #predictors
  
  boot_sample1 <- sample(1:nrow(dta_test),replace=T)
  subtest=dta_test[boot_sample1,]
  
  mod_mult <- multinom(step.both$call$formula, 
                       data = subpre, model = TRUE)
  pred_mult <- predict(mod_mult, newdata = subtest)
  #accuracy[i]=sum(as.vector(boot.dta$MOD)==as.vector(pred_mult))/length(pred_mult)
  c_matr <- confusionMatrix(pred_mult,subtest$MOD)
  accuracy[i] = c_matr$overall[1]
}
quantile(accuracy,c(0.025,0.975))
sd(accuracy)
mean(accuracy)
median(accuracy)
#CI=(0.30,  0.65)
# sd(accuracy)=0.090887
# mean(accuracy)= 0.04744741

##==============Random Forest Filter===============
#omit missing values in variable data set  
na.omit(otu_dta) 
na.omit(meta_dta)

#seperate training(80) and testing sets for the whole data set
set.seed(120)
train <- sample(1:120,80)
train_meta <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_meta <- meta_dta[-train,]
test_otu <- otu_dta[-train,]

# Divide training set into 8 folds
fold_size <-sample(1:80,80)
fold_1 <- fold_size[1:10]
fold_2 <- fold_size[11:20]
fold_3 <- fold_size[21:30]
fold_4 <- fold_size[31:40]
fold_5 <- fold_size[41:50]
fold_6 <- fold_size[51:60]
fold_7 <- fold_size[61:70]
fold_8 <- fold_size[71:80]
fold_samples = list(fold_1, fold_2, fold_3, fold_4, fold_5, fold_6, fold_7, fold_8)



## Now do CV. For each CV fold, we'll build and assess the accuracy of feature set sizes 
f_sizes = c(seq(10, 926, by = 10), 926) #yields 93 numbers
acc_f_b = matrix(NA, nrow = length(f_sizes), ncol = 8) #93*8


for(b in 1:8) { #loop for each fold
  cat(".")
  # meta & otu training and testing set for each fold
  fold_test_meta = train_meta[fold_samples[[b]],]
  fold_test_otu = train_otu[fold_samples[[b]],] 
  fold_train_meta = train_meta[-fold_samples[[b]],]
  fold_train_otu = train_otu[-fold_samples[[b]],] 
  
  #combine meta and otu training and testing for CV
  cross_train = data.frame(fold_train_meta, fold_train_otu)
  cross_test = data.frame(fold_test_meta,fold_test_otu)
  
  
  #Feature selection
  rf_mod = randomForest(MOD~., data=na.omit(cross_train), MaxNWts=4000)
  oo = order(importance(rf_mod), decreasing = T)
  
  for(f in 1:length(f_sizes)) {
    feature_set = oo[1:f_sizes[f]]
    
    #data frame of selected variables with top p-values
    select_train = na.omit(data.frame(fold_train_meta,fold_train_otu)[,c(3,feature_set)])
    select_test = na.omit(data.frame(fold_test_meta,fold_test_otu)[,c(3,feature_set)])
    
    fit_mod <- randomForest(MOD ~ ., data = select_train, MaxNWts=15000) #fit model using selected predictors 
    pred_mod <- predict(fit_mod, newdata = select_test)
    
    c_matr <- confusionMatrix(pred_mod, select_test$MOD) #one confusion matrix for one f
    
    acc_f_b[f, b] = c_matr$overall[1] #93*8 matrix
  }
}
acc_f = rowMeans(acc_f_b) #average accuracy for each f of all 8 folds 
f_select= f_sizes[which.max(acc_f)] #find maximum accuracy location
plot(f_sizes, acc_f, type = "both", xlab='feature size (f_sizes)', ylab='accuracy (acc_f)', 
     main='Feature size vs. Accuracy in Multinom (Wrapper) ')
#580 features selected

##------------------------extract selected variables and fit model----------------------------------------

# f = 58
# fit data to multinomial model

#extract selected variables
fit_train <- na.omit(data.frame(train_meta, train_otu)) #here use the train & test set for the whole data
fit_test <- na.omit(data.frame(test_meta, test_otu))
#dta=data.frame(meta_dta,otu_dta)

mod_rf = randomForest(MOD~., data=fit_train, MaxNWts=20000)
oo_fit = rownames(importance(mod_rf))[order(importance(mod_rf),decreasing=T)]

#select_var_fit = data.frame(dta[,oo_fit[1:f_select]])
rf_train = na.omit(fit_train[,c("MOD",oo_fit[1:f_select])])
rf_test = na.omit(fit_test[,c("MOD",oo_fit[1:f_select])])

mod_rf <- randomForest(MOD~., data = rf_train, MaxNWts=20000) #use selected variables to fit multinomial model
pred_rf <- predict(mod_rf, newdata = rf_test) #model predicts using test set
(c_matr <- confusionMatrix(pred_rf, rf_test$MOD))
# accuracy = 0.625, CI=(0.458, 0.7727)

##---------------------------bootstrap--------------------

s=10000
set.seed(120)
accuracy=numeric(s)
for(i in 1:s)
{
  boot_sample <- sample(1:nrow(rf_train),replace=T) #train data from selected varaibles
  #subres=rf_train[boot_sample,3] #responses
  subpre=rf_train[boot_sample,] #predictors
  boot_sample1 <- sample(1:nrow(rf_test),replace=T)
  subtest=rf_test[boot_sample1,]
  
  subpre$MOD=factor(subpre$MOD)
  subtest$MOD=factor(subtest$MOD)
  
  mod_rf <- randomForest(MOD~.,data=subpre, MaxNWts=20000)
  pred_rf <- predict(mod_rf, newdata = subtest)
  #c_matr <- confusionMatrix(pred_rf,subtest$MOD)
  #(accuracy[i] = c_matr$overall[1])
  accuracy[i]=sum(as.vector(subtest$MOD)==as.vector(pred_rf))/length(pred_rf)
}
quantile(accuracy,c(0.025,0.975))
sd(accuracy)
mean(accuracy)
median(accuracy)

#CI=(0.30, 0.65), sd=0.088629, mean(accuracy)=0.472575
