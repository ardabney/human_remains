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
    select_train = data.frame(fold_train_meta,fold_train_otu)[,c(3,feature_set)] 
    select_test = data.frame(fold_test_meta,fold_test_otu)[,c(3,feature_set)]
    
    fit_mod <- randomForest(MOD ~ ., data = na.omit(select_train), MaxNWts=15000) #fit model using selected predictors 
    pred_mod <- predict(fit_mod, newdata = na.omit(select_test))
    
    c_matr <- confusionMatrix(pred_mod, na.omit(select_test)$MOD) #one confusion matrix for one f
  
    acc_f_b[f, b] = c_matr$overall[1] #93*8 matrix
  }
}
acc_f = rowMeans(acc_f_b) #average accuracy for each f of all 8 folds 
which.max(acc_f) #find maximum accuracy location --> location(f) * 10 = how many features selected
#660 features selected

##------------------------extract selected variables and fit model----------------------------------------

# f = 66 
# fit data to multinomial model

#extract selected variables
fit_train <- na.omit(data.frame(train_meta, train_otu)) #here use the train & test set for the whole data
fit_test <- na.omit(data.frame(test_meta, test_otu))
dta=data.frame(meta_dta,otu_dta)

mod_rf = randomForest(MOD~., data=na.omit(fit_train), MaxNWts=20000)
oo_fit = order(importance(mod_rf),decreasing=T)

select_var_fit = data.frame(dta[,oo_fit[1:66]])
rf_train = fit_train[,c(3,oo_fit[1:66])]
rf_test = fit_test[,c(3,oo_fit[1:66])]

mod_rf <- randomForest(MOD~., data = rf_train, MaxNWts=20000) #use selected variables to fit multinomial model
pred_rf <- predict(mod_rf, newdata = rf_test) #model predicts using test set
(c_matr <- confusionMatrix(pred_rf, fit_test$MOD))
# accuracy = 0.525
##---------------------------bootstrap--------------------

s=100
accuracy=numeric(s)
for(i in 1:s)
{
  boot_sample <- sample(1:nrow(rf_train),replace=T) #train data from selected varaibles
  #subres=rf_train[boot_sample,3] #responses
  subpre=rf_train[boot_sample,] #predictors
  mod_rf <- randomForest(MOD~.,data=na.omit(subpre), MaxNWts=20000)
  pred_rf <- predict(mod_rf, newdata = rf_test)
  #accuracy[i]=sum(as.vector(fit_test$MOD)==as.vector(pred_mult))/length(pred_mult)
  c_matr <- confusionMatrix(pred_rf,na.omit(rf_test)$MOD)
  accuracy[i] = c_matr$overall[1]
}
quantile(accuracy,c(0.025,0.975))
sd(accuracy)
#CI = 0.36,0.575, s= 0.2495071