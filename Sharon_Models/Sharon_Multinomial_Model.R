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

##-------------------------------------------------------------------------
na.omit(otu_dta)
na.omit(meta_dta)

#seperate training and testing sets for the whole data set
set.seed(120)
train <- sample(1:120,80)
train_meta <- meta_dta[train,]
train_otu <-  otu_dta[train,]
test_meta <- meta_dta[-train,]
test_otu <- otu_dta[-train,]




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


# Cross Validation
f_sizes = c(seq(10, 926, by = 10), 926)
acc_f_b = matrix(NA, nrow = length(f_sizes), ncol = 8) # Misclassification Rate


for(b in 1:8) { #loop for each fold
  cat(".")
  fold_test_meta = train_meta[fold_samples[[b]],]
  fold_test_otu = train_otu[fold_samples[[b]],] 
  fold_train_meta = train_meta[-fold_samples[[b]],]
  fold_train_otu = train_otu[-fold_samples[[b]],] 
  
  cross_train = data.frame(fold_train_meta, fold_train_otu)
  cross_test = data.frame(fold_test_meta,fold_test_otu)

  ## Feature selection: rank p-value of OTU 
  # 

  p_val_otu <- NULL
  for(i in c(1:2,4:ncol(cross_train))){
    mod_fit <- multinom(cross_train[,3]~cross_train[,i],  model = TRUE)
    aov_mod_fit <- Anova(mod_fit)
    p_val_otu[i] <- aov_mod_fit$`Pr(>Chisq)`
  }
  p_val_otu[3]=2
  oo = order(p_val_otu) #rank p-values of OTU, ASCENDING
  
  for(f in 1:length(f_sizes)) {
    feature_set = oo[1:f_sizes[f]]

    select_train = data.frame(fold_train_meta,fold_train_otu)[,c(3,feature_set)]
    select_test = data.frame(fold_test_meta,fold_test_otu)[,c(3,feature_set)]
    
    fit_mod <- multinom(MOD ~ ., data = select_train, MaxNWts=15000) #fit model using p-value-ranked OTU 
    pred_mod <- predict(fit_mod, newdata = select_test)
    
    c_matr <- confusionMatrix(pred_mod, select_test$MOD) #one confusion matrix for one f
    acc_f_b[f, b] = c_matr$overall[1]
    
  }
}
acc_f = rowMeans(acc_f_b) #accuracy for each fold
which.max(acc_f) #find maximum accuracy location

#----------------------------------------------------------------
# f = 9 --> 90 features selected
# fit data to multinomial model

#Get selected variables
p_val <- NULL
fit_train <- data.frame(train_meta, train_otu)
fit_test <- data.frame(test_meta, test_otu)ma

for(i in c(1:2,4:ncol(fit_train))){
  multi_fit <- multinom(fit_train[,3]~fit_train[,i],  model = TRUE)
  aov_multi_fit <- Anova(multi_fit)
  p_val[i] <- aov_multi_fit$`Pr(>Chisq)`
}
p_val[3]=2
oo_fit = order(p_val)
select_var_fit = data.frame(fit_train[,oo_fit[1:90]])


mod_mult <- multinom(fit_train[,3]~., data = select_var_fit, MaxNWts=20000) #use selected variables to fit model
pred_mult <- predict(mod_mult, newdata = fit_test) #model predicts using test set
(c_matr <- confusionMatrix(pred_mult, fit_test$MOD))

