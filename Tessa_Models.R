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
otu_dta <-t(data.frame(ave_otu_dta))

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
train<-sample(1:120,round(120*0.66))
train_mt<-meta_dta[train,]
train_otu<-select_otu[train,]
test_mt <- meta_dta[-train,]
test_otu <- select_otu[-train,]

# Model creation and Testing
training <- data.frame(train_mt,train_otu)
pmi_mult <- multinom(Estimated_PMI ~ ., data = training, MaxNWts=1500)
testing <- data.frame(test_mt,test_otu)
pred_pmi <- predict(pmi_mult, newdata = testing)
table(pred_pmi, testing$Estimated_PMI) #Note: I'm not sure why 12 appears on this table
confusionMatrix(pred_pmi, testing$Estimated_PMI)
