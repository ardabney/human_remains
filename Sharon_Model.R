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


# select best predictors using forward and backwar
for(i in 1:926){
  selected <- data.frame(meta_dta,otu_dta[,i])
}
#fix missing values
selected2 <- na.omit(selected)

null <- multinom(MOD~ 1,selected2,model = TRUE)
full <- multinom(MOD~.,selected2, model = TRUE)
step.fwd <- step(null, scope=list(lower=null, upper=full), direction="forward")
step.bwd <- step(full, scope=list(lower=null, upper=full), direction="backward")

#both left only two predictors: age_group + race
step.fwd$anova 
step.bwd$anova
#fwd_select <- step.fwd$anova$Step


#Creating Test and Training Data
#I only dealt with meta data because only its predictos have been used 
set.seed(300)
train<-sample(1:120,80)
train_mt<-meta_dta[train,]
test_mt <- meta_dta[-train,]


# Multinomial Model creation and Testing
training <- data.frame(train_mt)
#porblem: code not robusted, I manually inputed the x variables, could be improved 
mod_mult <- multinom(MOD ~ Age_Group + Race, data = training, MaxNWts=1500)

testing <- data.frame(test_mt)
pred_mod <- predict(mod_mult, newdata = testing)
confusionMatrix(pred_mod, testing$MOD)
