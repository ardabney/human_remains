# Change Based on User's Preference 
setwd("~/Desktop")

#load in needed libraries
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)
library("vegan")

#import data
std_dta  <- import_biom("f120_r1k.biom") #whole data
m_dta <- data.frame(sample_data(std_dta)) #meta data
OTU_dta <- data.frame(otu_table(std_dta)) #OTU data
t_dta <- data.frame(tax_table(std_dta)) #Taxa data

# Meta data EDA
layout(matrix(c(1,2,3,4,5,6,7,7), 2, 4, byrow = TRUE))

#Adjust variables into correct formats
m_dta$BMI<-as.numeric(m_dta$BMI)
m_dta$Age<-as.numeric(m_dta$Age)
m_dta$Estimated_PMI<-factor(m_dta$Estimated_PMI,levels=c("12","<24",">24",">48",">72"),ordered=T)
m_dta$Race<-factor(m_dta$Race)
m_dta$Manner.of.Death<-factor(m_dta$Manner.of.Death,levels=c("Natural","Accident","Suicide","Homicide"),ordered=T)
m_dta$Season<-factor(m_dta$Season,levels=c("Spring","Summer","Autumn","Winter"),ordered=T)
m_dta$Sex<-factor(m_dta$Sex)
m_dta$Weight_Status<-factor(m_dta$Weight_Status,levels=c("Underweight","Normal Weight","Overweight","Obese","Severe Obesity","Morbid Obesity","Super Obese"))
m_dta$Event_Location<-factor(m_dta$Event_Location)
m_dta$Sample_Area <- factor(m_dta$Sample_Area,levels=c("Buccal","Ears","Eyes","Nares", "Rectum", "Umbilicus"))

# Create data Set - one with all the patients' samples combined so we can explore Race, Sex, etc,
mt_dta<-m_dta[!duplicated(m_dta$Pack_ID),] #patients
row.names(mt_dta)<-mt_dta$Pack_ID
meta_dta <-subset(mt_dta,select=c("Estimated_PMI","Race","Manner.of.Death","Season","Sex","Weight_Status","Event_Location","BMI","Age"))

# Boxplots - PMI - Split by Person
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Race)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Season)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Sex)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Weight_Status)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Event_Location)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Age)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$BMI)

#Boxplots - MOD - Split by Person 
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Race)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Season)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Sex)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Weight_Status)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Event_Location)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Age)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$BMI)

# Summary Statistics - Split by Person
summary(meta_dta)

layout(matrix(c(1), 1, 1, byrow = TRUE))
# Clustering - Split by Person
hc_ind <- hclust(dist(meta_dta))
plot(hc_ind)

# PCA - Split by Person
pca <- prcomp(dist(meta_dta))
summary(pca)
biplot(pca)

layout(matrix(c(1,2,3,4,5,6,7,7), 2, 4, byrow = TRUE))
#Sample Area Analysis
meta_dta <-m_dta[m_dta$Sample_Area == "Umbilicus",] 

# Boxplots - PMI - Umbilicus
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Race)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Season)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Sex)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Weight_Status)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Event_Location)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Age)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$BMI)

#Boxplots - MOD - Umbilicus
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Race)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Season)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Sex)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Weight_Status)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Event_Location)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Age)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$BMI)

# Summary Statistics - Umbilicus
summary(meta_dta)

meta_dta <-m_dta[m_dta$Sample_Area == "Nares",] 
# Boxplots - PMI - Nares
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Race)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Season)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Sex)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Weight_Status)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Event_Location)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Age)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$BMI)

#Boxplots - MOD - Nares
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Race)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Season)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Sex)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Weight_Status)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Event_Location)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Age)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$BMI)

# Summary Statistics - Nares
summary(meta_dta)

meta_dta <-m_dta[m_dta$Sample_Area == "Ears",] 
# Boxplots - PMI - Ears
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Race)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Season)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Sex)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Weight_Status)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Event_Location)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Age)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$BMI)

#Boxplots - MOD - Ears
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Race)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Season)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Sex)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Weight_Status)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Event_Location)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Age)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$BMI)

# Summary Statistics - Ears
summary(meta_dta)

meta_dta <-m_dta[m_dta$Sample_Area == "Eyes",] 
# Boxplots - PMI - Eyes
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Race)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Season)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Sex)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Weight_Status)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Event_Location)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Age)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$BMI)

#Boxplots - MOD -Eyes
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Race)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Season)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Sex)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Weight_Status)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Event_Location)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Age)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$BMI)

# Summary Statistics - Eyes
summary(meta_dta)

meta_dta  <- m_dta[m_dta$Sample_Area == "Buccal",]
# Boxplots - PMI - Buccal
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Race)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Season)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Sex)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Weight_Status)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Event_Location)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Age)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$BMI)

#Boxplots - MOD - Buccal
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Race)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Season)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Sex)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Weight_Status)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Event_Location)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Age)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$BMI)

# Summary Statistics - Buccal
summary(meta_dta)

meta_dta  <- m_dta[m_dta$Sample_Area == "Rectum",]
# Boxplots - PMI - Rectum
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Race)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Season)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Sex)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Weight_Status)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Event_Location)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Age)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$BMI)

#Boxplots - MOD - Rectum
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Race)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Season)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Sex)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Weight_Status)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Event_Location)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$Age)
boxplot(meta_dta$Manner.of.Death ~ meta_dta$BMI)

# Summary Statistics - Rectum
summary(meta_dta)

