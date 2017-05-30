# Change Based on User's Preference 
setwd("~/Desktop")

#load in needed libraries
library(phyloseq)

#import data
std_dta  <- import_biom("f120_r1k.biom") #whole data
m_dta <- data.frame(sample_data(std_dta)) #meta data
OTU_dta <- data.frame(otu_table(std_dta)) #OTU data
t_dta <- data.frame(tax_table(std_dta)) #Taxa data

# Adjust variables into correct formats
m_dta$BMI <- as.numeric(m_dta$BMI)
m_dta$Age <- as.numeric(m_dta$Age)
m_dta$Estimated_PMI <- factor(m_dta$Estimated_PMI,levels=c("12","<24",">24",">48",">72"), ordered=T)
m_dta$Race <- factor(m_dta$Race)
m_dta$Manner.of.Death <- factor(m_dta$Manner.of.Death,levels=c("Natural","Accident","Suicide","Homicide"),ordered=T)
m_dta$Season <- factor(m_dta$Season,levels=c("Spring","Summer","Autumn","Winter"),ordered=T)
m_dta$Sex <- factor(m_dta$Sex)
m_dta$Weight_Status <- factor(m_dta$Weight_Status,levels=c("Underweight","Normal Weight","Overweight","Obese","Severe Obesity","Morbid Obesity","Super Obese"))
m_dta$Event_Location <- factor(m_dta$Event_Location)
m_dta$Sample_Area <- factor(m_dta$Sample_Area,levels=c("Buccal","Ears","Eyes","Nares", "Rectum", "Umbilicus"))

# Data Processing - Meta
mt_dta<-m_dta[!duplicated(m_dta$Pack_ID),] #patients
row.names(mt_dta)<-mt_dta$Pack_ID
meta_dta <-subset(mt_dta,select=c("Estimated_PMI","Race","Manner.of.Death","Season","Sex","Weight_Status","Event_Location","BMI_Group","Age_Group"))

# Data Processing - Taxa
colnames(t_dta) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
t_dta$Kingdom <- gsub("k__", "", t_dta$Kingdom)
t_dta$Phylum <- gsub("p__", "", t_dta$Phylum)
t_dta$Class <- gsub("c__", "", t_dta$Class)
t_dta$Order <- gsub("o__", "", t_dta$Order)
t_dta$Family <- gsub("f__", "", t_dta$Family)
t_dta$Genus <- gsub("g__", "", t_dta$Genus)
t_dta$Species <- gsub("s__", "", t_dta$Species)