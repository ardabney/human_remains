# Change Based on User's Preference 
setwd("~/Desktop")

#load in needed libraries
library(phyloseq)
library("vegan")

#import data
std_dta  <- import_biom("f120_r1k.biom") #whole data
m_dta <- data.frame(sample_data(std_dta)) #meta data
OTU_dta <- data.frame(otu_table(std_dta)) #OTU data
t_dta <- data.frame(tax_table(std_dta)) #Taxa data

# Meta data EDA - Split By Individual

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

# Adjust Age & BMI into groups/factors for easier analysis 
m_dta$Age_Group <- cut(m_dta$Age,c(0,25,50,75,100))
m_dta$BMI_Group <- cut(m_dta$BMI,c(0,10,20,30,40,50,60))

# Create data Set - one with all the patients' samples combined so we can explore Race, Sex, etc,
mt_dta<-m_dta[!duplicated(m_dta$Pack_ID),] #patients
row.names(mt_dta)<-mt_dta$Pack_ID
meta_dta <-subset(mt_dta,select=c("Estimated_PMI","Race","Manner.of.Death","Season","Sex","Weight_Status","Event_Location","BMI_Group","Age_Group"))

# Boxplots - PMI - Split by Person (This is still meaningful, but the tables and barplots are probably
# better for looking at the categorical data)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Race)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Season)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Sex)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Weight_Status)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Event_Location)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$Age_Group)
boxplot(meta_dta$Estimated_PMI ~ meta_dta$BMI_Group)

#Barplots/Tables - PMI - Split by Person 
(PMI_dta <- table(meta_dta$Estimated_PMI, meta_dta$Race))
barplot(PMI_dta, beside=T, main="PMI vs. Race", xlab = "Race", ylab="Counts",col=c("lightgreen","white","grey","slategray1","pink"),legend = rownames(PMI_dta))
(PMI_dta <- table(meta_dta$Estimated_PMI, meta_dta$Season))
barplot(PMI_dta, beside=T, main="PMI vs. Season", xlab = "Season", ylab="Counts",col=c("lightgreen","white","grey","slategray1","pink"),legend = rownames(PMI_dta))
(PMI_dta <- table(meta_dta$Estimated_PMI, meta_dta$Sex))
barplot(PMI_dta, beside=T, main="PMI vs. Sex",xlab = "Sex", ylab="Counts",col=c("lightgreen","white","grey","slategray1","pink"),legend = rownames(PMI_dta))
(PMI_dta <- table(meta_dta$Estimated_PMI, meta_dta$Weight_Status))
barplot(PMI_dta, beside=T, main="PMI vs. Weight Status", xlab = "Weight Status", ylab="Counts",col=c("lightgreen","white","grey","slategray1","pink"),legend = rownames(PMI_dta))
(PMI_dta <- table(meta_dta$Estimated_PMI, meta_dta$Event_Location))
barplot(PMI_dta, beside=T, main="PMI vs. Event Location", xlab = "Event Location", ylab="Counts",col=c("lightgreen","white","grey","slategray1","pink"),legend = rownames(PMI_dta))
(PMI_dta <- table(meta_dta$Estimated_PMI, meta_dta$Age_Group))
barplot(PMI_dta, beside=T, main="PMI vs. Age", xlab = "Age Group", ylab="Counts",col=c("lightgreen","white","grey","slategray1","pink"),legend = rownames(PMI_dta))
(PMI_dta <- table(meta_dta$Estimated_PMI, meta_dta$BMI_Group))
barplot(PMI_dta, beside=T, main="PMI vs. BMI", xlab = "BMI Group", ylab="Counts",col=c("lightgreen","white","grey","slategray1","pink"),legend = rownames(PMI_dta))

#Barplots/Tables - MOD - Split by Person 
(MOD_dta <- table(meta_dta$Manner.of.Death, meta_dta$Race))
barplot(MOD_dta, beside=T, main="MOD vs. Race", xlab = "Race", ylab="Counts",col=c("lightgreen","white","grey","slategray1"),legend = rownames(MOD_dta))
(MOD_dta <- table(meta_dta$Manner.of.Death, meta_dta$Season))
barplot(MOD_dta, beside=T, main="MOD vs. Season", xlab = "Season", ylab="Counts",col=c("lightgreen","white","grey","slategray1"),legend = rownames(MOD_dta))
(MOD_dta <- table(meta_dta$Manner.of.Death, meta_dta$Sex))
barplot(MOD_dta, beside=T, main="MOD vs. Sex",xlab = "Sex", ylab="Counts",col=c("lightgreen","white","grey","slategray1"),legend = rownames(MOD_dta))
(MOD_dta <- table(meta_dta$Manner.of.Death, meta_dta$Weight_Status))
barplot(MOD_dta, beside=T, main="MOD vs. Weight Status", xlab = "Weight Status", ylab="Counts",col=c("lightgreen","white","grey","slategray1"),legend = rownames(MOD_dta))
(MOD_dta <- table(meta_dta$Manner.of.Death, meta_dta$Event_Location))
barplot(MOD_dta, beside=T, main="MOD vs. Event Location", xlab = "Event Location", ylab="Counts",col=c("lightgreen","white","grey","slategray1"),legend = rownames(MOD_dta))
(MOD_dta <- table(meta_dta$Manner.of.Death, meta_dta$Age_Group))
barplot(MOD_dta, beside=T, main="MOD vs. Age", xlab = "Age Group", ylab="Counts",col=c("lightgreen","white","grey","slategray1"),legend = rownames(MOD_dta))
(MOD_dta <- table(meta_dta$Manner.of.Death, meta_dta$BMI_Group))
barplot(MOD_dta, beside=T, main="MOD vs. BMI", xlab = "BMI Group", ylab="Counts",col=c("lightgreen","white","grey","slategray1"),legend = rownames(MOD_dta))

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

#Note: This following code for sample areas is not updated because it did not appear very meaningful
# orginally, but I'm leaving it in the code in case we later decide to use it, or come up with another
# method for analysing it. 

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

