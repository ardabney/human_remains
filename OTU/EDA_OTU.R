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

#Average OTU data
ave_otu_dta<-matrix(NA,nrow=dim(OTU_dta)[1],ncol=length(unique(m_dta$Pack_ID)))
row.names(ave_otu_dta)<-rownames(OTU_dta)
colnames(ave_otu_dta)<-unique(m_dta$Pack_ID)
for (i in unique(m_dta$Pack_ID)){
  temp_id<-rownames(m_dta[m_dta$Pack_ID%in%i,])
  ave_otu_dta[,i]<-rowMeans(OTU_dta[,temp_id])
}

#Average values of the columns(sample IDs) in the OTU data
otu_matrix<-data.matrix(ave_otu_dta)
colMeans(otu_matrix)
colMins(otu_matrix) #min of column values
colMaxs(otu_matrix) #max of column values

#EDA (Tentative - Rarecurve/Diversity,Clustering, PCA)
otu_dta <- round(as.matrix(ave_otu_dta))
rarecurve(otu_dta)
diversity(otu_dta, index="shannon")
diversity(otu_dta, index="simpson")

#Clustering
hc_ind <- hclust(dist(ave_otu_dta))
plot(hc_ind)

#PCA
pca <- prcomp(dist(ave_otu_dta))
summary(pca)
biplot(pca)
