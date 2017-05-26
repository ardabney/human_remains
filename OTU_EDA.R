#Import data

library(phyloseq)
library("vegan")
library(matrixStats)

std_dta <- import_biom("f120_r1k.biom") #access the whole data
meta_dta <- data.frame(sample_data(std_dta)) #Meta data
otu_dta <- data.frame(otu_table(std_dta)) #OTU data
taxa_dta <- data.frame(tax_table(std_dta)) #Taxa data

rarecurve(otu_dta)
diversity(otu_dta, index="shannon")
diversity(otu_dta, index="simpson")

#Average OTU data
ave_otu_dta<-matrix(NA,nrow=dim(otu_dta)[1],ncol=length(unique(meta_dta$Pack_ID)))
row.names(ave_otu_dta)<-rownames(otu_dta)
colnames(ave_otu_dta)<-unique(metaDta$Pack_ID)
for (i in unique(meta_dta$Pack_ID)){
	print(i)
	temp_id<-rownames(meta_dta[meta_dta$Pack_ID%in%i,])
	ave_otu_dta[,i]<-rowMeans(otu_dta[,temp_id])
}

#Average values of the columns(sample IDs) in the OTU data
otu_matrix<-data.matrix(ave_otu_dta)
colMeans(otu_matrix)
colMins(otu_matrix) #min of column values
colMaxs(otu_matrix) #max of column values
