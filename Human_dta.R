# Install Phyloseq package:
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

# Loading package: reference: https://joey711.github.io/phyloseq/install.html
library(phyloseq)

# Set directory
setwd("/Users/Le/Google Drive/Research/Human Remain Data")

# Load-in the OTU table
dta<-import_biom("HPMM_TAMU_Collab/f120.biom")

# Load-in the meta-data
meta_dta<-read.delim2("HPMM_TAMU_Collab/map_f120.txt",na.strings="na")
row.names(meta_dta)<-meta_dta$X.SampleID
meta_dta$BMI<-as.numeric(meta_dta$BMI)

# Integrate the metadata into OTU data
sample_data(dta)<-sample_data(meta_dta)


body_meta<-meta_dta[!duplicated(meta_dta$Pack_ID),]
attach(body_meta)
table(Sex)
table(Race)
table(Event_Location)
table(Season)
table(Estimated_PMI)
table(Weight_Status)
table(Season)
table(MoD)



summary(Age)
summary(BMI)

table(Sex, Race)
table(Sex, Event_Location)
table(Sex, Season)
table(Sex, Estimated_PMI)
table(Sex, MoD)
table(Sex, Weight_Status, useNA="always")
par(mfrow=c(1,2))
boxplot(Age~Sex,main="Age")
boxplot(BMI~Sex,main="BMI")

table(Race, Event_Location, useNA="always")
table(Race, Season, useNA="always")
table(Race, Estimated_PMI, useNA="always")
table(Race, MoD, useNA="always")
table(Race, Weight_Status, useNA="always")
par(mfrow=c(1,2))
boxplot(Age~ Race,main="Age")
boxplot(BMI~ Race,main="BMI")


table(Event_Location, Season)
table(Event_Location, Estimated_PMI)
table(Event_Location, MoD)
table(Event_Location, Weight_Status, useNA="always")
par(mfrow=c(1,2))
boxplot(Age~ Event_Location,main="Age")
boxplot(BMI~ Event_Location,main="BMI")



table(Estimated_PMI, MoD)
table(Estimated_PMI, Season)
table(Estimated_PMI, Weight_Status, useNA="always")
par(mfrow=c(1,2))
boxplot(Age~ Estimated_PMI,main="Age")
boxplot(BMI~ Estimated_PMI,main="BMI")

table(Season, MoD)
table(Season, Weight_Status, useNA="always")
par(mfrow=c(1,2))
boxplot(Age~ Season,main="Age")
boxplot(BMI~ Season,main="BMI")


table(MoD, Weight_Status, useNA="always")
par(mfrow=c(1,2))
boxplot(Age~ MoD,main="Age")
boxplot(BMI~ MoD,main="BMI")

par(mfrow=c(1,2),mar=c(8,6,4,4))
boxplot(Age~ Weight_Status,main="Age", las=2)
boxplot(BMI~ Weight_Status,main="BMI",las=2)

detach(body_meta)



# Use otu_table() function to access the OTU data
otu_raw<-as.matrix(otu_table(dta))


body<-levels(meta_dta$Pack_ID)

sum_otu<-matrix(NA, nrow=dim(otu_raw)[1],ncol=length(body))
colnames(sum_otu)<-body
row.names(sum_otu)<-rownames(otu_raw)

for (i in body){
	print(i)
	temp<-meta_dta[meta_dta$Pack_ID%in%i,]$X.SampleID
	temp_2<-rowMeans(otu_raw[,temp])
	sum_otu[,i]<-temp_2
	
}

zero_count<-data.frame(Pack_ID=body,count=rep(NA, length(body)))
row.names(zero_count)<-body
for (i in body){
	zero_count[i,2]<-sum(sum_otu[,i]==0)/length(sum_otu[,i])
}



sum_matrix<-t(t(colSums(sum_otu)))
sum_matrix<-data.frame(sum_matrix,Pack_ID=rownames(sum_matrix))

merged_meta<-merge(body_meta, sum_matrix, by="Pack_ID")
merged_meta_zero<-merge(merged_meta,zero_count, by="Pack_ID")
dim(merged_meta_zero)
colnames(merged_meta_zero)[c(17,18)]<-c("OTU_sum","Zero_prop")

attach(merged_meta_zero)

par(mfrow=c(1,2))
boxplot(OTU_sum~Sex,main="OTU_sum")
boxplot(Zero_prop~Sex,main="Zero_Prop")

par(mfrow=c(1,2))
boxplot(OTU_sum~Race,main="OTU_sum")
boxplot(Zero_prop~Race,main="Zero_Prop")

par(mfrow=c(1,2))
boxplot(OTU_sum~Event_Location,main="OTU_sum")
boxplot(Zero_prop~ Event_Location,main="Zero_Prop")

par(mfrow=c(1,2))
boxplot(OTU_sum~Estimated_PMI,main="OTU_sum")
boxplot(Zero_prop~ Estimated_PMI,main="Zero_Prop")

par(mfrow=c(1,2))
boxplot(OTU_sum~Season,main="OTU_sum")
boxplot(Zero_prop~ Season,main="Zero_Prop")

par(mfrow=c(1,2))
boxplot(OTU_sum~MoD,main="OTU_sum")
boxplot(Zero_prop~ MoD,main="Zero_Prop")

par(mfrow=c(1,2),mar=c(10,6,4,4))
boxplot(OTU_sum~Weight_Status,main="OTU_sum",las=2)
boxplot(Zero_prop~Weight_Status,main="Zero_Prop",las=2)

detach(merged_meta_zero)

temp_sum<-cbind(
apply(sum_otu,2,min),
apply(sum_otu,2,max),
apply(sum_otu,2,median))
temp_sum<-data.frame(Pack_ID=rownames(temp_sum),temp_sum)
colnames(temp_sum)<-c("Pack_ID","Min","Max","Median")

sum_merged_meta_zero<-merge(merged_meta_zero,temp_sum,by="Pack_ID")

attach(sum_merged_meta_zero)
par(mfrow=c(2,2))
plot(Min~BMI,main="Min OTU")
plot(Max~BMI,main="Max OTU")
plot(Median~BMI,main="Median OTU")
plot(Zero_prop~BMI,main="Zero_prop OTU")

par(mfrow=c(2,2))
plot(Min~Age,main="Min OTU")
plot(Max~Age,main="Max OTU")
plot(Median~Age,main="Median OTU")
plot(Zero_prop~Age,main="Zero_prop OTU")





# Use sample_data() function to access the metadata
sample_data(dta)

#Use the tax_table() function to access the taxaonomy
tax_table(dta)

