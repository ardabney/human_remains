# Install Phyloseq package:
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

#####################################################
# Initial exploration

# Loading package: reference: https://joey711.github.io/phyloseq/install.html
library(phyloseq)
library(vegan)

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

rarecurve(round(t(as.matrix(sum_otu))))

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







#####################################################
# Standarized data for diversity

# Loading package: reference: https://joey711.github.io/phyloseq/install.html
library(phyloseq)
library(vegan)

# Set directory
setwd("/Users/Le/Google Drive/Research/Human Remain Data")

# Load in the std data
std_dta<-import_biom("/Users/Le/Google Drive/Research/Human Remain Data/HPMM_TAMU_Collab/f120_r1k.biom")

# Reconsturct the OTU file by taking the average OTU of each body over sample_area
std_meta_dta<-data.frame(sample_data(std_dta))
std_otu_dta<-data.frame(otu_table(std_dta))
std_taxa_dta<-data.frame(tax_table(std_dta))

## Average all otu table
ave_otu_dta<-matrix(NA,nrow=dim(std_otu_dta)[1],ncol=length(unique(std_meta_dta$Pack_ID)))
row.names(ave_otu_dta)<-rownames(std_otu_dta)
colnames(ave_otu_dta)<-unique(std_meta_dta$Pack_ID)
for (i in unique(std_meta_dta$Pack_ID)){
	print(i)
	temp_id<-rownames(std_meta_dta[std_meta_dta$Pack_ID%in%i,])
	ave_otu_dta[,i]<-rowMeans(std_otu_dta[,temp_id])
}

## Adjust the meta data
ave_meta_dta<-std_meta_dta[!duplicated(std_meta_dta$Pack_ID),]
row.names(ave_meta_dta)<-ave_meta_dta$Pack_ID
ave_meta_dta$BMI<-as.numeric(ave_meta_dta$BMI)
ave_meta_dta$Age<-as.numeric(ave_meta_dta$Age)

## Adjust the taxa data
ave_taxa_dta<-data.frame(matrix(NA, nrow=dim(std_taxa_dta)[1],ncol=dim(std_taxa_dta)[2]))
colnames(ave_taxa_dta)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Speices")
row.names(ave_taxa_dta)<-rownames(std_taxa_dta)
ave_taxa_dta$Kingdom<-gsub("k__","",std_taxa_dta$Rank1)
ave_taxa_dta$Phylum<-gsub("p__","",std_taxa_dta$Rank2)
ave_taxa_dta$Class<-gsub("c__","",std_taxa_dta$Rank3)
ave_taxa_dta$Order<-gsub("o__","",std_taxa_dta$Rank4)
ave_taxa_dta$Family<-gsub("f__","",std_taxa_dta$Rank5)
ave_taxa_dta$Genus<-gsub("g__","",std_taxa_dta$Rank6)
ave_taxa_dta$Speices<-gsub("s__","",std_taxa_dta$Rank7)

## Combine into a phyloseq format
ave_dta<-phyloseq(
	otu_table(round(as.matrix(ave_otu_dta)),taxa_are_rows=T),
	tax_table(as.matrix(ave_taxa_dta)),
	sample_data(ave_meta_dta)
	)
# rm 0 count otu
ave_dta_clean<-filter_taxa(ave_dta,function(x) sum(x)>0,T)

# Rarefication curve on OTU
rarecurve(t(round(as.matrix(otu_table(ave_dta)))))

# alpha diversity
alpha_temp<-data.frame(estimate_richness(ave_dta,measures=c("Simpson","Shannon")))
alpha_diversity<-data.frame(alpha_temp,Pack_ID=rownames(alpha_temp))
alpha_data<-merge(ave_meta_dta,alpha_diversity,by="Pack_ID")
alpha_data$BMI<-as.numeric(alpha_data$BMI)
alpha_data$Age<-as.numeric(alpha_data$Age)

# Plot alpha diversity by univariate
plotting_factors<-(c("Estimated_PMI","Race","Manner.of.Death","Season","Sex","Weight_Status","Event_Location"))
pdf("Alpha_diversity.pdf",width=24,height=6)
for(i in plotting_factors){
	par(mfrow=c(1,2))
	boxplot(alpha_data$Simpson~alpha_data[,i],main=paste(i,"with","Simpson"))
	boxplot(alpha_data$Shannon~alpha_data[,i],main=paste(i,"with","Shannon"))
	
}

par(mfrow=c(1,2))
plot(Simpson~BMI, data=alpha_data)
plot(Shannon~BMI, data=alpha_data)

par(mfrow=c(1,2))
plot(Simpson~Age, data=alpha_data)
plot(Shannon~Age, data=alpha_data)
dev.off()

# Plot bivariate alpha

pdf("Bivariate_alpha.pdf")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Estimated_PMI",color="Race")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Estimated_PMI",color="Manner.of.Death")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Estimated_PMI",color="Season")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Estimated_PMI",color="Sex")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Estimated_PMI",color="Weight_Status")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Estimated_PMI",color="Event_Location")


plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Race",color="Manner.of.Death")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Race",color="Season")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Race",color="Sex")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Race",color="Weight_Status")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Race",color="Event_Location")

plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Manner.of.Death",color="Season")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Manner.of.Death",color="Sex")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Manner.of.Death",color="Weight_Status")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Manner.of.Death",color="Event_Location")

plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Season",color="Sex")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Season",color="Weight_Status")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Season",color="Event_Location")

plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Sex",color="Weight_Status")
plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Sex",color="Event_Location")

plot_richness(ave_dta,measures=c("Simpson","Shannon"),x="Weight_Status",color="Event_Location")

dev.off()


# Plotting beta diversity


pdf("beta_CAP.pdf")
temp_ord<-ordinate(ave_dta,method="CAP",distance="bray",formula=~Estimated_PMI)
plot_ordination(ave_dta,temp_ord,type="samples",color="Estimated_PMI",title="Estimated_PMI")

temp_ord<-ordinate(ave_dta,method="CAP",distance="bray",formula=~Race)
plot_ordination(ave_dta,temp_ord,type="samples",color="Race",title="Race")

temp_ord<-ordinate(ave_dta,method="CAP",distance="bray",formula=~Manner.of.Death)
plot_ordination(ave_dta,temp_ord,type="samples",color="Manner.of.Death",title="Manner.of.Death")

temp_ord<-ordinate(ave_dta,method="CAP",distance="bray",formula=~Season)
plot_ordination(ave_dta,temp_ord,type="samples",color="Season",title="Season")

temp_ord<-ordinate(ave_dta,method="CAP",distance="bray",formula=~Sex)
plot_ordination(ave_dta,temp_ord,type="samples",color="Sex",title="Sex")

temp_ord<-ordinate(ave_dta,method="CAP",distance="bray",formula=~Weight_Status)
plot_ordination(ave_dta,temp_ord,type="samples",color="Weight_Status",title="Weight_Status")

temp_ord<-ordinate(ave_dta,method="CAP",distance="bray",formula=~BMI,na.action=na.omit)
p1=plot_ordination(ave_dta,temp_ord,type="samples",title="BMI")
p1+geom_point(aes(colour=BMI))+scale_colour_gradient(low="red",high="blue")

temp_ord<-ordinate(ave_dta,method="CAP",distance="bray",formula=~Age,na.action=na.omit)
p1=plot_ordination(ave_dta,temp_ord,type="samples",title="Age")
p1+geom_point(aes(colour=Age))+scale_colour_gradient(low="red",high="blue")

dev.off()



#####################################################
# Standarized data for hypothesis testing, Dec 18, 2016
# Le

# load-in library
library(phyloseq)
library(vegan)
library(DESeq2)

# Use data from previous steps
setwd("/Users/Le/Google Drive/Research/Human Remain Data")
ave_dta_clean


meta_testing<-subset(data.frame(sample_data(ave_dta_clean)),select=c("Estimated_PMI","Race","Manner.of.Death","Season","Sex","Weight_Status","Event_Location","BMI","Age"))
meta_testing$Estimated_PMI<-factor(meta_testing$Estimated_PMI,levels=c("12","<24",">24",">48",">72"),ordered=T)
meta_testing$Race<-factor(meta_testing$Race)
meta_testing$Manner.of.Death<-factor(meta_testing$Manner.of.Death,levels=c("Natural","Accident","Suicide","Homicide",ordered=T))
meta_testing$Season<-factor(meta_testing$Season,levels=c("Spring","Summer","Autumn","Winter"),ordered=T)
meta_testing$Sex<-factor(meta_testing$Sex)
meta_testing$Weight_Status<-factor(meta_testing$Weight_Status,levels=c("Underweight","Normal Weight","Overweight","Obese","Severe Obesity","Morbid Obesity","Super Obese"))
meta_testing$Event_Location<-factor(meta_testing$Event_Location)


data_testing<-data.frame(otu_table(ave_dta_clean))
colnames(data_testing)<-gsub("X2","2",colnames(data_testing))
colnames(data_testing)<-gsub("[.]","-",colnames(data_testing))

fit<-adonis2(t(as.matrix(data_testing))~(Estimated_PMI+Race+Manner.of.Death+Season+Sex+Weight_Status+Event_Location+BMI+Age), data=meta_testing)

########### Univariate PerMANOVA
# PMI
fit_PMI<-adonis2(t(as.matrix(data_testing))~(Estimated_PMI), data=meta_testing)
fit_PMI
# Race
fit_Race<-adonis2(t(as.matrix(data_testing))~(Race), data=meta_testing)
fit_Race
# MoD
fit_MoD<-adonis2(t(as.matrix(data_testing))~(Manner.of.Death), data=meta_testing)
fit_MoD
# Season
fit_Season<-adonis2(t(as.matrix(data_testing))~(Season), data=meta_testing)
fit_Season
# Sex
fit_Sex<-adonis2(t(as.matrix(data_testing))~(Sex), data=meta_testing)
fit_Sex
# Event_Location
fit_location<-adonis2(t(as.matrix(data_testing))~(Event_Location), data=meta_testing)
fit_location
# Age
fit_Age<-adonis2(t(as.matrix(data_testing))~(Age),data=meta_testing)
fit_Age

#### Weight data NA removed
meta_weight<-meta_testing[!is.na(meta_testing$BMI),]
data_weight<-data_testing[,rownames(meta_weight)]

# BMI
fit_BMI<-adonis2(t(as.matrix(data_weight))~(BMI), data=meta_weight)
fit_BMI
# Weight_Status
fit_weight<-adonis2(t(as.matrix(data_weight))~(Weight_Status), data=meta_weight)
fit_weight



#### Without 5NA data
fit_1_way<-adonis2(t(as.matrix(data_weight))~(Estimated_PMI+Race+Manner.of.Death+Age+BMI+Weight_Status),data=meta_weight)
fit_1_way


# Backward model selection
adonis2(t(as.matrix(data_weight))~(Estimated_PMI+Race+Manner.of.Death+Age+BMI+Weight_Status)^2,data=meta_weight)
adonis2(t(as.matrix(data_weight))~(Estimated_PMI+Race+Manner.of.Death+Age+BMI+Weight_Status+Estimated_PMI:Race+Estimated_PMI:Manner.of.Death+Estimated_PMI:Age+Estimated_PMI:BMI+Estimated_PMI:Weight_Status+Race:BMI+Age:BMI+Age:Weight_Status),data=meta_weight)									
adonis2(t(as.matrix(data_weight))~(Estimated_PMI+Race+Manner.of.Death+Age+BMI+Weight_Status+Estimated_PMI:Race+Estimated_PMI:Manner.of.Death+Estimated_PMI:Age+Estimated_PMI:BMI+Estimated_PMI:Weight_Status+Age:Weight_Status),data=meta_weight)

