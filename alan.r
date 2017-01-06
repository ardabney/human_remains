#### Aggregation function of phyloseq object
phyloseq_summarize_taxa <- function(psdata, taxonomic.ranks = rank_names(psdata)) {
	require(plyr)
  if(length(taxonomic.ranks) > 1) {
    names(taxonomic.ranks) <- taxonomic.ranks
    llply(taxonomic.ranks, phyloseq_summarize_taxa, psdata = psdata)
  } else {
    taxa <- as(tax_table(psdata)[, taxonomic.ranks], 'character')
    sum_tax_table <- summarize_taxa(as(otu_table(psdata), 'matrix'), taxa)
    phyloseq(otu_table(sum_tax_table, taxa_are_rows = TRUE), sample_data(psdata, FALSE))    
  }
}

summarize_taxa <- function(counts, taxonomy) {
  if(is.matrix(taxonomy)) {
    #message('multiple taxonomies')
    alply(taxonomy, 2, summarize_taxa, counts = counts, .dims = TRUE)
  } else if(is.matrix(counts)) {
    #message('multiple counts')
    require('plyr')
    apply(counts, 2, summarize_taxa, taxonomy = taxonomy)
  } else {
    #message('summarize')
    tapply(counts, taxonomy, sum)
  }
}


####
#### Initial processing.
####

# Loading package: reference: https://joey711.github.io/phyloseq/install.html
library(phyloseq)
library(vegan)
library(glmnet)
library(DESeq2)
library(nnet)

# Load-in the rarefied OTU table
std_dta<-import_biom("../data/f120_r1k.biom")

# Load-in the meta-data
meta_dta<-read.delim2("../data/map_f120.txt",na.strings="na")
row.names(meta_dta)<-meta_dta$X.SampleID
meta_dta$BMI<-as.numeric(as.character(meta_dta$BMI))

# Reconstruct the OTU file by taking the average OTU of each body over sample_area
std_meta_dta<-data.frame(sample_data(std_dta))
std_otu_dta<-data.frame(otu_table(std_dta))
std_taxa_dta<-data.frame(tax_table(std_dta))

## The 'meta_dta' and 'std_meta_dta' data frames have differing values of BMI.
summary(meta_dta)
summary(std_meta_dta)
summary(meta_dta$BMI)
summary(as.numeric(std_meta_dta$BMI))

with(meta_dta, meta_dta[Pack_ID == "2015-S48", ])
with(std_meta_dta, std_meta_dta[Pack_ID == "2015-S48", ])

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
cn <- colnames(ave_meta_dta)
for(j in cn[cn != "BMI" & cn != "Age"])
  ave_meta_dta[, cn == j] <- factor(ave_meta_dta[, cn == j])

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
ave_taxa_dta$Species<-gsub("s__","",std_taxa_dta$Rank7)

## Combine into a phyloseq format
ave_dta<-phyloseq(
	otu_table(round(as.matrix(ave_otu_dta)),taxa_are_rows=T),
	tax_table(as.matrix(ave_taxa_dta)),
	sample_data(ave_meta_dta)
	)
# rm 0 count otu
ave_dta_clean<-filter_taxa(ave_dta,function(x) sum(x)>0,T)

# alpha diversity
alpha_temp<-data.frame(estimate_richness(ave_dta,measures=c("Simpson","Shannon")))
alpha_diversity<-data.frame(alpha_temp,Pack_ID=rownames(alpha_temp))
alpha_data<-merge(ave_meta_dta,alpha_diversity,by="Pack_ID")
alpha_data$BMI<-as.numeric(alpha_data$BMI)
alpha_data$Age<-as.numeric(alpha_data$Age)

##
## Prepare for permanova testing.
##

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

## Recode PMI by merging the earliest two time intervals. Also replace instances of '<' 
## and '>' with text.
meta_testing$Estimated_PMI[meta_testing$Estimated_PMI == "12"] <- "<24"
meta_testing$Estimated_PMI <- factor(as.character(meta_testing$Estimated_PMI))
levels(meta_testing$Estimated_PMI)<-c("Less24","More24","More48","More72")

## Drop samples with missing weight information.
meta_weight<-meta_testing[!is.na(meta_testing$BMI),]
data_weight<-data_testing[,rownames(meta_weight)]

## Drop empty factor levels for PMI and MoD.
meta_weight$Estimated_PMI <- droplevels(meta_weight$Estimated_PMI)
meta_weight$Manner.of.Death <- droplevels(meta_weight$Manner.of.Death)

##
## Multinomial regression for predicting PMI and MoD after DESeq for OTU selection.
##

## Split dataset into training and testing. Choose shrinkage factor based on training 
## data.
set.seed(101)
ID_train<-sample(1:115,round(115*0.6))
data_weight_train_temp<-data_weight[,ID_train]
# Filtering
data_weight_train<-data_weight_train_temp[rowSums(data_weight_train_temp)>=(dim(data_weight_train_temp)[2]),]
meta_weight_train<-meta_weight[ID_train,]
data_weight_test<-data_weight[,-ID_train]
meta_weight_test<-meta_weight[-ID_train,]
testing_set<-data.frame(meta_weight_test,t(data_weight_test))

##
## Manual feature selection using DESeq on training data.
##

## DESeq testing for differentially expressed OTUs with respect to PMI.
de_pmi<-DESeqDataSetFromMatrix(data_weight_train,meta_weight_train,design=~Estimated_PMI)
temp_pmi<-DESeq(de_pmi,test="LRT",reduced=~1,fitType="parametric",minReplicatesForReplace=Inf)
padj_pmi<-subset(data.frame(results(temp_pmi)),select=padj)
hist(results(temp_pmi)$pvalue)

## Extract the top 15 OTUs.
pmi_selected <- padj_pmi[order(padj_pmi), , drop = FALSE]
selected_otu_pmi<-rownames(pmi_selected)[1:15]
X_selected_pmi_otu<-t((data_weight_train[selected_otu_pmi,]))
X_selected_pmi<-data.frame(subset(meta_weight_train,select=c("Estimated_PMI","Age","Sex","Season","BMI")),X_selected_pmi_otu)

## Fit multinomial model without OTUs. Predict 'More24' every time.
fit_pmi_0 <- multinom(Estimated_PMI ~ ., data = subset(meta_weight_train, select = 
  c("Estimated_PMI","Age","Sex","Season","BMI")), maxit = 10000)
pred_pmi_0 <- predict(fit_pmi_0, newdata = testing_set)
table(pred_pmi_0, testing_set$Estimated_PMI)

## Fit multinomial model using selected features.
fit_pmi_selected_otu<-multinom(Estimated_PMI~.,data=X_selected_pmi,maxit=10000)
pred_pmi_selected_otu <- predict(fit_pmi_selected_otu,newdata=testing_set)
table(pred_pmi_selected_otu,testing_set$Estimated_PMI)


## DESeq testing for differentially expressed OTUs with respect to MoD.
de_mod<-DESeqDataSetFromMatrix(data_weight_train,meta_weight_train,design=~Manner.of.Death)
temp_mod<-DESeq(de_pmi,test="LRT",reduced=~1,fitType="parametric",minReplicatesForReplace=Inf)
padj_mod<-subset(data.frame(results(temp_mod)),select=padj)
hist(results(temp_mod)$pvalue)

## Extract the top 15 OTUs.
mod_selected <- padj_mod[order(padj_mod), , drop = FALSE]
selected_otu_mod<-rownames(mod_selected)[1:15]
X_selected_mod_otu<-t((data_weight_train[selected_otu_mod,]))
X_selected_mod<-data.frame(subset(meta_weight_train,select=c("Manner.of.Death","Age","Sex","Season","BMI")),X_selected_mod_otu)

## Fit multinomial model without OTUs.
fit_mod_0 <- multinom(Manner.of.Death ~ ., data = subset(meta_weight_train, select = 
  c("Manner.of.Death","Age","Sex","Season","BMI")), maxit = 10000)
pred_mod_0 <- predict(fit_mod_0, newdata = testing_set)
table(pred_mod_0, testing_set$Estimated_PMI)

## Fit multinomial model using selected features.
fit_mod_selected_otu<-multinom(Manner.of.Death~.,data=X_selected_mod,maxit=10000)
pred_mod_selected_otu <- predict(fit_mod_selected_otu,newdata=testing_set)
table(pred_mod_selected_otu,testing_set$Estimated_PMI)



##
## Testing for sampling_area
##
std_dta
filter_dta<-filter_taxa(std_dta,function(x)sum(x)>10,T)
# add 1 to each cell for DESeq size estimation
otu_table(filter_dta)<-otu_table(filter_dta)+1
filter_dta_deseq<-phyloseq_to_deseq2(filter_dta,~Sample_Area)
sample_area_deseq<-DESeq(filter_dta_deseq,test="LRT",reduced=~1,minReplicatesForReplace=Inf)
padj_area<-subset(data.frame(results(sample_area_deseq)),select=padj)
hist(results(sample_area_deseq)$pvalue)


##
## Attempted lasso. Errors thrown. Ignore below for now.
##

temp<-data.frame(subset(meta_weight,select=c("Race","Age","Sex","BMI","Weight_Status")),t(data_weight))
X.Mat<-model.matrix(~-1+(Race+Age+Sex+BMI+Weight_Status+Weight_Status*Age)+.,data=temp)
X.Mat_train <- X.Mat[ID_train, ]
X.Mat_test <- X.Mat[-ID_train, ]
Y.PMI_train <- meta_weight$Estimated_PMI[ID_train]
Y.PMI_test <- meta_weight$Estimated_PMI[-ID_train]
Y.MoD_train <- meta_weight$Manner.of.Death[ID_train]
Y.MoD_test <- meta_weight$Manner.of.Death[-ID_train]

# PMI
PMI_lambda<-cv.glmnet(x=X.Mat_train,y=Y.PMI_train,family="multinomial",penalty.factor=c(rep(0,11),rep(1,922),rep(0,6)),type="class")
temp_coef<-(coef(PMI_lambda,s="lambda.min"))
coef_PMI<-cbind(temp_coef$"Less24"[,1],temp_coef$"More24"[,1],temp_coef$"More48"[,1],temp_coef$"More72"[,1])

PMI_non_zero_otu <- (1:922)[apply(coef_PMI[13:934, ], 1, function(x) { any(x != 0) })]
PMI_predictors <- X.Mat[, c(1:11, PMI_non_zero_otu + 11, 934:939)]

PMI_fit <- glmnet(X.Mat, Y.PMI, family = "multinomial", penalty.factor = c(rep(0, 11), 
  rep(1, 922), rep(0, 6)), lambda = PMI_lambda$lambda.min)
PMI_hat <- predict(PMI_fit, newx = X.Mat, type = "class")

table(PMI_hat, Y.PMI)

# MoD
MoD_lambda<-cv.glmnet(x=X.Mat,y=Y.MoD,family="multinomial",penalty.factor=c(rep(0,11),rep(1,922),rep(0,6)),type="class")
