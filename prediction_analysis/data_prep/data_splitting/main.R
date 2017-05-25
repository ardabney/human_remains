attach(readRDS('../data_cleaning/clean_data.rds'))

meta_data$Manner.of.Death = factor(meta_data$Manner.of.Death)
levels(meta_data$Manner.of.Death) = c('non_natural', 'non_natural', 'natural', 'non_natural')
set.seed(101)
train_id = with(meta_data, sample(unique(Pack_ID), 0.6*length(unique(Pack_ID))))

train_meta = subset(meta_data, Pack_ID%in%train_id)
train_otu = cbind(
  Manner.of.Death=train_meta$Manner.of.Death
  , Sample_Area=train_meta$Sample_Area
  , data.frame(t(subset(otu_data, select=rownames(train_meta)))))

table(train_meta$Manner.of.Death)
saveRDS(
  list(train_meta=train_meta
    , train_otu=train_otu)
  , 'train_set.rds')

test_meta = subset(meta_data, !Pack_ID%in%train_id)
test_otu = data.frame(t(subset(otu_data, select=rownames(test_meta))))
test_data = cbind(test_meta, test_otu)

table(test_data$Manner.of.Death)
saveRDS(test_data, 'test_data.rds')
