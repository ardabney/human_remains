attach(readRDS('../data_prep/data_splitting/train_set.rds'))
train_meta$Age = as.numeric(train_meta$Age)
train_meta$BMI = as.numeric(train_meta$BMI)

meta_pred = readRDS('../predictor_preselection/meta_selection/ind_meta.rds')
otu_pred = readRDS('../predictor_preselection/otu_selection/ind_otu.rds')

pred = sapply(unique(train_meta$Sample_Area),
  function(x){
    meta_tmp = subset(train_meta, Sample_Area==x, select=c(meta_pred[,x],'Manner.of.Death'))
    otu_tmp = subset(train_otu, select=otu_pred[[x]])[rownames(meta_tmp),]
    tmp_dta = cbind(meta_tmp, otu_tmp)
    fit = glm(Manner.of.Death~., family=binomial(logit), data=tmp_dta)
    return(fit)
  }, simplify=F
  )

saveRDS(pred, 'pred.rds')    

