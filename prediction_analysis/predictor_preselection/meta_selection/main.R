attach(readRDS('../../data_prep/data_splitting/train_set.rds'))

train_meta$BMI = as.numeric(train_meta$BMI)
train_meta$Age = as.numeric(train_meta$Age)

ind_meta = sapply(unique(train_meta$Sample_Area)
  , function(x){
    cat(x)
    dta = subset(train_meta, Sample_Area==x, select=c('Manner.of.Death', 'Event_Location', 'BMI', 'Race', 'Season', 'Sex', 'Weight_Status', 'Age'))

    tmp1 = glm(Manner.of.Death~., family=binomial(logit), data=dta)
    rownames(subset(data.frame(anova(tmp1, test='Chisq')), Pr..Chi.<=0.15))
  }
  )

ind_meta
saveRDS(ind_meta, 'ind_meta.rds')

