attach(readRDS('../../data_prep/data_splitting/train_set.rds'))

library(indicspecies)
ind_otu = sapply(levels(train_otu$Sample_Area)
  , function(x){
    cat(x)
    dta = subset(train_otu, Sample_Area==x)
    tmp = multipatt(dta[,-c(1:2)], cluster=dta$Manner.of.Death, duleg=T)
    buccal_ind = rownames(subset(tmp$sign, p.value<=0.05))
  }
  )

ind_otu
saveRDS(ind_otu, 'ind_otu.rds')

