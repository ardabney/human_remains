library(phyloseq)
dta = import_biom('../../../data/f120_r1k.biom')

pre_clean = subset(data.frame(sample_data(dta)), select=c('Pack_ID', 'Sample_Area'))
pre_clean$id = rownames(pre_clean)

clean_tmp = reshape2::dcast(pre_clean, Pack_ID~Sample_Area, value.var='id', drop=F)
clean_tmp$na_val = apply(is.na(clean_tmp), 1, sum)
selected_pack = clean_tmp[clean_tmp$na_val==0, 'Pack_ID']
length(selected_pack)

meta_dta = subset(data.frame(sample_data(dta)), Pack_ID%in%selected_pack, select=c('Manner.of.Death', 'Pack_ID', 'Event_Location', 'BMI', 'Race', 'Season', 'Sex', 'Weight_Status', 'Sample_Area', 'Age'))
otu_dta = subset(data.frame(otu_table(dta)), select=rownames(meta_dta))

saveRDS(
  list(meta_data=meta_dta
    , otu_data=otu_dta)
  , 'clean_data.rds')
