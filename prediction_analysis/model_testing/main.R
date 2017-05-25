mod = readRDS('../model_training/pred.rds')

test_dta = readRDS('../data_prep/data_splitting/test_data.rds')
test_dta$BMI = as.numeric(test_dta$BMI)
test_dta$Age = as.numeric(test_dta$Age)

sapply(unique(test_dta$Sample_Area),
  function(x){
    tmp_dta = subset(test_dta, Sample_Area==x)
    tmp_mod = mod[[x]]
    predict.glm(tmp_mod, newdata=tmp_dta, type='response')
  }
  )
