######################################################################################################################################
##Remove non-used pheno, ZS49 is 1st to try, leave one at a time
#modelling_data = phen.env %>% 
#  select(-c(trial_index, Reason, DOP, BOM.Stn, ZS49PlHt, ZS91.days, PPd.days,
#            PrdGrYld.kg.ha, PrdGrYld.t.ha, SEGrYld, TlrNo, PlntNo)) #%>% 
#  filter(!is.na(ZS49.days))

#alternatively, leave the phenotypes to be predicted in modelling_data
modelling_data = phen.env %>% 
  select(-c(trial_index, Reason, DOP, BOM.Stn, PPd.days,
            PrdGrYld.kg.ha, PrdGrYld.t.ha, SEGrYld, TlrNo, PlntNo)) #%>% 
  
##################################################################
#check variables, filter zero variance features
near0variance = nearZeroVar(modelling_data, saveMetrics = T)
modelling_data = modelling_data[, -which(near0variance$nzv)]

######################################################################################################################################
##join marker data with filered pheno data
markers = modelling_data %>% select(Variety) %>% 
  left_join(setDT(as.data.frame(MAF), keep.rownames = "Variety")) %>% 
  select(-Variety)

modelling_data = modelling_data %>%
  unclass() %>% 
  data.frame(stringsAsFactors = T)

######################################################################################################################################
##Dummy var set up
#dummies = dummyVars(ZS49.days ~ ., data = modelling_data %>% select(-Variety))
#modelling_data = predict(dummies, newdata = modelling_data) 