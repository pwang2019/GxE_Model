######################################################################################################################################
##Create dataframe with all variables to be used --------------- use all env variables of averaged staged ones, and all markers
ML_dataframe = modelling_data %>% select(-c(Year, Trial, Latitude, Longitude, soiltype, ZS49PlHt, ZS91.days)) %>% 
  cbind(markers)
prep = preProcess(ML_dataframe)
ML_dataframe = predict(prep, ML_dataframe)

#############################################################
#Create a data frame for training ZS49 and GrYld
zs49.na = which(ML_dataframe$ZS49.days %>% is.na())
ML_df_ZS49 = ML_dataframe[-zs49.na,] %>% select(-GrYld.kg.ha)

gryld.na = which(ML_dataframe$GrYld.kg.ha %>% is.na())
ML_df_GrYld = ML_dataframe[-gryld.na,] %>% select(-ZS49.days)

#############################################################
#Initiate RF tuning param
mtry.EC = round(ncol(ML_dataframe)/3)
rf.tunegrid.EC = expand.grid(.mtry = mtry.EC)

######################################################################################################################################
##Do RF, CV1 -------------------------------- 
rf.control = trainControl(method = 'cv',
                           number = 5,
                           verboseIter = T,
                           allowParallel = T,
                           savePredictions = "final")

#############################################################
#Option 1: 80-20 all sites training and testing
set.seed(123)
ZS49.trn = sample(nrow(ML_df_ZS49), 0.8*nrow(ML_df_ZS49))
GrYld.trn = sample(nrow(ML_df_GrYld), 0.8*nrow(ML_df_GrYld))

#############################################################
#Option 2: using a single site as testing
ZS49.trn = which(ML_df_ZS49$Location == "KAT")
GrYld.trn = which(ML_df_GrYld$Location == "KAT")

#############################################################
#Do RF now. use multiple thread
CPU_Num = 20
cl = makePSOCKcluster(CPU_Num)
registerDoParallel(cl)

rf.ZS49.AMarker.AEnv = train(ZS49.days ~ . - Variety,
                 data = ML_df_ZS49[ZS49.trn,],
                 method = "rf",
                 metric = "RMSE",
                 tuneGrid = rf.tunegrid.EC,
                 trControl = rf.control
)

rf.GrYld.AMarker.AEnv = train(GrYld.kg.ha ~ . - Variety,
                      data = ML_df_GrYld[GrYld.trn,],
                      method = "rf",
                      metric = "RMSE",
                      tuneGrid = rf.tunegrid.EC,
                      trControl = rf.control
)

stopCluster(cl)
registerDoSEQ()

######################################################################################################################################
##ZS49 RF model ---------------------
ZS49.rf.pred.AMarker.AEnv = predict(rf.ZS49.AMarker.AEnv, ML_df_ZS49[-ZS49.trn,])
cor(ZS49.rf.pred.AMarker.AEnv, ML_df_ZS49$ZS49.days[-ZS49.trn])
postResample(ZS49.rf.pred.AMarker.AEnv, ML_df_ZS49[-ZS49.trn, ]$ZS49.days)

#######################################################
##GrYld RF model---------------------
GrYld.rf.pred.AMarker.AEnv = predict(rf.GrYld.AMarker.AEnv, ML_df_GrYld[-GrYld.trn,])
cor(GrYld.rf.pred.AMarker.AEnv, ML_df_GrYld$GrYld.kg.ha[-GrYld.trn])
postResample(GrYld.rf.pred.AMarker.AEnv, ML_df_GrYld[-GrYld.trn, ]$GrYld.kg.ha)


######################################################################################################################################
#Save models
ZS49_trn.AMarker.AEnv = ZS49.trn
GrYld.trn.AMarker.AEnv = GrYld.trn
rf.ZS49.AMarker.AEnv.CV = rf.ZS49.AMarker.AEnv
rf.GrYld.AMarker.AEnv.CV = rf.GrYld.AMarker.AEnv

ZS49_trn.AMarker.AEnv.NewSite = ZS49.trn
GrYld.trn.AMarker.AEnv.NewSite = GrYld.trn
rf.ZS49.AMarker.AEnv.NewSite = rf.ZS49.AMarker.AEnv
rf.GrYld.AMarker.AEnv.NewSite = rf.GrYld.AMarker.AEnv

save(rf.ZS49.AMarker.AEnv.CV, rf.GrYld.AMarker.AEnv.CV, ML_df_GrYld, ML_df_ZS49, ZS49_trn.AMarker.AEnv, GrYld.trn.AMarker.AEnv, file = "outputs/models/RF_FMarker_AEnv_80-20.RData")
save(rf.ZS49.AMarker.AEnv.NewSite, rf.GrYld.AMarker.AEnv.NewSite, ML_df_GrYld, ML_df_ZS49, ZS49_trn.AMarker.AEnv.NewSite, GrYld.trn.AMarker.AEnv.NewSite, file = "outputs/models/RF_FMarker_AEnv_site.RData")

######################################################################################################################################
