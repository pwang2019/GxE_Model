######################################################################################################################################
#Feature selection -----------------
features_df = modelling_data %>%
  select(-c(Year, Location, Trial, soiltype, Longitude, Latitude))

##################################################################
##Normalise predictors for rfe and importance ranking
#Option 1: normalised predictors to get better results? 
predictors = features_df[, -c(1:5)]
normalization = preProcess(predictors)
predictors = predict(normalization, predictors)
predictors = as.data.frame(predictors)

##################################################################
#Option 2: no normalisation on the predictors
#predictors = features_df[, -c(1:5)]

######################################################################################################################################
##Do RFE now ------------------- RF based RFE, all soil infos, Averaged Env factors
subsetSizes = unique(round(1.3^(0:16)))
set.seed(123)
seeds = vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]] = sample.int(1000, length(subsetSizes) + 1)
seeds[[11]] = sample.int(1000, 1)

#getting ready parameters
recursive_feature_elim_control = rfeControl(functions = rfFuncs, 
                                             method = "cv", 
                                             number = 10,
                                             verbose = T,
                                             allowParallel = T,
                                             seeds = seeds)

#register parallel
cl = makePSOCKcluster(20)
registerDoParallel(cl)

############################################
#ZS49 RFE --------------- normalised Predictors, 67 top var; un-normalised predictors,  top vars
na.rows = which(is.na(features_df$ZS49.days))
rfe_results_ZS49_AMarker_BOM = rfe(predictors[-na.rows,], 
                        y = features_df$ZS49.days[-na.rows], 
                        sizes = subsetSizes,
                        rfeControl = recursive_feature_elim_control)
						
############################################						
#GrYld RFE --------------- normalised Predictors, 6 top var; un-normalised predictors, top vars
na.rows = which(is.na(features_df$GrYld.kg.ha))
rfe_results_GrYld_AMarker_BOM = rfe(predictors[-na.rows,], 
                         y = features_df$GrYld.kg.ha[-na.rows],
                         sizes = c(unique(round(1.3^(0:16)))),
                         rfeControl = recursive_feature_elim_control)
																		 
stopCluster(cl)
registerDoSEQ()

######################################################################################################################################
###Explore Feature Selection for each tested trait
####################
##ZS49
#RFE plot
rfe_AMarker_BOM_ZS49 = ggplot(rfe_results_ZS49_AMarker_BOM) + theme_cowplot(12)
ggsave("outputs/plots/rfe_AMarker_BOM_ZS49.pdf", rfe_env_ZS49, width = 6, height = 4)

#Importance
imp_env_ZS49 = ggplot(data = rfe_results_ZS49_AMarker_BOM$fit$importance %>% as.data.frame() %>% rownames_to_column()) +
  geom_bar(mapping = aes(y = `%IncMSE`, x = reorder(rowname, `%IncMSE`)), stat = "identity") +
  coord_flip() + ylab("Importance (% Increase in MSE When Omitted)") + xlab("Environmental Covariate") + theme_cowplot(12)
ggsave("outputs/plots/Imp_AMarker_BOM_ZS49.pdf", imp_env_ZS49, width = 6, height = 4)

####################
##Yld
#RFE plot
rfe_env_GrYld = ggplot(rfe_results_GrYld_AMarker_BOM) + theme_cowplot(12) #+
ggtitle("Recursive Feature Elimination for Environmental Variables")
ggsave("outputs/plots/rfe_AMarker_BOM_GrYld.pdf", rfe_env_GrYld, width = 6, height = 4)

#Importance
imp_env_GrYld = ggplot(data = rfe_results_GrYld_AMarker_BOM$fit$importance %>% as.data.frame() %>% rownames_to_column() %>% top_n(25, `%IncMSE`)) +
  geom_bar(mapping = aes(y = `%IncMSE`, x = reorder(rowname, `%IncMSE`)), stat = "identity") +
  coord_flip() + ylab("Importance (% Increase in MSE When Omitted)") + xlab("Environmental Covariate") + theme_cowplot(12)
ggsave("outputs/plots/Imp_AMarker_BOM_GrYld.pdf", imp_env_GrYld, width = 6, height = 4)

######################################################################################################################################
#now combine key features for the tested traits
shared_predictors = rfe_results_ZS49_AMarker_BOM$optVariables %in% rfe_results_GrYld_AMarker_BOM$optVariables
combined_predictors = c(rfe_results_ZS49_AMarker_BOM$optVariables, rfe_results_GrYld_AMarker_BOM$optVariables) %>% unique()

#################################
#Option 1: use only RFE picked top vars for the ZS49 and GrYld phenotypes
refined_features = features_df %>% select(ZS49.days, GrYld.kg.ha, combined_predictors)

#################################
#Option 2: use all the features, ignore RFE results
refined_features = features_df

######################################################################################################################################
