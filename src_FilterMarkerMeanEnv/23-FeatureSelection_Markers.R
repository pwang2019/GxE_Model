######################################################################################################################################
#Feature selection ----------------- RFE on markers and/or Env
#Option 1: Only RFE on markers, combine phenotype data and weather data matrix with marker matrix
#features_df_geno = cbind(phen.env.f, MAF)		#Marker_Data is the filtered 4260 markers
MAF = as.data.frame(MAF)
MAF$Variety = rownames(MAF)

features_df_geno = inner_join(modelling_data[, c(4,6,7)], MAF)			#only selection first 19 columns without average Env, only site info

MAF = MAF %>% select(-"Variety")

#############################################
#Option 2: RFE on marker and env
MAF = as.data.frame(MAF)
MAF$Variety = rownames(MAF)

features_df_geno_env = inner_join(modelling_data[, -c(1:3,5,7)], MAF)		

MAF = MAF %>% select(-"Variety")

######################################################################################################################################
##Normalise predictors for rfe and importance ranking
#Option 1: use normalisation
predictors_geno = features_df_geno[, -c(1:3)]
normalization_geno = preProcess(predictors_geno)
predictors_geno = predict(normalization_geno, predictors_geno)
predictors_geno = as.data.frame(predictors_geno)

#Option 2: no normalisation
#predictors_geno = features_df_geno[, -c(1:3)]

######################################################################################################################################
##do RFE
#subsetSizes = unique(round(1.3^(0:16)))
subsetSizes = seq(1, ncol(MAF), 100)
set.seed(123)
seeds = vector(mode = "list", length = 11)
for(i in 1:10) seeds[[i]] = sample.int(1000, length(subsetSizes) + 1)
seeds[[11]] = sample.int(1000, 1)

##################################################################
#getting ready parameters
recursive_feature_elim_control = rfeControl(functions = rfFuncs, 
                                             method = "cv", 
                                             number = 10,
                                             verbose = T,
                                             allowParallel = T)#,
                                             #seeds = seeds)

##################################################################
#register parallel
CPU_Num = 20
cl = makePSOCKcluster(CPU_Num)
registerDoParallel(cl)

##############################
#ZS49 RFE: Normalised predictors: 8 SNP, 
na.rows_geno = which(is.na(features_df_geno$ZS49.days))
rfe_results_ZS49_FMarker_AEnv_geno = rfe(predictors_geno[-na.rows_geno,], 
                        y = features_df_geno$ZS49.days[-na.rows_geno], 
                        sizes = subsetSizes,
                        rfeControl = recursive_feature_elim_control)

##############################						
#GrYld RFE: Normalised predictors: 14 SNP, 
na.rows_geno = which(is.na(features_df_geno$GrYld.kg.ha))
rfe_results_GrYld_FMarker_AEnv_geno = rfe(predictors_geno[-na.rows_geno,], 
                         y = features_df_geno$GrYld.kg.ha[-na.rows_geno],
                         sizes = c(unique(round(1.3^(0:16)))),
                         rfeControl = recursive_feature_elim_control)
						 
##############################
stopCluster(cl)
registerDoSEQ()

######################################################################################################################################
##explore Feature Selection for each tested trait
##############################
##ZS49
#RFE plot, 14 then 51, then 4358
rfe_ZS49 = ggplot(rfe_results_ZS49_geno) + 
  theme_cowplot(12)
ggsave("outputs/plots/rfe_FMarker_AEnv_ZS49_geno.png", rfe_ZS49, width = 6, height = 4)

#Importance
imp_rfe_ZS49 = ggplot(data = rfe_results_ZS49_geno$fit$importance %>% as.data.frame() %>% rownames_to_column()) +
  geom_bar(mapping = aes(y = `%IncMSE`, x = reorder(rowname, `%IncMSE`)), stat = "identity") +
  coord_flip() +
  ylab("Importance (% Increase in MSE When Omitted)") +
  xlab("Environmental Covariate") +
  theme_cowplot(12)
ggsave("outputs/plots/Imp_FMarker_AEnv_ZS49_geno.png", imp_rfe_ZS49, width = 6, height = 4)

####################
##Yld
#RFE, 8, 23, 4358
rfe_GrYld = ggplot(rfe_results_GrYld_geno) +
  theme_cowplot(12) #+
ggtitle("Recursive Feature Elimination for Environmental Variables")
ggsave("outputs/plots/rfe_FMarker_AEnv_GrYld_geno.png", rfe_GrYld, width = 6, height = 4)

#Importance
imp_GrYld = ggplot(data = rfe_results_GrYld_geno$fit$importance %>% as.data.frame() %>% rownames_to_column() %>% top_n(25, `%IncMSE`)) +
  geom_bar(mapping = aes(y = `%IncMSE`, x = reorder(rowname, `%IncMSE`)), stat = "identity") +
  coord_flip() +
  ylab("Importance (% Increase in MSE When Omitted)") +
  xlab("Environmental Covariate") +
  theme_cowplot(12)
ggsave("outputs/plots/Imp_FMarker_AEnv_GrYld_geno.png", imp_GrYld, width = 6, height = 4)

######################################################################################################################################
#now combine key features for the tested traits
#Option 1: use RFE selected markers
shared_predictors_geno = rfe_results_ZS49_geno$optVariables %in% rfe_results_GrYld_geno$optVariables
combined_predictors_geno = c(rfe_results_ZS49_geno$optVariables, rfe_results_GrYld_geno$optVariables) %>% unique()
refined_features_geno = features_df_geno %>% select(ZS49.days, GrYld.kg.ha, combined_predictors_geno)

#Option 2: use all 4260 filtered markers
refined_features_geno = features_df_geno

###################################################################################################################################### 
