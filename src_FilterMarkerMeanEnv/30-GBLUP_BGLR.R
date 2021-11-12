######################################################################################################################################
#function of CV split data for training, work on 1 feature, e.g. zs49
partition = function(pheno = NULL, test.size = 0.2, Do.Scale=TRUE){
  #Vector of row numbers with missing values
  na.omitter <<- which(pheno %>% is.na())
  #Exclude missing values in pheno
  Y = pheno[-na.omitter]
  #Scale pheno
  if(Do.Scale){
    y <<- (Y-mean(Y))/sd(Y)
  }
  #no scale
  if(!Do.Scale){
    y <<- Y
  }
  #Create testing set
  yNA <<- y
  test.index <<- sample(1:length(y), round(0.2*length(y)))
  yNA[test.index] <<- NA
  #  print(list("training y" = yNA,
  #             "na.omitter" = na.omitter))
}
######################################################################################################################################
######################################################################################################################################
#function of CV split data for training, work on 1 feature, e.g. zs49, leave one Site as testing
partition_site = function(pheno = NULL, model_data = NULL, test.site = NULL, Do.Scale=TRUE){
  #Vector of row numbers with missing values
  na.omitter <<- which(pheno %>% is.na())
  #Exclude missing values in pheno
  Y = pheno[-na.omitter]
  tmat = model_data[-na.omitter, ]
  
  #Scale pheno
  if(Do.Scale){
    y <<- (Y-mean(Y))/sd(Y)
  }
  #no scale
  if(!Do.Scale){
    y <<- Y
  }
  
  #Create testing set
  yNA <<- y
  test.index <<- which(tmat$Location == test.site)
  yNA[test.index] <<- NA
}
######################################################################################################################################
#now get training and testing ZS49, filter NA ones
#NOTE: DO one trait at a time!
set.seed(123)
#Option 1: use 80-20 CV
partition(refined_features$ZS49.days, 0.2)
partition(refined_features$GrYld.kg.ha, 0.2)

#################################################
#Option 2: reserve a site for testing, model will not see that site for training.
#5 sites: GER ESP SoP KAT MER
Test_Site = "KAT"
partition_site(refined_features$ZS49.days, modelling_data, Test_Site)
partition_site(refined_features$GrYld.kg.ha, modelling_data, Test_Site)

######################################################################################################################################
##Create Genetic Relationship Matrix, and then get eigenvalues
GRM = snpReady::G.matrix(MAF, method = "VanRaden", format = "wide")$Ga
Zb = predict(dummyVars(~ Variety, data = modelling_data[-na.omitter,]),modelling_data[-na.omitter,]) %>% as.matrix()
ZbGZbp = Zb %*% GRM %*% t(Zb)

EVD = eigen(ZbGZbp)

WRM = tcrossprod(as.matrix(refined_features[-na.omitter,-(1:5)]))
EV.WRM = eigen(WRM)

GWHa = hadamard.prod(ZbGZbp, WRM)
EV.GWHa = eigen(GWHa)

#############################################################
#y = G			no E
ETA_GBLUP = list(
  list(V = EVD$vectors, d = EVD$values, model = 'RKHS'),
  list(~ as.factor(Trial), 
       data = cbind(modelling_data[,1:3], refined_features[,-(1:5)])[-na.omitter,], 
       model = "FIXED")
)
#############################################################
#y = G + E is stratified
ETA_GBLUP_G = list(
  list(V = EVD$vectors, d = EVD$values, model = 'RKHS'),
  list(~ as.factor(Year):as.factor(Location):as.factor(Trial), 
       data = cbind(modelling_data[,1:3], refined_features[,-(1:5)])[-na.omitter,], 
       model = "FIXED")
)
#############################################################
#y = G x E
ETA_GBLUP_GEGxE = list(
  list(V = EVD$vectors, d = EVD$values, model = 'RKHS'),
  list(V = EV.WRM$vectors, d = EV.WRM$values, model = 'RKHS'),
  list(V = EV.GWHa$vectors, d = EV.GWHa$values, model = 'RKHS')
)

######################################################################################################################################
#Here be models::::::::::::::::::::::::::::::
ITER = 1200
BURN = 200
VERB = FALSE

############################################
#register parallel
CPU_Num = 20
cl = makePSOCKcluster(CPU_Num)
registerDoParallel(cl)

#############################################################
#y = G
GBLUP = BGLR(y = yNA, 
                ETA = ETA_GBLUP, 
                nIter = ITER, 
                burnIn = BURN,
                saveAt = "outputs/models/GBLUP",
                verbose = VERB)

#############################################################
#y = G + stratified E
GBLUP_G = BGLR(y = yNA, 
                ETA = ETA_GBLUP_G, 
                nIter = ITER, 
                burnIn = BURN,
                saveAt = "outputs/models/GBLUP_G",
                verbose = VERB)

#############################################################
#y = G + E + GxE
GBLUP_GEGxE = BGLR(y = yNA, 
                  ETA = ETA_GBLUP_GEGxE, 
                  nIter = ITER, 
                  burnIn = BURN,
                  saveAt = "outputs/models/GBLUP_GEGxE",
                  verbose = VERB)

############################################																			 
#stop parallel process
stopCluster(cl)
registerDoSEQ()

######################################################################################################################################
#compare results
cor(GBLUP$yHat[test.index], y[test.index])
postResample(GBLUP_G$yHat[test.index], y[test.index])
 
######################################################################################################################################
#compare results
cor(GBLUP_G$yHat[test.index], y[test.index])
postResample(GBLUP_G$yHat[test.index], y[test.index])

#############################################################
cor(GBLUP_GEGxE$yHat[test.index], y[test.index])
postResample(GBLUP_GEGxE$yHat[test.index], y[test.index])

######################################################################################################################################
##Save models and y var since it will be changed every time of training
GBLUP_G_FMarker_AEnv_ZS49 = GBLUP
GBLUP_GE_FMarker_AEnv_ZS49 = GBLUP_G
GBLUP_GEGxE_FMarker_AEnv_ZS49 = GBLUP_GEGxE
y_FMarker_AEnv_ZS49 = y
test_index.FMarker_AEnv_ZS49 = test.index

GBLUP_G_FMarker_AEnv_GrYld = GBLUP
GBLUP_GE_FMarker_AEnv_GrYld = GBLUP_G
GBLUP_GEGxE_FMarker_AEnv_GrYld = GBLUP_GEGxE
y_FMarker_AEnv_GrYld = y
test_index.FMarker_AEnv_GrYld = test.index

GBLUP_G_FMarker_AEnv_ZS49_NewSite = GBLUP
GBLUP_GE_FMarker_AEnv_ZS49_NewSite = GBLUP_G
GBLUP_GEGxE_FMarker_AEnv_ZS49_NewSite = GBLUP_GEGxE
y_FMarker_AEnv_ZS49_NewSite = y
test_index.FMarker_AEnv_ZS49_NewSite = test.index

GBLUP_G_FMarker_AEnv_GrYld_NewSite = GBLUP
GBLUP_GE_FMarker_AEnv_GrYld_NewSite = GBLUP_G
GBLUP_GEGxE_FMarker_AEnv_GrYld_NewSite = GBLUP_GEGxE
y_FMarker_AEnv_GrYld_NewSite = y
test_index.FMarker_AEnv_GrYld_NewSite = test.index

######################################################################################################################################
##save models
save(GBLUP_G_FMarker_AEnv_ZS49, GBLUP_GE_FMarker_AEnv_ZS49, GBLUP_GEGxE_FMarker_AEnv_ZS49, y_FMarker_AEnv_ZS49,
     GBLUP_G_FMarker_AEnv_GrYld,GBLUP_GE_FMarker_AEnv_GrYld,GBLUP_GEGxE_FMarker_AEnv_GrYld,y_FMarker_AEnv_GrYld,
     file = "outputs/models/GBLUP_FMarker_AEnv_80-20.RData")

save(GBLUP_G_FMarker_AEnv_ZS49_NewSite, GBLUP_GE_FMarker_AEnv_ZS49_NewSite, GBLUP_GEGxE_FMarker_AEnv_ZS49_NewSite, y_FMarker_AEnv_ZS49_NewSite,
     GBLUP_G_FMarker_AEnv_GrYld_NewSite,GBLUP_GE_FMarker_AEnv_GrYld_NewSite,GBLUP_GEGxE_FMarker_AEnv_GrYld_NewSite,y_FMarker_AEnv_GrYld_NewSite,
     file = "outputs/models/GBLUP_FMarker_AEnv_site.RData")

#############################################################
#Can extract test string from which(model$y %>% is.na)

######################################################################################################################################
