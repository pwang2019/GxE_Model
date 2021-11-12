######################################################################################################################################
#Create dataframe with all variables to be used
ML_dataframe = modelling_data %>% select(-c(Year, Trial, Location, Latitude, Longitude, soiltype)) %>% 
  cbind(markers)
prep = preProcess(ML_dataframe)
ML_dataframe = predict(prep, ML_dataframe)

#############################################################
zs49.na = which(ML_dataframe$ZS49.days %>% is.na())
ML_df_ZS49 = ML_dataframe[-zs49.na,] %>% 
  select(-GrYld.kg.ha)

gryld.na = which(ML_dataframe$GrYld.kg.ha %>% is.na())
ML_df_GrYld = ML_dataframe[-gryld.na,] %>% 
  select(-ZS49.days)

#############################################################
mtry.EC = round(ncol(ML_dataframe)/3)
rf.tunegrid.EC = expand.grid(.mtry = mtry.EC)

######################################################################################################################################
#explore ZS49 and GrYld
boxplot(modelling_data$ZS49.days)
hist(modelling_data$ZS49.days, breaks=50)

boxplot(modelling_data$GrYld.kg.ha)
hist(modelling_data$GrYld.kg.ha, breaks=50)

######################################################################################################################################
# set up RF control parameteres
rf.control = trainControl(method = 'cv',
                           number = 5,
                           verboseIter = T,
                           allowParallel = T,
                           savePredictions = "final")

######################################################################################################################################
#now predict G from E and P
library(randomForestSRC)
library(xlsx)

############################################################
#save some key variables for ML prediction, MEM!
save(rfe_results_ZS49_geno, rfe_results_GrYld_geno, modelling_data, ML_dataframe, markers, file ="Key_DataFrame.RData")

#get all Markers' name
tnames = colnames(ML_dataframe)[95:ncol(ML_dataframe)]

######################################################################################################################################
#Build model on top ranking cluster of markers, zs 10 cluster, yld 13 cluster
#get the best cluster SNP lists 
top_snp_zs49 = read.xlsx2("Genotypes_KeyHaplotype_zs49.xlsx", sheetIndex=1)		#87 SNPs
rownames(top_snp_zs49) = top_snp_zs49$ID

top_snp_yld = read.xlsx2("Genotypes_KeyHaplotype_yld.xlsx", sheetIndex=1)		#60 SNPs
rownames(top_snp_yld) = top_snp_yld$ID

############################################
#get training and test data
G.trn = sample(nrow(ML_dataframe), 0.8*nrow(ML_dataframe))
#G.trn = sample(nrow(ML_df_GrYld), 0.8*nrow(ML_df_GrYld))

########################################################################################
#do yld
#register parallel
CPU_Num = 20
cl = makePSOCKcluster(CPU_Num)
registerDoParallel(cl)

tvar = paste(top_snp_yld$ID, collapse=",")
tform = as.formula(paste("c(",tvar, ")~.", sep=""))

#G.yld.var60 = rfsrc(cbind(top_snp_yld$ID)~., data = ML_dataframe[G.trn, ], na.action = "na.omit")
#G.yld.var60 = rfsrc( top_snp_yld$ID ~ ., data = ML_dataframe[G.trn, ], na.action = "na.impute")
G.yld.var60 = rfsrc(cbind(snp_1_359703846, snp_1_359703848, snp_1_359703985, snp_1_421515942,
    snp_1_421516158, snp_1_421516318, snp_1_421516359, snp_2_80405222,
    snp_2_80405291, snp_2_80405293, snp_2_160353179, snp_2_160354963,
    snp_2_160354964, snp_2_160355023, snp_2_160355086, snp_2_160355232,
    snp_3_112633626, snp_3_112634446, snp_3_112634929, snp_3_634078704,
    snp_3_634079937, snp_3_634079978, snp_3_634189000, snp_3_634189663,
    snp_3_634190153, snp_3_634928088, snp_3_634931110, snp_4_627080516,
    snp_4_627080614, snp_4_627082280, snp_4_627082555, snp_5_131870539,
    snp_5_131998865, snp_5_132001788, snp_5_243851853, snp_5_243996248,
    snp_5_244021572, snp_5_244024104, snp_5_559673235, snp_5_560195991,
    snp_5_560570138, snp_5_560570331, snp_5_560587395, snp_5_560588247,
    snp_5_560588251, snp_5_560588278, snp_5_560588280, snp_5_560588323,
    snp_5_598231913, snp_5_598232043, snp_5_599122947, snp_5_599329482,
    snp_5_599329864, snp_5_599331920, snp_7_35129549, snp_7_35488946,
    snp_7_631947468, snp_7_631947469, snp_7_631948152, snp_7_631949363) ~., data = ML_dataframe[G.trn, ], na.action = "na.impute")
	
stopCluster(cl)
registerDoSEQ()

plot.variable(G.yld.var)

G.yld.pred = predict(G.yld.var60, newdata = ML_dataframe[-G.trn, ], na.action = "na.impute")
G.yld.pred = predict(G.yld.var60, newdata = ML_dataframe[-G.trn, ] %>% select(-top_snp_yld$ID), na.action = "na.impute")
G.yld.y = ML_dataframe[-G.trn, top_snp_yld$ID]
sum(G.yld.pred$yvar == G.yld.y)

########################################################################################
#do zs49

############################################
#get training and test data
G.trn = sample(nrow(ML_dataframe), 0.8*nrow(ML_dataframe))
#G.trn = sample(nrow(ML_df_ZS49), 0.8*nrow(ML_df_ZS49))

tvar = paste(top_snp_zs49$ID, collapse=",")
tform = as.formula(paste("c(",tvar, ")~.", sep=""))

#register parallel
CPU_Num = 20
cl = makePSOCKcluster(CPU_Num)
registerDoParallel(cl)

tvar = paste(top_snp_zs49$ID, collapse=",")
tform = as.formula(paste("c(",tvar, ")~.", sep=""))

G.zs49.var87 = rfsrc(cbind(snp_1_556899940, snp_1_556899955, snp_1_556902615, snp_1_556902668,
    snp_1_558226167, snp_1_558226415, snp_1_558227579, snp_1_558229862,
    snp_2_523377147, snp_2_523378047, snp_2_523378077, snp_2_523378515,
    snp_2_523379223, snp_2_523379290, snp_2_523379293, snp_2_523379296,
    snp_2_523379368, snp_2_592334612, snp_2_592335927, snp_2_592336146,
    snp_3_117874739, snp_3_117876749, snp_3_119252382, snp_3_119253152,
    snp_4_645064549, snp_4_645064767, snp_4_645065364, snp_4_645065411,
    snp_4_645065921, snp_5_559673560, snp_5_559687961, snp_5_560569994,
    snp_5_560570097, snp_5_560570331, snp_5_560571241, snp_5_560571491,
    snp_5_560587293, snp_5_598229924, snp_5_598231895, snp_5_598560301,
    snp_5_599019952, snp_5_599069784, snp_5_599112501, snp_5_599122947,
    snp_5_599122952, snp_5_599329482, snp_5_599331920, snp_5_599332514,
    snp_5_599332925, snp_5_599333006, snp_5_603612478, snp_5_603612743,
    snp_5_603613102, snp_5_603613468, snp_5_603613647, snp_5_603613740,
    snp_5_603613940, snp_5_603614147, snp_5_603614684, snp_5_603614823,
    snp_5_603614892, snp_5_603615056, snp_5_603615493, snp_5_603615644,
    snp_5_603615648, snp_5_603616317, snp_5_603616894, snp_5_603616898,
    snp_5_603617202, snp_5_603617654, snp_5_603618083, snp_6_133137227,
    snp_6_133170553, snp_6_133172173, snp_6_133173262, snp_6_133174007,
    snp_7_37608111, snp_7_37609741, snp_7_37610330, snp_7_37611274,
    snp_7_37611446, snp_7_37888799, snp_7_37891608, snp_7_37903679,
    snp_7_37903766, snp_7_37904484, snp_7_37904889) ~ ., data = ML_dataframe[G.trn, ], na.action = "na.impute")
	
stopCluster(cl)
registerDoSEQ()

plot.variable(G.zs49.var87)

G.zs49.pred = predict(G.zs49.var87, newdata = ML_dataframe[-G.trn, ], na.action = "na.impute")
G.zs49.y = ML_dataframe[-G.trn, top_snp_zs49$ID]
sum(G.zs49.pred$yvar == G.zs49.y)

######################################################################################################################################
######################################################################################################################################
G.zs.var = rfsrc( cbind(snp_5_599331920, snp_1_558227579, snp_6_267736262, snp_5_599333006, snp_5_599332925, snp_6_156668468, snp_2_341083604, snp_2_432061917, snp_2_432062544, snp_2_341083143) ~ ., data = ML_dataframe[G.trn, ], na.action = "na.impute")

G.yld.var = rfsrc( cbind(snp_3_634079937, snp_3_117876749, snp_1_359703846, snp_1_359703848, snp_3_634078704) ~ ., data = ML_dataframe[G.trn, ], na.action = "na.impute")

######################################################################################################################################
######################################################################################################################################
#plot G prediction accuracy, 60 SNP yld, 87 SNP zs49
############################################################
#zs49
acc.zs49 = G.zs49.y[1, ]
for(m in 1:ncol(G.zs49.y)){
    tpred = round(G.zs49.pred$regrOutput[[m]][[1]])
	acc.zs49[m] = sum(tpred == G.zs49.y[, m])/length(tpred)*100
}
acc.zs49 = t(acc.zs49)
acc.zs49 =as.data.frame(acc.zs49)
colnames(acc.zs49) = "Accuracy"
acc.zs49$ID = rownames(acc.zs49)

acc.zs49$Chr = as.numeric(substring(acc.zs49$ID, 5, 5))
acc.zs49$Pos = as.numeric(substring(acc.zs49$ID, 7, nchar(acc.zs49$ID)))
acc.zs49$Cluster = as.factor(top_snp_zs49[acc.zs49$ID, "Cluster"])

############################################################
#yld
acc.yld = G.yld.y[1, ]
for(m in 1:ncol(G.yld.y)){
    tpred = round(G.yld.pred$regrOutput[[m]][[1]])
	acc.yld[m] = sum(tpred == G.yld.y[, m])/length(tpred)*100
}
acc.yld = t(acc.yld)
acc.yld =as.data.frame(acc.yld)
colnames(acc.yld) = "Accuracy"
acc.yld$ID = rownames(acc.yld)

acc.yld$Chr = as.numeric(substring(acc.yld$ID, 5, 5))
acc.yld$Pos = as.numeric(substring(acc.yld$ID, 7, nchar(acc.yld$ID)))
acc.yld$Cluster = as.factor(top_snp_yld[acc.yld$ID, "Cluster"])
levels(acc.yld$Cluster) = as.character(1:13)

######################################################################################################################################
######################################################################################################################################
