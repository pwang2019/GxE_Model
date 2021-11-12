######################################################################################################################################
#this is the good one, 30543 SNPs, in ACGT coding --------
Marker_Gaofeng_filter = read_excel("compiled_data/Genotype-homozygous_wo_9chromosome.hmp.xlsx", trim_ws = TRUE)
Marker_Gaofeng_filter = as.data.frame(Marker_Gaofeng_filter)

############################################
#register parallel
CPU_Num = 20
cl = makePSOCKcluster(CPU_Num)
registerDoParallel(cl)

ref_al = substring(Marker_Gaofeng_filter$alleles, 1, 1)
alt_al = substring(Marker_Gaofeng_filter$alleles, 3, 3)

tgeno_data = Marker_Gaofeng_filter
for(m in 1:nrow(Marker_Gaofeng_filter)){
    tgeno_data[m, tgeno_data[m, ]==ref_al[m]] = 0
    tgeno_data[m, tgeno_data[m, ]==alt_al[m]] = 2
	tgeno_data[m, tgeno_data[m, ]=="N"] = NA
}
for(m in 1:nrow(Marker_Gaofeng_filter)){
    tgeno_data[m, which(tgeno_data[m, ]==ref_al[m])] = 0
    tgeno_data[m, which(tgeno_data[m, ]==alt_al[m])] = 2
	tgeno_data[m, which(tgeno_data[m, ]=="N")] = NA
}
Marker_Gaofeng_filter = tgeno_data

############################################
#Save this time consuming process 
save(Marker_Gaofeng_filter, file="outputs/Marker_all_processed.RData")

############################################																			 
#stop parallel process
stopCluster(cl)
registerDoSEQ()

############################################
#clean up
rm(ref_al, alt_al)

######################################################################################################################################
# Filter for varieties that have both genotype and phenotype data--------------
#marker_varieties = colnames(Marker_Data[,8:ncol(Marker_Data)])
marker_varieties = colnames(Marker_Gaofeng_filter[,14:ncol(Marker_Gaofeng_filter)])

phenotype_varieties = phenology$Variety %>% unique() %>% as.vector()

#Phen vars in marker vars
union.vars = phenotype_varieties[phenotype_varieties %in% marker_varieties]

all_varieties = unique(c(marker_varieties, phenotype_varieties))
removed_varieties = all_varieties[which(!(all_varieties %in% union.vars))]

#Marker_Data = Marker_Data %>% select_if(!(colnames(Marker_Data) %in% removed_varieties))
Marker_Gaofeng_filter = Marker_Gaofeng_filter %>% select_if(!(colnames(Marker_Gaofeng_filter) %in% removed_varieties))

for(m in 14:ncol(Marker_Gaofeng_filter)){
    Marker_Gaofeng_filter[, m] = as.numeric(Marker_Gaofeng_filter[, m])
}
######################################################################################################################################
#Exclude phenology varieties with no Marker data
phenology = phenology %>% filter(!(phenology$Variety %in% removed_varieties))

#########################################################
rm(marker_varieties, phenotype_varieties, union.vars, all_varieties, removed_varieties)

######################################################################################################################################
#Transpose, set colnames to SNP
#MAF = Marker_Data[,-c(1:7)] %>% t()
#colnames(MAF) = Marker_Data$ID
MAF = Marker_Gaofeng_filter[, -c(1:13)] %>% t() 
colnames(MAF) = Marker_Gaofeng_filter$rs#

#Imputing missing genotypes
MAF = raw.data(MAF, 
                    frame = "wide", # Markers as cols, samples as rows
                    base = FALSE, # Already coded as 0, 1, 2
                    sweep.sample = 1,  
                    maf = 0.001, #Frequency of second most common allele
                    call.rate = 0,
                    imput = TRUE,
                    imput.type = "wright",
                    outfile = "012"
                )$M.clean
######################################################################################################################################