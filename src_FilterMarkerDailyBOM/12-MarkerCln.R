######################################################################################################################################
##Load in main Marker data --------
Marker_Data = read_excel("compiled_data/Marker_Data.xlsx", trim_ws = TRUE)
Marker_Data[8:ncol(Marker_Data)] = lapply(Marker_Data[8:ncol(Marker_Data)], as.numeric)

Marker_Data = Marker_Data %>% rename("chrom" = '#CHROM')

######################################################################################################################################
##Filter for varieties that have both genotype and phenotype data--------------
#create Vectors of varieties for intercepting
marker_varieties = colnames(Marker_Data[,8:ncol(Marker_Data)])
phenotype_varieties = phenology$Variety %>% unique() %>% as.vector()

#Phen vars in marker vars
union.vars = phenotype_varieties[phenotype_varieties %in% marker_varieties]

#Trying to bring back misspelled varieties
#marker_only = marker_varieties[which(!marker_varieties %in% union.vars)]
#phenology_only = phenotype_varieties[which(!phenotype_varieties %in% union.vars)]
#matches = sapply(marker_only, FUN = function(x) phenology_only[agrep(x, phenology_only, max.distance = 0.01)]) %>% unlist()

all_varieties = unique(c(marker_varieties, phenotype_varieties))

removed_varieties = all_varieties[which(!(all_varieties %in% union.vars))]

Marker_Data = Marker_Data %>% select_if(!(colnames(Marker_Data) %in% removed_varieties))

######################################################################################################################################
#Exclude phenology varieties with no Marker data
phenology = phenology %>% filter(!(phenology$Variety %in% removed_varieties))

rm(marker_varieties, phenotype_varieties, union.vars, all_varieties, removed_varieties)

######################################################################################################################################
##Transpose, set colnames to SNP
MAF = Marker_Data[,-c(1:7)] %>% t()
colnames(MAF) = Marker_Data$ID

###############################################
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
