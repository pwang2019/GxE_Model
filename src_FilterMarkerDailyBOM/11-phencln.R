######################################################################################################################################
##Load in Phenology main data ----
phenology = read_excel("compiled_data/Reuben_Phenology_All Sites_ 14-16_FIELD Data.xlsx",
                        col_types = c("numeric", "text", "text", "text", "text", 
                                      "text", "date", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric"),
                        trim_ws = TRUE)

colnames(phenology) = c("Sort", "Variety", "Reason", "Year", "Trial", "Location", "DOP", "ZS49PlHt", 
                         "ZS49.days", "HrvPlHt", "ZS91.days", "PPd.days", 
                         "RN", "GrYld.kg.ha", "PrdGrYld.kg.ha", "PrdGrYld.t.ha", 
                         "SEGrYld", "TlrNo", "PlntNo")
						 
######################################################################################################################################
##Filter out trials missing environmental data, and duplicates
##2014 trial no records, TOP trials not well recorded
phenology = phenology %>% 
  select(-Sort) %>% 
  distinct() %>% 
  filter(Year != "2014", 
         !(Trial %in% c("TOP1", "TOP2", "TOP3")))
phenology$Trial = phenology$Trial %>% 
  recode("1NI" = "NI",
         "2Ir" = "Ir")
phenology$Trial = fct_explicit_na(phenology$Trial, "Normal") %>% as.character()
  
######################################################################################################################################
##Make Year numeric for later merges
#phenology$Year = as.numeric(phenology$Year)
#phenology$Trial[is.na(phenology$Trial)] = "NA"
#phenology$ZS49.days = as.integer(phenology$ZS49.days)

phenology$DOP = as.Date(phenology$DOP)

#Add the BOM.stn variable
field = phenology %>% 
  select(Year, Location, Trial) %>% 
  unique()

field$BOM.Stn =  c("009225", # 2016 SoP
                   "009225", # 2016 SoP
                   "008315", # 2015 GER
                   "008315", # 2016 GER
                   "010916", # 2015 KAT
                   "010916", # 2016 KAT
                   "009542", # 2015 ESP
                   "009542", # 2016 ESP                   
                   "010092", # 2016 MER
                   "010092") # 2016 MER
field$trial_index = 1:length(field$BOM.Stn)
field$irrigation = !(field$Trial == "Normal" | field$Trial == "Ir")
#setDT(field, keep.rownames = "trial_index")

phenology = phenology %>% 
  merge(field, by = c("Year", "Location", "Trial"))
######################################################################################################################################
