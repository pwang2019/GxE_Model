######################################################################################################################################
##BOM data ----------------------
BOM.temp = list.files(path = "compiled_data/Weather/Temperature", full.names = TRUE) %>% 
  lapply(read_csv) %>%
  bind_rows()

B.min = BOM.temp %>% 
  filter(`Product code` == "IDCJAC0011") %>%
  select(-c(`Product code`, `Maximum temperature (Degree C)`, 
            `Days of accumulation of maximum temperature`, Quality))

B.max = BOM.temp %>% 
  filter(`Product code` == "IDCJAC0010") %>%
  select(-c(`Product code`, `Minimum temperature (Degree C)`, 
            `Days of accumulation of minimum temperature`, Quality)) 

BOM.temp = full_join(B.min, B.max)
BOM.temp$Year = as.character(BOM.temp$Year)

rm(B.min)
rm(B.max)

###################################################################
BOM.solar = list.files(path = "compiled_data/Weather/Solar", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows
BOM.solar$Year = as.character(BOM.solar$Year)

BOM.rain = list.files(path = "compiled_data/Weather/Rainfall", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows
BOM.rain$Year = as.character(BOM.rain$Year)

######################################################################################################################################
#SLGA data pull ----------------------------
load("rdata/data/env.comb.RData")

env.comb$Trial.year = as.character(env.comb$Trial.year)
env.comb$Trial[is.na(env.comb$Trial)] = "Normal"
env.comb$DOP = as.Date(env.comb$DOP)
env.comb$Trial = env.comb$Trial %>% 
  recode("NL" = "NatL")
env.comb$Location = env.comb$Location %>% 
  recode("SoP" = "STHPER")

env.comb = env.comb %>% 
  select(-State, -Seeding.rate..kg.ha.) %>% 
  rename("Year" = Trial.year)

#Soil data comes out as text for some reason
env.comb[, 11:19] =  sapply(env.comb[, 11:19], as.numeric)

#Create soil classes vector:
soildata = data.frame(
    "CLAY" = env.comb$CLY,
    "SILT" = env.comb$SLT,
    "SAND" = env.comb$SND,
    "OC"   = NA)

#Creates factor variable for soil type
env.comb$soiltype = soiltexture::TT.points.in.classes(soildata,
                                                       class.sys = "AU.TT",
                                                       PiC.type = "t",
                                                       collapse = "_")
													   
######################################################################################################################################													   
													   