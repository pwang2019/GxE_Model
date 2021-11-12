####################################################################################################################################################
##Trial site details ----  Date of planting is important!
Trial15_16 = read_excel("compiled_data/Trials.xlsx", trim_ws = TRUE)

trial_data = Trial15_16

###############################################################################
#get DOP into phenology
colnames(trial_data)[1:2] = c("Year", "BOM.Stn")
trial_data$Year = as.character(trial_data$Year)
trial_data$DOP = as.Date(trial_data$DOP)
phenology$DOP = as.Date(phenology$DOP)

temp = phenology %>% left_join(trial_data[,1:11],
		by = c("Year", "BOM.Stn", "Location", "DOP"))

temp[temp$Location == "STHPER", "Location"] = "SoP"

phen.env = temp

######################################################################################################################################
##Merge phenology and comb.env into phen.env--------------------
#Note: phenology filtered by var genotype ^ phenotype; env.comb contains site soil
phen.env = phen.env %>% inner_join(env.comb)

######################################################################################################################################
##BOM data --------
#recode station 0090901 to 009225
BOM.solar$`Bureau of Meteorology station number` = plyr::revalue(BOM.solar$`Bureau of Meteorology station number`,
                                                                  c("009091" = "009225"))
#####################################################
#join temp, solar, rain to a single matrix
weather.BOM = full_join(BOM.rain %>% select(-`Product code`), 
                         BOM.solar %>% select(-`Product code`))
weather.BOM = full_join(weather.BOM, BOM.temp)

rm(BOM.rain, BOM.solar, BOM.temp)

#####################################################
weather.BOM$date = as.Date(paste(weather.BOM$Day, weather.BOM$Month, as.character(weather.BOM$Year), sep = '-'), "%d-%m-%Y")

weather.BOM = weather.BOM %>% select(-c("Year", "Month", "Day", "Quality", 
                                         "Period over which rainfall was measured (days)",
                                         "Days of accumulation of minimum temperature", 
                                         "Days of accumulation of maximum temperature"))
weather.BOM = as.data.frame(weather.BOM)

#modelling_data will have NA by BOM some measurements
na.index = rep(FALSE, nrow(weather.BOM))
for(m in 1:nrow(weather.BOM)){
    if(sum(is.na(weather.BOM[m, ])>0))
	    na.index[m]=TRUE
}
weather.BOM.f = weather.BOM[!na.index, ]

rm(na.index, m)

##############################################################################################################################################################
###############################################################################
#Summarise weather over the period of growth for an individual plant. Including days prior to planting as may be relevant filter weather data to trial and year
#get weather summary on DAILY basis not stage here, write summary to variety at each year each location
#firstly, need to generate an empty data.frame holder for daily weather data, from DOP to 11-30
#weather var per day is: "Rainfall amount (millimetres)", "Daily global solar exposure (MJ/m*m)", "Minimum temperature (Degree C)", "Maximum temperature (Degree C)", "temp.range", "date"
temp = as.data.frame(trial_data)
as.Date("2015-11-30") - as.Date(temp[1:3, "DOP"])	#186 194 201
as.Date("2016-11-30") - as.Date(temp[4:10, "DOP"])	#203 182 198 198 195 196 196

###############################################################################
#generate an empty matrix for full BOM ENV
#NOTE: to avoid NA problems, use 180 days as the complete grow season !!!!!!! maybe 210 days
temp = matrix(rep(NA, 180*4), nrow=1)
colnames(temp) = as.character(1:(180*4))
for(m in 1:180){
    colnames(temp)[((m-1)*4+1):(m*4)] = paste("D", as.character(m), 
	  c("Rainfall amount (millimetres)", "Daily global solar exposure (MJ/m*m)", "Minimum temperature (Degree C)", "Maximum temperature (Degree C)"), #, "temp.range"), #, "date"),
	  sep="__")
}

tt = colnames(temp)
gs.weather_results = data.frame(temp)
colnames(gs.weather_results) = tt

rm(temp,tt)

###############################################################################
#now assign each BOM record to the matrix
phen.env.f = phen.env

for(individual in 1:nrow(phen.env.f)){
  #Filter weather data for individual
  individual_weather = weather.BOM %>% filter(`Bureau of Meteorology station number` == phen.env.f$BOM.Stn[individual],
                                               date > as.Date(phen.env.f$DOP[individual]), 
											   date <= as.Date(phen.env.f$DOP[individual] + 180))	#NOTE, specific row for traitID row can be NA!
  
  gs_weather_df = individual_weather[1, ] %>% select(-c("Bureau of Meteorology station number", "date"))
  colnames(gs_weather_df)[1:4] = paste("1", colnames(gs_weather_df)[1:4], sep="__")
  
  for(m in 2:nrow(individual_weather)){
      tt = individual_weather[m, ] %>% select(-c("Bureau of Meteorology station number", "date"))
      colnames(tt)[1:4] = paste(as.character(m), colnames(tt)[1:4], sep="__")
      gs_weather_df = cbind(gs_weather_df, tt)
  }
  
  gs.weather_results[individual, 1:ncol(gs_weather_df)] = gs_weather_df 
}

###############################################################################
#NOTE: BOM weather data still has missing values, plus mapping everyday weather data create new NA
#use KNN imputation (VIM package) to fill these NA 
tt = apply(gs.weather_results, 2, is.na)
tt = apply(tt, 2, sum)

#filter out those columns where missingness > 20%, one column corresponds to one weather feature of a specific day
gs.weather_results.knn = gs.weather_results[, which(tt<=as.integer(nrow(gs.weather_results)*0.2))]

#now VIM knn imputation
gs.weather_results.knn = kNN(gs.weather_results.knn)[1:ncol(gs.weather_results)]	

###############################################################################
#now append phenotype and env all BOM vars
phen.env.f = data.frame(phen.env.f, gs.weather_results.knn) #%>% 
  #filter(!is.na(traitID))

###############################################################################
#save these time-consuming variables
save(phen.env, phen.env.f, gs.weather_results, gs.weather_results.knn, file="WeatherVar_AMarker_BOM.RData")

###############################################################################
#some cleanup
rm(gs_weather_df, individual, individual_weather, tt)
rm(BOM.rain, BOM.solar, BOM.temp)

##############################################################################################################################################################


















if(FALSE){
######################################################################################################################################
#Summarise weather over the period of growth for an individual plant.
#Including days prior to planting as may be relevant

#for sample
# filter weather data to trial and year
# get stage transition dates (DOP + stage*ZS49days/3)
# for stage
#  filter weather to between current and previous stage
#  get weather summary
#write summary to sample

load(file = "rdata/data/gs_columns.RData")
gs.weather_results = as.data.frame(matrix(NA, nrow = 1, ncol = length(gs_columns)))
colnames(gs.weather_results) = gs_columns

DAYS_AFTER_DOP = 120
GROWTH_STAGES = 4

for(individual in 1:nrow(phen.env)){
  #Filter weather data for individual
  individual_weather = weather.BOM %>% filter(`Bureau of Meteorology station number` == phen.env$BOM.Stn[individual],
                                               date > as.Date(phen.env$DOP[individual]),
                                               date < as.Date(phen.env$DOP[individual]) + DAYS_AFTER_DOP)
  #Stage transition dates
  date_seed = phen.env$DOP[individual]
  date_gs1 = as.Date(phen.env$DOP[individual]) + round(1 * DAYS_AFTER_DOP/GROWTH_STAGES)
  date_gs2 = as.Date(phen.env$DOP[individual]) + round(2 * DAYS_AFTER_DOP/GROWTH_STAGES)
  date_gs3 = as.Date(phen.env$DOP[individual]) + round(3 * DAYS_AFTER_DOP/GROWTH_STAGES)
  date_gs4 = as.Date(phen.env$DOP[individual]) + round(4 * DAYS_AFTER_DOP/GROWTH_STAGES)

  
  individual_weather$growth_stage = cut(individual_weather$date,
                                         breaks = c(date_seed, date_gs1, date_gs2, date_gs3, date_gs4),
                                         labels = c("gs1", "gs2", "gs3", "gs4")
  )
  
  
  gs_weather_df = individual_weather %>%
    group_by(growth_stage) %>% 
    summarise("rainfall_median"  = median(`Rainfall amount (millimetres)` %>% na.omit()),
              "rainfall_mean"    = mean(`Rainfall amount (millimetres)` %>% na.omit()),
              "rainfall_min"     = min(`Rainfall amount (millimetres)` %>% na.omit()),
              "rainfall_max"     = max(`Rainfall amount (millimetres)` %>% na.omit()),
              "rainfall_sum"     = sum(`Rainfall amount (millimetres)` %>% na.omit()),
              "solar_mean"       = mean(`Daily global solar exposure (MJ/m*m)`%>% na.omit()),
              "solar_min"        = min(`Daily global solar exposure (MJ/m*m)` %>% na.omit()),
              "solar_max"        = max(`Daily global solar exposure (MJ/m*m)` %>% na.omit()),
              "solar_sum"        = max(`Daily global solar exposure (MJ/m*m)` %>% na.omit()),
              "mintemp_min"      = min(`Minimum temperature (Degree C)` %>% na.omit()),
              "mintemp_max"      = max(`Minimum temperature (Degree C)` %>% na.omit()),
              "mintemp_mean"     = mean(`Minimum temperature (Degree C)` %>% na.omit()),
              "mintemp_range"    = max(`Minimum temperature (Degree C)` %>% na.omit()) - min(`Minimum temperature (Degree C)` %>% na.omit()),
              "maxtemp_min"      = min(`Maximum temperature (Degree C)` %>% na.omit()),
              "maxtemp_max"      = max(`Maximum temperature (Degree C)` %>% na.omit()),
              "maxtemp_mean"     = mean(`Maximum temperature (Degree C)` %>% na.omit()),
              "maxtemp_range"    = max(`Maximum temperature (Degree C)` %>% na.omit()) - min(`Maximum temperature (Degree C)` %>% na.omit()),
              "temp_under4"    = sum(c(max.temp.cut %>% na.omit() == "< 4",       min.temp.cut %>% na.omit() == "< 4"))       / (2*length(max.temp.cut %>% na.omit())),
              "temp_Temperate" = sum(c(max.temp.cut %>% na.omit() == "Temperate", min.temp.cut %>% na.omit() == "Temperate")) / (2*length(max.temp.cut %>% na.omit())),
              "temp_over30"    = sum(c(max.temp.cut %>% na.omit() == "> 30",      min.temp.cut %>% na.omit() == "> 30"))      / (2*length(max.temp.cut %>% na.omit())),
              "temp_mean_range"  = mean(temp.range %>% na.omit())
    )
  
  gs_values = as.matrix(gs_weather_df[,-1]) %>% as.vector() %>% as.numeric()
  gs.weather_results[individual,] = gs_values 
  
}

#############################################################
#get the final phenology matrix ready, with weather data combined
phen.env = data.frame(phen.env, gs.weather_results)

  #IF YOU CHANGE ANYTHING, RERUN THIS TO REBUILD THE DATAFRAME (reuben is a bad coder)
  rn = as.vector(gs_weather_df$growth_stage)
  cn = as.vector(colnames(gs_weather_df[,-1]))
  gs_columns = outer(rn, cn, paste, sep="_") %>% as.vector()

#Add beginning photoperiod for weather data
phen.env$ppd = daylength(phen.env$Latitude, phen.env$DOP)
phen.env$ppd[which(phen.env$Trial == "18Hrs")] = 18

}

