######################################################################################################################################
##Merge phenology and comb.env into phen.env -------------------- add soil info to phenotype matrix
##Note: phenology filtered by var genotype ^ phenotype; env.comb contains site soil
phen.env = phenology %>% inner_join(env.comb)

######################################################################################################################################
##BOM data preprocessing -------------------- make all BOM things into weather.BOM matrix
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
#organise temp to 3 classes: <4, >30, temperate, 
weather.BOM$min.temp.cut = cut(weather.BOM$`Minimum temperature (Degree C)`,
                                breaks = c(-Inf, 4, 30, Inf),
                                labels = c("< 4", "Temperate", "> 30")
)

weather.BOM$max.temp.cut = cut(weather.BOM$`Maximum temperature (Degree C)`,
                                breaks = c(-Inf, 4, 30, Inf),
                                labels = c("< 4", "Temperate", "> 30")
)

weather.BOM$temp.range = weather.BOM$`Maximum temperature (Degree C)` - weather.BOM$`Minimum temperature (Degree C)`

weather.BOM$date = as.Date(paste(weather.BOM$Day, weather.BOM$Month, as.character(weather.BOM$Year), sep = '-'), "%d-%m-%Y")

weather.BOM = weather.BOM %>% select(-c("Year", "Month", "Day", "Quality", 
                                         "Period over which rainfall was measured (days)",
                                         "Days of accumulation of minimum temperature", 
                                         "Days of accumulation of maximum temperature"))

######################################################################################################################################
##Summarise weather over the period of growth for an individual plant.
#Including days prior to planting as may be relevant
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

######################################################################################################################################
##get the final phenology matrix ready, with weather data combined ------- phen.env contains phenotype x staged env x soil
phen.env = data.frame(phen.env, gs.weather_results)

############################################################
#Add beginning photoperiod for weather data
phen.env$ppd = daylength(phen.env$Latitude, phen.env$DOP)
phen.env$ppd[which(phen.env$Trial == "18Hrs")] = 18

######################################################################################################################################
