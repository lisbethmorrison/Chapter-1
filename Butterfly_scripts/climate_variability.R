#########################################################
## Title: Climate variability between seasons for early & late comparison butterflies
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: March 2018
##########################################################

rm(list=ls()) # clear R

## read in climate data
met_office_temp <- read.csv("../Data/Temp_data/Met_Office_Seasonal_Temps.csv", header=T)
rainfall_data <- read.table("../Data/Temp_data/UK_rainfall_data.txt", header=TRUE) # add climate data

## create 2 data frames: one with early years (1980-2004) and one with late years (1995-2016) and rename columns
climate_var_early <- met_office_temp[met_office_temp$year<=2004,]
climate_var_early2 <- met_office_temp[(met_office_temp$year<=2000) & (met_office_temp$year>=1985),]
climate_var_early$year<- "Early"
climate_var_early2$year<- "Early"
climate_var_late <- met_office_temp[met_office_temp$year>=1995,]
climate_var_late2 <- met_office_temp[(met_office_temp$year<=2012) & (met_office_temp$year>=2000),]
climate_var_late$year<- "Late"
climate_var_late2$year<- "Late"

## rbind the two together
climate_var <- rbind(climate_var_early, climate_var_late)
## NOTE: years are duplicated between early and late (10 years)
climate_var2 <- rbind(climate_var_early2, climate_var_late2)
## these years have no duplication

## F test between early and late for each variable
var.test(mean_MAM ~ year, data=climate_var) # mean spring
var.test(mean_JJA ~ year, data=climate_var) # mean summer
var.test(mean_SON ~ year, data=climate_var) # mean autumn
var.test(mean_DJF ~ year, data=climate_var) # mean winter
### all non-significant
var.test(max_MAM ~ year, data=climate_var) # max spring
var.test(max_JJA ~ year, data=climate_var) # max summer
var.test(max_SON ~ year, data=climate_var) # max autumn
var.test(max_DJF ~ year, data=climate_var) # max winter
### all non-significant 

## therefore variances of the two variables are approximately equal (NOT significantly different)

## same tests as above for the no overlap years
## F test between early and late for each variable
var.test(mean_MAM ~ year, data=climate_var2) # mean spring
var.test(mean_JJA ~ year, data=climate_var2) # mean summer
var.test(mean_SON ~ year, data=climate_var2) # mean autumn
var.test(mean_DJF ~ year, data=climate_var2) # mean winter
### all non-significant
var.test(max_MAM ~ year, data=climate_var2) # max spring
var.test(max_JJA ~ year, data=climate_var2) # max summer
var.test(max_SON ~ year, data=climate_var2) # max autumn
var.test(max_DJF ~ year, data=climate_var2) # max winter
### all non-significant 

##################### RAINFALL DATA #######################
## remove monthly columns, just leave seasons
rainfall_data <- rainfall_data[-c(2:13,18)]

## create 2 data frames: one with early years (1980-2004) and one with late years (1995-2016) and rename columns
rain_var_early <- rainfall_data[(rainfall_data$Year<=2004) & (rainfall_data$Year>=1980),] # 1980-2004
rain_var_early$Year<- "Early"
rain_var_early2 <- rainfall_data[(rainfall_data$Year<=2000) & (rainfall_data$Year>=1985),] # 1985-2000
rain_var_early2$Year<- "Early"
rain_var_late <- rainfall_data[(rainfall_data$Year>=1995) & (rainfall_data$Year<=2016),] # 1995-2016
rain_var_late$Year<- "Late"
rain_var_late2 <- rainfall_data[(rainfall_data$Year<=2012) & (rainfall_data$Year>=2000),] # 2000-2012
rain_var_late2$Year<- "Late"

## rbind the two together
rain_var <- rbind(rain_var_early, rain_var_late)
## NOTE: years are duplicated between early and late (10 years)
rain_var2 <- rbind(rain_var_early2, rain_var_late2)
## these years have no duplication

## for some reason winter is a factor
rain_var$WIN <- as.numeric(rain_var$WIN)
rain_var2$WIN <- as.numeric(rain_var2$WIN)

## F test between early and late for each variable
var.test(SPR ~ Year, data=rain_var) # mean spring
var.test(SUM ~ Year, data=rain_var) # mean summer
var.test(AUT ~ Year, data=rain_var) # mean autumn
var.test(WIN ~ Year, data=rain_var) # mean winter
### all non-significant
### same but with the non-overlapping years
var.test(SPR ~ Year, data=rain_var2) # mean spring
var.test(SUM ~ Year, data=rain_var2) # mean summer
var.test(AUT ~ Year, data=rain_var2) # mean autumn
var.test(WIN ~ Year, data=rain_var2) # mean winter
### all non-significant 

# qqnorm(rain_var2$WIN, pch = 1, frame = FALSE)
# qqline(rain_var2$WIN, col = "steelblue", lwd = 2)
# qqnorm(rain_var2$SPR, pch = 1, frame = FALSE)
# qqline(rain_var2$SPR, col = "steelblue", lwd = 2)
# qqnorm(rain_var2$SUM, pch = 1, frame = FALSE)
# qqline(rain_var2$SUM, col = "steelblue", lwd = 2)
# qqnorm(rain_var2$AUT, pch = 1, frame = FALSE)
# qqline(rain_var2$AUT, col = "steelblue", lwd = 2)

## therefore variances of the two variables are approximately equal (NOT significantly different)


####################### EXTREME CLIMATIC EVENTS ######################

##### MEAN TEMPERATURE

mean_winter_median <- median(met_office_temp$mean_DJF) # 4.6
mean_spring_median <- median(met_office_temp$mean_MAM) # 9
mean_summer_median <- median(met_office_temp$mean_JJA) # 15.8
mean_autumn_median <- median(met_office_temp$mean_SON) # 10.7

### calculate median absolute deviation for each mean temperature
mean_winter_MAD <- mad(met_office_temp$mean_DJF, center = median(met_office_temp$mean_DJF), constant = 1.4826,
                       na.rm = FALSE, low = FALSE, high = FALSE) # 1.18
mean_spring_MAD <- mad(met_office_temp$mean_MAM, center = median(met_office_temp$mean_MAM), constant = 1.4826,
                       na.rm = FALSE, low = FALSE, high = FALSE) # 0.89
mean_summer_MAD <- mad(met_office_temp$mean_JJA, center = median(met_office_temp$mean_JJA), constant = 1.4826,
                       na.rm = FALSE, low = FALSE, high = FALSE) # 0.74
mean_autumn_MAD <- mad(met_office_temp$mean_SON, center = median(met_office_temp$mean_SON), constant = 1.4826,
                       na.rm = FALSE, low = FALSE, high = FALSE) # 0.59

#### calculate extreme years for mean temperature for each season
extreme_winter <- met_office_temp[(met_office_temp$mean_DJF>=(mean_winter_median+(2*mean_winter_MAD))) | 
                    (met_office_temp$mean_DJF<=(mean_winter_median-(2*mean_winter_MAD))),] ## no years
extreme_spring <- met_office_temp[(met_office_temp$mean_MAM>=(mean_spring_median+(2*mean_spring_MAD))) | 
                    (met_office_temp$mean_MAM<=(mean_spring_median-(2*mean_spring_MAD))),] ## 2013 (negative)
extreme_summer <- met_office_temp[(met_office_temp$mean_JJA>=(mean_summer_median+(2*mean_summer_MAD))) | 
                    (met_office_temp$mean_JJA<=(mean_summer_median-(2*mean_summer_MAD))),] ## 1995 and 2003 (both positive)
extreme_autumn <- met_office_temp[(met_office_temp$mean_SON>=(mean_autumn_median+(2*mean_autumn_MAD))) | 
                    (met_office_temp$mean_SON<=(mean_autumn_median-(2*mean_autumn_MAD))),] ## 1992, 1993 (both negative), 2006, 2011, 2014 (all positive)

##### MAX TEMPERATURE

max_winter_median <- median(met_office_temp$max_DJF) # 7.4
max_spring_median <- median(met_office_temp$max_MAM) # 13.1
max_summer_median <- median(met_office_temp$max_JJA) # 20.2
max_autumn_median <- median(met_office_temp$max_SON) # 14.1

### calculate median absolute deviation for each max temperature
max_winter_MAD <- mad(met_office_temp$max_DJF, center = median(met_office_temp$max_DJF), constant = 1.4826,
                       na.rm = FALSE, low = FALSE, high = FALSE) # 1.48
max_spring_MAD <- mad(met_office_temp$max_MAM, center = median(met_office_temp$max_MAM), constant = 1.4826,
                       na.rm = FALSE, low = FALSE, high = FALSE) # 1.04
max_summer_MAD <- mad(met_office_temp$max_JJA, center = median(met_office_temp$max_JJA), constant = 1.4826,
                       na.rm = FALSE, low = FALSE, high = FALSE) # 1.33
max_autumn_MAD <- mad(met_office_temp$max_SON, center = median(met_office_temp$max_SON), constant = 1.4826,
                       na.rm = FALSE, low = FALSE, high = FALSE) # 0.89

extreme_max_winter <- met_office_temp[(met_office_temp$max_DJF>=(max_winter_median+(2*max_winter_MAD))) | 
                        (met_office_temp$max_DJF<=(max_winter_median-(2*max_winter_MAD))),] ## no years
extreme_max_spring <- met_office_temp[(met_office_temp$max_MAM>=(max_spring_median+(2*max_spring_MAD))) | 
                        (met_office_temp$max_MAM<=(max_spring_median-(2*max_spring_MAD))),] ## 1986, 2013 (negative)
extreme_max_summer <- met_office_temp[(met_office_temp$max_JJA>=(max_summer_median+(2*max_summer_MAD))) | 
                        (met_office_temp$max_JJA<=(max_summer_median-(2*max_summer_MAD))),] ## no years
extreme_max_autumn <- met_office_temp[(met_office_temp$max_SON>=(max_autumn_median+(2*max_autumn_MAD))) | 
                        (met_office_temp$max_SON<=(max_autumn_median-(2*max_autumn_MAD))),] ## 1993 (negative), 2006, 2011 (both positive)

##### MIN TEMPERATURE
min_winter_median <- median(met_office_temp$min_DJF) # 1.9
min_spring_median <- median(met_office_temp$min_MAM) # 4.9
min_summer_median <- median(met_office_temp$min_JJA) # 11.3
min_autumn_median <- median(met_office_temp$min_SON) # 7.3

### calculate median absolute deviation for each min temperature
min_winter_MAD <- mad(met_office_temp$min_DJF, center = median(met_office_temp$min_DJF), constant = 1.4826,
                      na.rm = FALSE, low = FALSE, high = FALSE) # 1.03
min_spring_MAD <- mad(met_office_temp$min_MAM, center = median(met_office_temp$min_MAM), constant = 1.4826,
                      na.rm = FALSE, low = FALSE, high = FALSE) # 0.74
min_summer_MAD <- mad(met_office_temp$min_JJA, center = median(met_office_temp$min_JJA), constant = 1.4826,
                      na.rm = FALSE, low = FALSE, high = FALSE) # 0.44
min_autumn_MAD <- mad(met_office_temp$min_SON, center = median(met_office_temp$min_SON), constant = 1.4826,
                      na.rm = FALSE, low = FALSE, high = FALSE) # 0.74

extreme_min_winter <- met_office_temp[(met_office_temp$min_DJF>=(min_winter_median+(2*min_winter_MAD))) | 
                        (met_office_temp$min_DJF<=(min_winter_median-(2*min_winter_MAD))),] ## 1982 (negative)
extreme_min_spring <- met_office_temp[(met_office_temp$min_MAM>=(min_spring_median+(2*min_spring_MAD))) | 
                        (met_office_temp$min_MAM<=(min_spring_median-(2*min_spring_MAD))),] ## 1984 and 2013 (both negative)
extreme_min_summer <- met_office_temp[(met_office_temp$min_JJA>=(min_summer_median+(2*min_summer_MAD))) | 
                        (met_office_temp$min_JJA<=(min_summer_median-(2*min_summer_MAD))),] ## 2003, 2006 (both positive), 2011 (negative)
extreme_min_autumn <- met_office_temp[(met_office_temp$min_SON>=(min_autumn_median+(2*min_autumn_MAD))) | 
                        (met_office_temp$min_SON<=(min_autumn_median-(2*min_autumn_MAD))),] ## 1993 (negative), 2006, 2011 (both positive)

### 2013 (negative extreme event) is common for spring mean min and max
### 1982 (negative extreme event) for winter min (down to -0.9)



