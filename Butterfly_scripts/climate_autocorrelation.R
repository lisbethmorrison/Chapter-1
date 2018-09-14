###########################################################
## Title: Temporal trend in spatial autocorrelation in climate
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: September 2017
##########################################################


#### getwd() = ("C:/Users/wh890746/Documents/Mini_project/R_script") ####

rm(list=ls()) # clear R

###### CLIMATE DATA CLEANING ######

## add data
climate.data <- read.table("../Data/Temp_data/all.sites.CIP.mean.temp.AND.rainfall.1971_2012_running.3month.means.txt", header=TRUE) # add climate data
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header = TRUE) # pair attr site data

## remove years before 1980 and remove occasional 2013 data
climate.data <- climate.data[climate.data$year>1979,]
climate.data <- climate.data[climate.data$year<2013,]

## remove mean.temp and rainfall columns (i.e. only keep the running 3 month values)
climate.data = subset(climate.data, select = -c(mean.temp, rainfall))

## remove months 2,3,5,6,8,9,11 (1=winter, 4=spring, 7=summer, 10=autumn) 
climate.data <- subset(climate.data,!month%in%c(2,3,5,6,8,9,11,12))

## change month to season ##
colnames(climate.data)[3] <- "season"

## 1=winter, 4=spring, 7=summer, 10=autumn
climate.data$season[climate.data$season==1]<-"winter"
climate.data$season[climate.data$season==4]<-"spring"
climate.data$season[climate.data$season==7]<-"summer"
climate.data$season[climate.data$season==10]<-"autumn"

## produce table with spring_temp, spring_rainfall, summer_temp, summer_rainfall etc...
## create empty dataframe
climate.autocorr <- data.frame(site=character(68805), year=character(68805),spring_temp=numeric(68805))

## add site and year data 
climate.autocorr$site <- (climate.data$site[climate.data$season=="spring"])
climate.autocorr$year <- (climate.data$year[climate.data$season=="spring"])

## create seasonal mean temperatures
climate.autocorr$spring_temp <- (climate.data$running.3m.mean.temp[climate.data$season=="spring"])
climate.autocorr$summer_temp <- (climate.data$running.3m.mean.temp[climate.data$season=="summer"])
climate.autocorr$autumn_temp <- (climate.data$running.3m.mean.temp[climate.data$season=="autumn"])
climate.autocorr$winter_temp <- (climate.data$running.3m.mean.temp[climate.data$season=="winter"])

## create seasonal mean rainfall 
climate.autocorr$spring_rain <- (climate.data$running.3m.mean.rain[climate.data$season=="spring"])
climate.autocorr$summer_rain <- (climate.data$running.3m.mean.rain[climate.data$season=="summer"])
climate.autocorr$autumn_rain <- (climate.data$running.3m.mean.rain[climate.data$season=="autumn"])
climate.autocorr$winter_rain <- (climate.data$running.3m.mean.rain[climate.data$season=="winter"])

###### ADD SITE EASTING AND NORTHING DATA ######
site_data <- pair_attr[,c(2,3,9,10,11,12)] ## include site data columns (site names, easting and northing)
site_names <- as.data.frame(unique(c(site_data$site1, site_data$site2)))
names(site_names) <- "site_code"

site_names <- merge(site_names, site_data, by.x="site_code", by.y="site1")
## keep columns needed
site_names <- site_names[,c(1,3,4)]
site_names <- unique(site_names) # 673 sites with climate data

## merge site_names with climate.autocorr
climate.autocorr <- merge(climate.autocorr, site_names, by.x="site", by.y="site_code")

## remove any rows with NAs (i.e sites without data)
climate.autocorr<-na.omit(climate.autocorr)

## change column names
colnames(climate.autocorr)[11] <- "easting"
colnames(climate.autocorr)[12] <- "northing"


#########################
## calculate moran's I ##
#########################

library(ape)

############################
#### SPRING TEMPERATURE ####
############################

## loop through year and calculate moran's i spring temp for each year
results.temp <- NULL
moran.spring.temp <- NULL
moran.final <- NULL
for (i in sort(unique(climate.autocorr$year))){
  climate.autocorr.year <- climate.autocorr[climate.autocorr$year==i,]
 
  climate.dist <- as.matrix(dist(cbind(climate.autocorr.year$easting, climate.autocorr.year$northing)))
  
  climate.dist.inv <- 1/climate.dist
  diag(climate.dist.inv) <- 0
  climate.dist.inv[is.infinite(climate.dist.inv)] <- 0

moran.i.result <- Moran.I(climate.autocorr.year$spring_temp, climate.dist.inv)

results.temp <- data.frame(moran.i.result$observed,i)
colnames(results.temp)[1] <- "moran.i.spring.temp"
colnames(results.temp)[2] <- "year"
moran.spring.temp <- rbind(moran.spring.temp, results.temp)

}

moran.final <- rbind(moran.final, moran.spring.temp)

############################
#### SUMMER TEMPERATURE ####
############################

## loop through year and calculate moran's i summer temp for each year
results.temp <- NULL
moran.summer.temp <- NULL
for (i in sort(unique(climate.autocorr$year))){
  climate.autocorr.year <- climate.autocorr[climate.autocorr$year==i,]
  
  climate.dist <- as.matrix(dist(cbind(climate.autocorr.year$easting, climate.autocorr.year$northing)))
  
  climate.dist.inv <- 1/climate.dist
  diag(climate.dist.inv) <- 0
  climate.dist.inv[is.infinite(climate.dist.inv)] <- 0
  
  moran.i.result <- Moran.I(climate.autocorr.year$summer_temp, climate.dist)
  moran.i.result<-as.data.frame(moran.i.result)
  
  results.temp <- data.frame(moran.i.result$observed,i)
  colnames(results.temp)[1] <- "moran.i.summer.temp"
  colnames(results.temp)[2] <- "year"
  moran.summer.temp <- rbind(moran.summer.temp, results.temp)
  
}

moran.final <- merge(moran.final, moran.summer.temp, by.x="year", by.y="year")

############################
#### AUTUMN TEMPERATURE ####
############################

## loop through year and calculate moran's i autumn temp for each year
results.temp <- NULL
moran.autumn.temp <- NULL
for (i in sort(unique(climate.autocorr$year))){
  climate.autocorr.year <- climate.autocorr[climate.autocorr$year==i,]
  
  climate.dist <- as.matrix(dist(cbind(climate.autocorr.year$easting, climate.autocorr.year$northing)))
  
  climate.dist.inv <- 1/climate.dist
  diag(climate.dist.inv) <- 0
  climate.dist.inv[is.infinite(climate.dist.inv)] <- 0
  
  moran.i.result <- Moran.I(climate.autocorr.year$autumn_temp, climate.dist)
  moran.i.result<-as.data.frame(moran.i.result)
  
  results.temp <- data.frame(moran.i.result$observed,i)
  colnames(results.temp)[1] <- "moran.i.autumn.temp"
  colnames(results.temp)[2] <- "year"
  moran.autumn.temp <- rbind(moran.autumn.temp, results.temp)
  
}

moran.final <- merge(moran.final, moran.autumn.temp, by.x="year", by.y="year")

############################
#### WINTER TEMPERATURE ####
############################

## loop through year and calculate moran's i winter temp for each year
results.temp <- NULL
moran.winter.temp <- NULL
for (i in sort(unique(climate.autocorr$year))){
  climate.autocorr.year <- climate.autocorr[climate.autocorr$year==i,]
  
  climate.dist <- as.matrix(dist(cbind(climate.autocorr.year$easting, climate.autocorr.year$northing)))
  
  climate.dist.inv <- 1/climate.dist
  diag(climate.dist.inv) <- 0
  climate.dist.inv[is.infinite(climate.dist.inv)] <- 0
  
  moran.i.result <- Moran.I(climate.autocorr.year$winter_temp, climate.dist)
  moran.i.result<-as.data.frame(moran.i.result)
  
  results.temp <- data.frame(moran.i.result$observed,i)
  colnames(results.temp)[1] <- "moran.i.winter.temp"
  colnames(results.temp)[2] <- "year"
  moran.winter.temp <- rbind(moran.winter.temp, results.temp)
  
}

moran.final <- merge(moran.final, moran.winter.temp, by.x="year", by.y="year")


#########################
#### SPRING RAINFALL ####
#########################

## loop through year and calculate moran's i spring rainfall for each year
results.temp <- NULL
moran.spring.rain <- NULL

for (i in sort(unique(climate.autocorr$year))){
  climate.autocorr.year <- climate.autocorr[climate.autocorr$year==i,]
  
  climate.dist <- as.matrix(dist(cbind(climate.autocorr.year$easting, climate.autocorr.year$northing)))
  
  climate.dist.inv <- 1/climate.dist
  diag(climate.dist.inv) <- 0
  climate.dist.inv[is.infinite(climate.dist.inv)] <- 0 
  
  moran.i.result <- Moran.I(climate.autocorr.year$spring_rain, climate.dist)
  moran.i.result<-as.data.frame(moran.i.result)
  
  results.temp <- data.frame(moran.i.result$observed,i)
  colnames(results.temp)[1] <- "moran.i.spring.rain"
  colnames(results.temp)[2] <- "year"
  moran.spring.rain <- rbind(moran.spring.rain, results.temp)
  
}

moran.final <- merge(moran.final, moran.spring.rain, by.x="year", by.y="year")

#########################
#### SUMMER RAINFALL ####
#########################

## loop through year and calculate moran's i summer rainfall for each year
results.temp <- NULL
moran.summer.rain <- NULL
for (i in sort(unique(climate.autocorr$year))){
  climate.autocorr.year <- climate.autocorr[climate.autocorr$year==i,]
  
  climate.dist <- as.matrix(dist(cbind(climate.autocorr.year$easting, climate.autocorr.year$northing)))
  
  climate.dist.inv <- 1/climate.dist
  diag(climate.dist.inv) <- 0
  climate.dist.inv[is.infinite(climate.dist.inv)] <- 0
  
  moran.i.result <- Moran.I(climate.autocorr.year$summer_rain, climate.dist)
  moran.i.result<-as.data.frame(moran.i.result)
  
  results.temp <- data.frame(moran.i.result$observed,i)
  colnames(results.temp)[1] <- "moran.i.summer.rain"
  colnames(results.temp)[2] <- "year"
  moran.summer.rain <- rbind(moran.summer.rain, results.temp)
  
}

moran.final <- merge(moran.final, moran.summer.rain, by.x="year", by.y="year")

#########################
#### AUTUMN RAINFALL ####
#########################

## loop through year and calculate moran's i autumn rainfall for each year
results.temp <- NULL
moran.autumn.rain <- NULL
for (i in sort(unique(climate.autocorr$year))){
  climate.autocorr.year <- climate.autocorr[climate.autocorr$year==i,]
  
  climate.dist <- as.matrix(dist(cbind(climate.autocorr.year$easting, climate.autocorr.year$northing)))
  
  climate.dist.inv <- 1/climate.dist
  diag(climate.dist.inv) <- 0
  climate.dist.inv[is.infinite(climate.dist.inv)] <- 0
  
  moran.i.result <- Moran.I(climate.autocorr.year$autumn_rain, climate.dist)
  moran.i.result<-as.data.frame(moran.i.result)
  
  results.temp <- data.frame(moran.i.result$observed,i)
  colnames(results.temp)[1] <- "moran.i.autumn.rain"
  colnames(results.temp)[2] <- "year"
  moran.autumn.rain <- rbind(moran.autumn.rain, results.temp)
  
}

moran.final <- merge(moran.final, moran.autumn.rain, by.x="year", by.y="year")

#########################
#### WINTER RAINFALL ####
#########################

## loop through year and calculate moran's i winter rainfall for each year
results.temp <- NULL
moran.winter.rain <- NULL
for (i in sort(unique(climate.autocorr$year))){
  climate.autocorr.year <- climate.autocorr[climate.autocorr$year==i,]
  
  climate.dist <- as.matrix(dist(cbind(climate.autocorr.year$easting, climate.autocorr.year$northing)))
  
  climate.dist.inv <- 1/climate.dist
  diag(climate.dist.inv) <- 0
  climate.dist.inv[is.infinite(climate.dist.inv)] <- 0 
  
  moran.i.result <- Moran.I(climate.autocorr.year$winter_rain, climate.dist)
  moran.i.result<-as.data.frame(moran.i.result)
  
  results.temp <- data.frame(moran.i.result$observed,i)
  colnames(results.temp)[1] <- "moran.i.winter.rain"
  colnames(results.temp)[2] <- "year"
  moran.winter.rain <- rbind(moran.winter.rain, results.temp)
  
}

moran.final <- merge(moran.final, moran.winter.rain, by.x="year", by.y="year")

write.csv(moran.final, file = "../Results/Butterfly_results/moran.i.final.data.csv", row.names=FALSE)

######## LINEAR MODELS ##########

lm_spring_temp <- lm(moran.i.spring.temp ~ year, data = moran.final)
summary(lm_spring_temp)

lm_summer_temp <- lm(moran.i.summer.temp ~ year, data = moran.final)
summary(lm_summer_temp)

lm_autumn_temp <- lm(moran.i.autumn.temp ~ year, data = moran.final)
summary(lm_autumn_temp)

lm_winter_temp <- lm(moran.i.winter.temp ~ year, data = moran.final)
summary(lm_winter_temp)

lm_spring_rain <- lm(moran.i.spring.rain ~ year, data = moran.final)
summary(lm_spring_rain)

lm_summer_rain <- lm(moran.i.summer.rain ~ year, data = moran.final)
summary(lm_summer_rain)

lm_autumn_rain <- lm(moran.i.autumn.rain ~ year, data = moran.final)
summary(lm_autumn_rain)

lm_winter_rain <- lm(moran.i.winter.rain ~ year, data = moran.final)
summary(lm_winter_rain)
###### NONE OF THE LINEAR MODELS ARE SIGNIFICANT ######

###### polynomial models of year against moran's i
pm_spring_temp <- lm(moran.i.spring.temp ~ I(year^2), data=moran.final)
summary(pm_spring_temp)

pm_summer_temp <- lm(moran.i.summer.temp ~ I(year^2), data=moran.final)
summary(pm_summer_temp)

pm_autumn_temp <- lm(moran.i.autumn.temp ~ I(year^2), data=moran.final)
summary(pm_autumn_temp)

pm_winter_temp <- lm(moran.i.winter.temp ~ I(year^2), data=moran.final)
summary(pm_winter_temp)

pm_spring_rain <- lm(moran.i.spring.rain ~ I(year^2), data=moran.final)
summary(pm_spring_rain)

pm_summer_rain <- lm(moran.i.summer.rain ~ I(year^2), data=moran.final)
summary(pm_summer_rain)

pm_autumn_rain <- lm(moran.i.autumn.rain ~ I(year^2), data=moran.final)
summary(pm_autumn_rain)

pm_winter_rain <- lm(moran.i.winter.rain ~ I(year^2), data=moran.final)
summary(pm_winter_rain)
##### NONE OF THE ABOVE ARE SIGNIFICANT

############################################
###### SPATIAL AUTOCORRELATION PLOTS  ######
############################################

moran.final <- read.csv("../Results/Butterfly_results/moran.i.final.data.csv", header=TRUE)

library(ggplot2)

## spring temperature ##
spring_temp <- ggplot(moran.final, aes(x=year, y=moran.i.spring.temp)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name="Moran's I spring temperature",breaks=seq(0.2,0.4,0.01)) +
  scale_x_continuous(name="Year",breaks=seq(1980,2012,2)) +
  theme_bw()
spring_temp

ggsave("../Graphs/Autcorrelation_plots/moran_i_spring_temp.png", plot=spring_temp, width=7, height=5)

## summer temperature ##
summer_temp <- ggplot(moran.final, aes(x=year, y=moran.i.summer.temp)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name="Moran's I summer temperature",breaks=seq(-0.3,-0.1,0.01)) +
  scale_x_continuous(name="Year",breaks=seq(1980,2012,2)) +
  theme_bw()
summer_temp

ggsave("../Graphs/Autcorrelation_plots/moran_i_summer_temp.png", plot=summer_temp, width=7, height=5)

## autumn temperature ##
autumn_temp <- ggplot(moran.final, aes(x=year, y=moran.i.autumn.temp)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name="Moran's I autumn temperature",breaks=seq(-0.3,-0.1,0.01)) +
  scale_x_continuous(name="Year",breaks=seq(1980,2012,2)) +
  theme_bw()
autumn_temp

ggsave("../Graphs/Autcorrelation_plots/moran_i_autumn_temp.png", plot=autumn_temp, width=7, height=5)

## winter temperature ##
winter_temp <- ggplot(moran.final, aes(x=year, y=moran.i.winter.temp)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name="Moran's I winter temperature",breaks=seq(-0.22,0,0.02)) +
  scale_x_continuous(name="Year",breaks=seq(1980,2012,2)) +
  theme_bw()
winter_temp

ggsave("../Graphs/Autcorrelation_plots/moran_i_winter_temp.png", plot=winter_temp, width=7, height=5)

## spring rainfall ##
spring_rain <- ggplot(moran.final, aes(x=year, y=moran.i.spring.rain)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name="Moran's I spring rainfall",breaks=seq(-0.2,0,0.01)) +
  scale_x_continuous(name="Year",breaks=seq(1980,2012,2)) +
  theme_bw()
spring_rain

ggsave("../Graphs/Autcorrelation_plots/moran_i_spring_rain.png", plot=spring_rain, width=7, height=5)

## summer rainfall ##
summer_rain <- ggplot(moran.final, aes(x=year, y=moran.i.summer.rain)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name="Moran's I summer rainfall",breaks=seq(-0.3,0,0.02)) +
  scale_x_continuous(name="Year",breaks=seq(1980,2012,2)) +
  theme_bw()
summer_rain 

ggsave("../Graphs/Autcorrelation_plots/moran_i_summer_rain.png", plot=summer_rain, width=7, height=5)

## autumn rainfall ##
autumn_rain <- ggplot(moran.final, aes(x=year, y=moran.i.autumn.rain)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name="Moran's I autumn rainfall",breaks=seq(-0.3,0,0.02)) +
  scale_x_continuous(name="Year",breaks=seq(1980,2012,2)) +
  theme_bw()
autumn_rain

ggsave("../Graphs/Autcorrelation_plots/moran_i_autumn_rain.png", plot=autumn_rain, width=7, height=5)

## winter rainfall ##
winter_rain <- ggplot(moran.final, aes(x=year, y=moran.i.winter.rain)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name="Moran's I winter rainfall",breaks=seq(-0.2,0,0.01)) +
  scale_x_continuous(name="Year",breaks=seq(1980,2012,2)) +
  theme_bw()
winter_rain

ggsave("../Graphs/Autcorrelation_plots/moran_i_winter_rain.png", plot=winter_rain, width=7, height=5)

################################################################
## CALCULATE MEAN SPATIAL AUTOCORRELATION USING MOVING WINDOW ##
################################################################

moran.final <- read.csv("../Results/Butterfly_results/moran.i.final.data.csv", header=TRUE)

## take mean moran's i across all 8 climate variables 
moran.final$moran.i.climate <- rowMeans(moran.final[,2:9])

## remove 8 climate varibles, leave year and moran.i.cliamte columns ## 
moran.final = subset(moran.final, select = -c(moran.i.spring.temp, moran.i.summer.temp, moran.i.autumn.temp, moran.i.winter.temp, moran.i.spring.rain, moran.i.summer.rain, moran.i.autumn.rain, moran.i.winter.rain))

## calculate mean moran's i for each moving winow ##
year.list <- 1980:2003
results.moran <- NULL

# loop through each year 
for (i in year.list){
    start.year <- i 
    mid.year <- i+4.5
    print(paste("mid.year=",mid.year))
    end.year <- i+9
    
    moran.i.10.year.data <- moran.final[moran.final$year>=start.year&moran.final$year<=end.year,]
  
    mean.moran.i <- mean(moran.i.10.year.data$moran.i.climate)
    
    results.temp <- data.frame(mid.year,mean.moran.i)
    results.moran <- rbind(results.moran,results.temp)
  
} # end i in year.list

## add 0.5 to mid year
results.moran$mid.year <- results.moran$mid.year+0.5

write.csv(results.moran, file = "../Results/Butterfly_results/moran.i.mean.results.csv", row.names=FALSE)


