###########################################################
## Title: Moving window climate synchrony calculation  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: March 2019
##########################################################

rm(list=ls()) # clear R
library(dplyr)
library(maptools)

#############################################################
###################### MEAN TEMPERATURE ##################### 
#############################################################

### add all files 
file_names <- list.files("../Data/MetOffice_data/Monthly_mean_temp_1980_2016", pattern=".txt") # where all the files are stored
final_table<-NULL
## read in ascii files and collate all data
for (f in file_names[]){
  print(f)
  temp_table<-NULL
    temp_table<- data.frame(readAsciiGrid(paste("../Data/MetOffice_data/Monthly_mean_temp_1980_2016",f,sep="/")))
    names(temp_table)[1]<-"mean_temp" ## change first column name to mean_temp
    temp_table$year<-substr(f,54,57) ## take characteres 54-57 from file name and paste into year
    temp_table$month<-substr(f,58,59) ## take characteres 58-59 from file name and paste into month
      final_table<-rbind(final_table,temp_table)
}

length(unique(final_table$year)) ## 37 years (1980-2016)
length(unique(final_table$month)) ## 12 months
## change colnames to easting and northing
names(final_table)[2:3] <- c("easting","northing")

head(final_table)
summary(final_table)

final_table$month <- as.numeric(final_table$month)
final_table$ag.month<-final_table$month
final_table$ag.month[final_table$month==1|final_table$month==2|final_table$month==3]<-"a"   # split months into 4 groups
final_table$ag.month[final_table$month==4|final_table$month==5|final_table$month==6]<-"b"
final_table$ag.month[final_table$month==7|final_table$month==8|final_table$month==9]<-"c"
final_table$ag.month[final_table$month==10|final_table$month==11|final_table$month==12]<-"d"
final_table$ag.month<-as.factor(final_table$ag.month)

agg.data<-aggregate(final_table$mean_temp, by = list(easting = final_table$easting, northing = final_table$northing, 
                year = final_table$year, ag = final_table$ag.month), mean)
names(agg.data)[5]<-"mean_temp"
head(agg.data)

## save file
write.csv(agg.data, file="../Data/MetOffice_data/Mean_temp_1980_2016_final.csv", row.names=FALSE)


########################################################
#################### MEAN RAINFALL ##################### 
########################################################

rm(list=ls()) # clear R

### add all files 
file_names <- list.files("../Data/MetOffice_data/Monthly_rainfall_1980_2016", pattern=".txt") # where all the files are stored
final_table<-NULL
## read in ascii files and collate all data
for (f in file_names[]){
  print(f)
  temp_table<-NULL
  temp_table<- data.frame(readAsciiGrid(paste("../Data/MetOffice_data/Monthly_rainfall_1980_2016",f,sep="/")))
  names(temp_table)[1]<-"mean_rainfall" ## change first column name to mean_temp
  temp_table$year<-substr(f,46,49) ## take characteres 54-57 from file name and paste into year
  temp_table$month<-substr(f,50,51) ## take characteres 58-59 from file name and paste into month
  final_table<-rbind(final_table,temp_table)
}

length(unique(final_table$year)) ## 37 years (1980-2016)
length(unique(final_table$month)) ## 12 months
## change colnames to easting and northing
names(final_table)[2:3] <- c("easting","northing")

head(final_table)
summary(final_table)

final_table$month <- as.numeric(final_table$month)
final_table$ag.month<-final_table$month
final_table$ag.month[final_table$month==1|final_table$month==2|final_table$month==3]<-"a"   # split months into 4 groups
final_table$ag.month[final_table$month==4|final_table$month==5|final_table$month==6]<-"b"
final_table$ag.month[final_table$month==7|final_table$month==8|final_table$month==9]<-"c"
final_table$ag.month[final_table$month==10|final_table$month==11|final_table$month==12]<-"d"
final_table$ag.month<-as.factor(final_table$ag.month)

agg.data<-aggregate(final_table$mean_rainfall, by = list(easting = final_table$easting, northing = final_table$northing, 
                                                     year = final_table$year, ag = final_table$ag.month), mean)
names(agg.data)[5]<-"mean_rainfall"
head(agg.data)

## save file
write.csv(agg.data, file="../Data/MetOffice_data/Mean_rainfall_1980_2016_final.csv", row.names=FALSE)


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

##### CALCULATE SEASONAL MEAN TEMPERATURE SYNCHRONY

## read in data
mean_temp <- read.csv("../Data/MetOffice_data/Mean_temp_1980_2016_final.csv", header=TRUE)
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE)

## find unique pair_attr sites to merge with mean_temp sites
site1 <- unique(subset(pair_attr[c(2,9,10)]))
site2 <- unique(subset(pair_attr[c(3,11,12)]))
colnames(site1) <- c("site", "east", "north")
colnames(site2) <- c("site", "east", "north")
site_list <- rbind(site1, site2)
site_list <- unique(site_list) ## 701 sites

## change coordinates from 1km to 5km
nrow(site_list)    # 2131 sites
for (i in 1:length(site_list$east)){
  ifelse(site_list$east[i]%%10000<5000,y<-2500,y<-7500)
  site_list$east.5k[i]<-(site_list$east[i]%/%10000)*10000+y
}
for (i in 1:length(site_list$north)){
  ifelse(site_list$north[i]%%10000<5000,y<-2500,y<-7500)
  site_list$north.5k[i]<-(site_list$north[i]%/%10000)*10000+y
}
site_list<-site_list[,c(1,4,5)]
site_list[site_list$site==3,]   # 517500,282500
nrow(site_list)      # 701 sites

## merge pair_attr site data with mean_temp (left with only UKBMS sites)
mean_temp <- merge(mean_temp, site_list, by.x=c("easting", "northing"), by.y=c("east.5k", "north.5k"))
length(unique(mean_temp$site)) # 686 sites

#### calculate syncrhony
final_pair_data <- NULL
season <- unique(mean_temp$ag)
### split based on season ###
for (g in season){ # loop through each season
  season_data <- mean_temp[mean_temp$ag==g,]
  total_comp <- NULL
  print(paste("season",g))    

year.list<-1980:2007

for (i in year.list){ # loop through year.list
  start.year<-i
  mid.year<-i+4.5
  print(paste("mid.year=",mid.year)) 
  end.year<-i+9
  subset_10yr_data<-season_data[season_data$year>=start.year&season_data$year<=end.year,]
  
  ################################
  # create a matrix to be filled #
  site.list <- unique(subset_10yr_data$site)
  year.list.temp <- min(subset_10yr_data$year):max(subset_10yr_data$year)
  Temp.matrix<-matrix(c(subset_10yr_data$mean_temp), nrow=length(year.list.temp))
  ncol(Temp.matrix)
  rownames(Temp.matrix)<-year.list.temp
  colnames(Temp.matrix)<-site.list
  
  length(site.list)
  
  num.ts <- length(site.list)   # number of time series
  TS <- matrix(Temp.matrix, ncol=num.ts)    # stores Grow data in a matrix     (simply removes names!)
  # create site match table used to correct sites names after calculating synchrony
  site_match <- data.frame(TS_name = 1:ncol(TS), site_name = colnames(Temp.matrix))
  
  ###############################
  # cross-correlation functions #
  pair.list <- t(combn(1:num.ts,2))
  nrow(pair.list)
  
  colnames(pair.list) <- c("site1","site2")
  
  ## calculate cross-correlation functions at different time lags
  ## (NB - we are only be interested in the current year, i.e. lag = 0)
  max.lag <- 0 # maximum time lag over which cross-correlations are assessed
  CCF <- data.frame(matrix(NA, ncol = 2*max.lag+1, nrow = nrow(pair.list)))	## matrix to store cross-correlation coefficients
  names(CCF) <- c(paste("lag", -max.lag:max.lag, sep="")) ## adds column names
  
  for(k in 1:nrow(pair.list)){
    ## calculation of CCF...
    try(CCF[k,ncol(CCF)] <- ccf(TS[,pair.list[k,1]], TS[,pair.list[k,2]], lag.max=max.lag, na.action=na.exclude, plot=F, type="correlation")$acf, silent=T) ## use try to prevent crashes if no data
  } ## this will take a long time!!! ## end k in nrow(pair.list)
  
  pair.attr <- pair.list ## matrix to hold attributes of pairs...
  colnames(pair.attr) <- c("site1","site2")
  pair.attr <- cbind(pair.attr, CCF)     ### add in correlation scores and number of comparisons at each site.
  head(pair.attr)
  
  ###  correct site numbers ### 
  pair.attr <- merge(pair.attr, site_match, by.x="site1", by.y="TS_name")  
  pair.attr <- pair.attr[,c("site_name", "site2", "lag0")]
  names(pair.attr) <- gsub("site_name", "site1", names(pair.attr))
  pair.attr <- merge(pair.attr, site_match, by.x="site2", by.y="TS_name")  
  pair.attr <- pair.attr[,c("site1", "site_name", "lag0")]
  names(pair.attr) <- gsub("site_name", "site2", names(pair.attr))
  
  pair.attr <- pair.attr[!is.na(pair.attr$lag0),]
  pair.attr$mid.year<-mid.year
  pair.attr$start.year<-start.year
  pair.attr$end.year<-end.year
  
  all_pair_attr <- NULL
  all_pair_attr <- rbind(all_pair_attr, pair.attr) 
  all_pair_attr$season <- g
  final_pair_data <- rbind(final_pair_data, all_pair_attr)
  
    } ## end in year.list
  } ## end in season

head(final_pair_data)

## save data
write.csv(final_pair_data, file="../Data/MetOffice_data/final_pair_data_mean_temp.csv", row.names=FALSE)

##### CALCULATE SEASONAL RAINFALL SYNCHRONY

rm(list=ls()) # clear R

## read in data
rainfall <- read.csv("../Data/MetOffice_data/Mean_rainfall_1980_2016_final.csv", header=TRUE)
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr_no_zeros2.csv", header=TRUE)

## find unique pair_attr sites to merge with mean_temp sites
site1 <- unique(subset(pair_attr[c(2,9,10)]))
site2 <- unique(subset(pair_attr[c(3,11,12)]))
colnames(site1) <- c("site", "east", "north")
colnames(site2) <- c("site", "east", "north")
site_list <- rbind(site1, site2)
site_list <- unique(site_list) ## 701 sites

## change coordinates from 1km to 5km
nrow(site_list)    # 2131 sites
for (i in 1:length(site_list$east)){
  ifelse(site_list$east[i]%%10000<5000,y<-2500,y<-7500)
  site_list$east.5k[i]<-(site_list$east[i]%/%10000)*10000+y
}
for (i in 1:length(site_list$north)){
  ifelse(site_list$north[i]%%10000<5000,y<-2500,y<-7500)
  site_list$north.5k[i]<-(site_list$north[i]%/%10000)*10000+y
}
site_list<-site_list[,c(1,4,5)]
site_list[site_list$site==3,]   # 517500,282500
nrow(site_list)      # 701 sites

## merge pair_attr site data with rainfall (left with only UKBMS sites)
rainfall <- merge(rainfall, site_list, by.x=c("easting", "northing"), by.y=c("east.5k", "north.5k"))
length(unique(rainfall$site)) # 686 sites

#### calculate syncrhony
final_pair_data <- NULL
season <- unique(rainfall$ag)
### split based on season ###
for (g in season){ # loop through each season
  season_data <- rainfall[rainfall$ag==g,]
  total_comp <- NULL
  print(paste("season",g))    
  
  year.list<-1980:2007
  
  for (i in year.list){ # loop through year.list
    start.year<-i
    mid.year<-i+4.5
    print(paste("mid.year=",mid.year)) 
    end.year<-i+9
    subset_10yr_data<-season_data[season_data$year>=start.year&season_data$year<=end.year,]
    
    ################################
    # create a matrix to be filled #
    site.list <- unique(subset_10yr_data$site)
    year.list.temp <- min(subset_10yr_data$year):max(subset_10yr_data$year)
    Rain.matrix<-matrix(c(subset_10yr_data$mean_rainfall), nrow=length(year.list.temp))
    ncol(Rain.matrix)
    rownames(Rain.matrix)<-year.list.temp
    colnames(Rain.matrix)<-site.list
    
    length(site.list)
    
    num.ts <- length(site.list)   # number of time series
    TS <- matrix(Rain.matrix, ncol=num.ts)    # stores Grow data in a matrix     (simply removes names!)
    # create site match table used to correct sites names after calculating synchrony
    site_match <- data.frame(TS_name = 1:ncol(TS), site_name = colnames(Rain.matrix))
    
    ###############################
    # cross-correlation functions #
    pair.list <- t(combn(1:num.ts,2))
    nrow(pair.list)
    
    colnames(pair.list) <- c("site1","site2")
    
    ## calculate cross-correlation functions at different time lags
    ## (NB - we are only be interested in the current year, i.e. lag = 0)
    max.lag <- 0 # maximum time lag over which cross-correlations are assessed
    CCF <- data.frame(matrix(NA, ncol = 2*max.lag+1, nrow = nrow(pair.list)))	## matrix to store cross-correlation coefficients
    names(CCF) <- c(paste("lag", -max.lag:max.lag, sep="")) ## adds column names
    
    for(k in 1:nrow(pair.list)){
      ## calculation of CCF...
      try(CCF[k,ncol(CCF)] <- ccf(TS[,pair.list[k,1]], TS[,pair.list[k,2]], lag.max=max.lag, na.action=na.exclude, plot=F, type="correlation")$acf, silent=T) ## use try to prevent crashes if no data
    } ## this will take a long time!!! ## end k in nrow(pair.list)
    
    pair.attr <- pair.list ## matrix to hold attributes of pairs...
    colnames(pair.attr) <- c("site1","site2")
    pair.attr <- cbind(pair.attr, CCF)     ### add in correlation scores and number of comparisons at each site.
    head(pair.attr)
    
    ###  correct site numbers ### 
    pair.attr <- merge(pair.attr, site_match, by.x="site1", by.y="TS_name")  
    pair.attr <- pair.attr[,c("site_name", "site2", "lag0")]
    names(pair.attr) <- gsub("site_name", "site1", names(pair.attr))
    pair.attr <- merge(pair.attr, site_match, by.x="site2", by.y="TS_name")  
    pair.attr <- pair.attr[,c("site1", "site_name", "lag0")]
    names(pair.attr) <- gsub("site_name", "site2", names(pair.attr))
    
    pair.attr <- pair.attr[!is.na(pair.attr$lag0),]
    pair.attr$mid.year<-mid.year
    pair.attr$start.year<-start.year
    pair.attr$end.year<-end.year
    
    all_pair_attr <- NULL
    all_pair_attr <- rbind(all_pair_attr, pair.attr) 
    all_pair_attr$season <- g
    final_pair_data <- rbind(final_pair_data, all_pair_attr)
    
  } ## end in year.list
} ## end in season

head(final_pair_data)

## save file
write.csv(final_pair_data, file="../Data/MetOffice_data/final_pair_data_mean_rainfall.csv", row.names=FALSE)




#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

## read in final pair data for temp and rainfall
final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp.csv", header=TRUE) ## need to re-run 
final_pair_data_rain <- read.csv("../Data/MetOffice_data/final_pair_data_mean_rainfall.csv", header=TRUE)

###### Rainfall ###### 

## create PairID column
final_pair_data_rain$pair.id <- paste("ID", final_pair_data_rain$site1, final_pair_data_rain$site2, sep = "_")

final_pair_data_rain$mid.year <- as.factor(final_pair_data_rain$mid.year)
final_pair_data_rain$pair.id <- as.character(final_pair_data_rain$pair.id)
final_pair_data_rain$season <- as.factor(final_pair_data_rain$season)

## split into 4 dataframes (one for each season)
winter_rainfall <- final_pair_data_rain[final_pair_data_rain$season=="a",]
winter_rainfall <- droplevels(winter_rainfall)
spring_rainfall <- final_pair_data_rain[final_pair_data_rain$season=="b",]
spring_rainfall <- droplevels(spring_rainfall)
summer_rainfall <- final_pair_data_rain[final_pair_data_rain$season=="c",]
summer_rainfall <- droplevels(summer_rainfall)
autumn_rainfall <- final_pair_data_rain[final_pair_data_rain$season=="d",]
autumn_rainfall <- droplevels(autumn_rainfall)

## run mixed effects model to extract coefficients for each dataframe
library(lme4)

###### WINTER RAINFALL ######
winter_rainfall_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=winter_rainfall)
summary(winter_rainfall_model)
## save results
results_table_winter_rain <- data.frame(summary(winter_rainfall_model)$coefficients[,1:3])
## change names and add in parameter column ##
names(results_table_winter_rain) <- c("synchrony", "SD", "t")
results_table_winter_rain$parameter <- paste(row.names(results_table_winter_rain))
rownames(results_table_winter_rain) <- 1:nrow(results_table_winter_rain)
## change parameter names to year
results_table_winter_rain$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
results_table_winter_rain$rescaled_sync <- results_table_winter_rain$synchrony*(100/results_table_winter_rain$synchrony[1])
results_table_winter_rain$rescaled_sd <- results_table_winter_rain$SD*(100/results_table_winter_rain$synchrony[1])
results_table_winter_rain$rescaled_ci <- results_table_winter_rain$rescaled_sd*1.96

## save final results table ##
write.csv(results_table_winter_rain, file = "../Results/Climate_results/winter_rainfall_synchrony.csv", row.names=FALSE)

###### SPRING RAINFALL ######
spring_rainfall_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=spring_rainfall)
summary(spring_rainfall_model)
## save results
results_table_spring_rain <- data.frame(summary(spring_rainfall_model)$coefficients[,1:3])
## change names and add in parameter column ##
names(results_table_spring_rain) <- c("synchrony", "SD", "t")
results_table_spring_rain$parameter <- paste(row.names(results_table_spring_rain))
rownames(results_table_spring_rain) <- 1:nrow(results_table_spring_rain)
## change parameter names to year
results_table_spring_rain$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
results_table_spring_rain$rescaled_sync <- results_table_spring_rain$synchrony*(100/results_table_spring_rain$synchrony[1])
results_table_spring_rain$rescaled_sd <- results_table_spring_rain$SD*(100/results_table_spring_rain$synchrony[1])
results_table_spring_rain$rescaled_ci <- results_table_spring_rain$rescaled_sd*1.96

## save final results table ##
write.csv(results_table_spring_rain, file = "../Results/Climate_results/spring_rainfall_synchrony.csv", row.names=FALSE)

###### SUMMER RAINFALL ######
summer_rainfall_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=summer_rainfall)
summary(summer_rainfall_model)
## save results
results_table_summer_rain <- data.frame(summary(summer_rainfall_model)$coefficients[,1:3])
## change names and add in parameter column ##
names(results_table_summer_rain) <- c("synchrony", "SD", "t")
results_table_summer_rain$parameter <- paste(row.names(results_table_summer_rain))
rownames(results_table_summer_rain) <- 1:nrow(results_table_summer_rain)
## change parameter names to year
results_table_summer_rain$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
results_table_summer_rain$rescaled_sync <- results_table_summer_rain$synchrony*(100/results_table_summer_rain$synchrony[1])
results_table_summer_rain$rescaled_sd <- results_table_summer_rain$SD*(100/results_table_summer_rain$synchrony[1])
results_table_summer_rain$rescaled_ci <- results_table_summer_rain$rescaled_sd*1.96

## save final results table ##
write.csv(results_table_summer_rain, file = "../Results/Climate_results/summer_rainfall_synchrony.csv", row.names=FALSE)

###### AUTUMN RAINFALL ######
autumn_rainfall_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=autumn_rainfall)
summary(autumn_rainfall_model)
## save results
results_table_autumn_rain <- data.frame(summary(autumn_rainfall_model)$coefficients[,1:3])
## change names and add in parameter column ##
names(results_table_autumn_rain) <- c("synchrony", "SD", "t")
results_table_autumn_rain$parameter <- paste(row.names(results_table_autumn_rain))
rownames(results_table_autumn_rain) <- 1:nrow(results_table_autumn_rain)
## change parameter names to year
results_table_autumn_rain$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
results_table_autumn_rain$rescaled_sync <- results_table_autumn_rain$synchrony*(100/results_table_autumn_rain$synchrony[1])
results_table_autumn_rain$rescaled_sd <- results_table_autumn_rain$SD*(100/results_table_autumn_rain$synchrony[1])
results_table_autumn_rain$rescaled_ci <- results_table_autumn_rain$rescaled_sd*1.96

## save final results table ##
write.csv(results_table_autumn_rain, file = "../Results/Climate_results/autumn_rainfall_synchrony.csv", row.names=FALSE)

par(mfrow=c(2,2))
plot(rescaled_sync ~ parameter, data=results_table_winter_rain)
plot(rescaled_sync ~ parameter, data=results_table_spring_rain)
plot(rescaled_sync ~ parameter, data=results_table_summer_rain)
plot(rescaled_sync ~ parameter, data=results_table_autumn_rain)

