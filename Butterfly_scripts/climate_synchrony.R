###########################################################
## Title: Moving window climate synchrony calculation  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: March 2019
##########################################################

rm(list=ls()) # clear R
library(dplyr)
library(maptools)
library(ggplot2)
options(scipen=999)

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
unique(final_table$year)
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

unique(agg.data$year)

## save file
write.csv(agg.data, file="../Data/MetOffice_data/Mean_temp_1980_2016_final.csv", row.names=FALSE)


########################################################
#################### MEAN RAINFALL ##################### 
########################################################

rm(list=ls()) # clear R

### add all files 
file_names <- list.files("../Data/MetOffice_data/Monthly_rainfall_1980_2016", pattern=".txt") # where all the files are stored
final_table2<-NULL
## read in ascii files and collate all data
for (f in file_names[]){
  print(f)
  temp_table<-NULL
  temp_table<- data.frame(readAsciiGrid(paste("../Data/MetOffice_data/Monthly_rainfall_1980_2016",f,sep="/")))
  names(temp_table)[1]<-"mean_rainfall" ## change first column name to mean_temp
  temp_table$year<-substr(f,46,49) ## take characteres 54-57 from file name and paste into year
  temp_table$month<-substr(f,50,51) ## take characteres 58-59 from file name and paste into month
  final_table2<-rbind(final_table2,temp_table)
}

length(unique(final_table2$year)) ## 37 years (1980-2016)
length(unique(final_table2$month)) ## 12 months
## change colnames to easting and northing
names(final_table2)[2:3] <- c("easting","northing")

head(final_table2)
summary(final_table2)

final_table2$month <- as.numeric(final_table2$month)
final_table2$ag.month<-final_table2$month
final_table2$ag.month[final_table2$month==1|final_table2$month==2|final_table2$month==3]<-"a"   # split months into 4 groups
final_table2$ag.month[final_table2$month==4|final_table2$month==5|final_table2$month==6]<-"b"
final_table2$ag.month[final_table2$month==7|final_table2$month==8|final_table2$month==9]<-"c"
final_table2$ag.month[final_table2$month==10|final_table2$month==11|final_table2$month==12]<-"d"
final_table2$ag.month<-as.factor(final_table2$ag.month)

agg.data2<-aggregate(final_table2$mean_rainfall, by = list(easting = final_table2$easting, northing = final_table2$northing, 
                                                     year = final_table2$year, ag = final_table2$ag.month), mean)
names(agg.data2)[5]<-"mean_rainfall"
head(agg.data2)

## save file
write.csv(agg.data2, file="../Data/MetOffice_data/Mean_rainfall_1980_2016_final.csv", row.names=FALSE)

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

rm(list=ls()) # clear R

## read in data
mean_temp <- read.csv("../Data/MetOffice_data/Mean_temp_1980_2016_final.csv", header=TRUE)
rainfall <- read.csv("../Data/MetOffice_data/Mean_rainfall_1980_2016_final.csv", header=TRUE)
pair_attr <- read.csv("../Data/Butterfly_sync_data/pair_attr.csv", header=TRUE)

## find unique UKBMS pair_attr sites to merge with mean_temp sites
site1 <- unique(subset(pair_attr[c(2,9,10)]))
site2 <- unique(subset(pair_attr[c(3,11,12)]))
colnames(site1) <- c("site", "east", "north")
colnames(site2) <- c("site", "east", "north")
site_list <- rbind(site1, site2)
site_list <- unique(site_list) ## 701 sites

## change UKBMS coordinates from 1km to 5km
nrow(site_list)    # 701
for (i in 1:length(site_list$east)){
  ifelse(site_list$east[i]%%10000<5000,y<-2500,y<-7500)
  site_list$east.5k[i]<-(site_list$east[i]%/%10000)*10000+y
}
for (i in 1:length(site_list$north)){
  ifelse(site_list$north[i]%%10000<5000,y<-2500,y<-7500)
  site_list$north.5k[i]<-(site_list$north[i]%/%10000)*10000+y
}
site_list<-site_list[,c(1,4,5)] ## 701 sites
## remove duplicates where two sites have the same 5km easting and northing
site_list <- site_list[!duplicated(t(apply(site_list[2:3], 1, sort))),] ## 499 sites

## merge pair_attr site data with mean_temp (left with only UKBMS sites)
mean_temp <- merge(mean_temp, site_list, by.x=c("easting", "northing"), by.y=c("east.5k", "north.5k"))
length(unique(mean_temp$site)) # 488 sites which match with climate and UKBMS data
## merge pair_attr site data with rainfall (left with only UKBMS sites)
rainfall <- merge(rainfall, site_list, by.x=c("easting", "northing"), by.y=c("east.5k", "north.5k"))
length(unique(rainfall$site)) # 486 sites which match with climate and UKBMS data

##### CALCULATE SEASONAL MEAN TEMPERATURE SYNCHRONY

final_pair_data <- NULL
season <- unique(mean_temp$ag)
### split based on season ###
for (g in season[3]){ # loop through each season
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
  subset_10yr_data <- droplevels(subset_10yr_data)
  
  ################################
  # create a matrix to be filled #
  site.list <- unique(subset_10yr_data$site)
  year.list.temp <- min(subset_10yr_data$year):max(subset_10yr_data$year)
  # Temp.matrix<-matrix(c(subset_10yr_data$mean_temp), nrow=length(year.list.temp))
  # ncol(Temp.matrix)
  # rownames(Temp.matrix)<-year.list.temp
  # colnames(Temp.matrix)<-site.list

  Temp.matrix <- reshape2::acast(subset_10yr_data, year ~ site, value.var="mean_temp")

  num.ts <- length(site.list)   # number of time series
  TS <- matrix(Temp.matrix, ncol=num.ts)    # stores Grow data in a matrix     (simply removes names!)
  # create site match table used to correct sites names after calculating synchrony
  site_match <- data.frame(TS_name = 1:ncol(TS), site_name = colnames(Temp.matrix))
  
  ###############################
  # cross-correlation functions #
  pair.list <- t(combn(1:num.ts,2))
  nrow(pair.list)
  
  colnames(pair.list) <- c("site1","site2")

  ###  correct site numbers ###
  pair.list <- merge(pair.list, site_match, by.x="site1", by.y="TS_name")
  pair.list <- pair.list[,c("site_name", "site2")] ## get rid of site1 (which are not real site names)
  names(pair.list) <- gsub("site_name", "site1", names(pair.list))
  pair.list <- merge(pair.list, site_match, by.x="site2", by.y="TS_name")
  pair.list <- pair.list[,c("site1", "site_name")] ## get rid of site2 (which are not real site names)
  names(pair.list) <- gsub("site_name", "site2", names(pair.list)) ## now site1 and site2 are REAL site names
  
  ## merge pair.list with site data
  pair_attr2 <- unique(pair_attr[,c(2,3,7)]) ## site1, site2 and start year
  pair_attr2[] <- lapply(pair_attr2, factor)
  pair_attr2 <- pair_attr2[pair_attr2$start.year==i,] ## subset by start year i
  site_data_pair1 <- merge(pair.list, pair_attr2, by=c("site1", "site2"), all=FALSE)
  pair_attr_reverse <- pair_attr2
  names(pair_attr_reverse)[1:3] <- c("site2", "site1", "start.year")
  site_data_pair2 <- merge(pair.list, pair_attr_reverse, by=c("site1", "site2"), all=FALSE)
  site_data_pair <- rbind(site_data_pair1, site_data_pair2)
  site_data_pair <- unique(site_data_pair)
  site_data_pair <- site_data_pair[,c(1,2)]
  
  ## merge pair.list with site_data pair so pair.list only contains sites which are in the correct 10 year moving window
  #pair.list <- site_data_pair
  pair.list<-merge(pair.list,site_data_pair,by.x=c("site1", "site2"),by.y=c("site1", "site2")) ## pair list becomes same length as site_data_pair
  
  # change site names back to numbers
  pair.list <- merge(pair.list, site_match, by.x="site1", by.y="site_name")
  pair.list <- pair.list[,c("TS_name", "site2")]
  names(pair.list) <- gsub("TS_name", "site1", names(pair.list))
  pair.list <- merge(pair.list, site_match, by.x="site2", by.y="site_name")
  pair.list <- pair.list[,c("site1", "TS_name")]
  names(pair.list) <- gsub("TS_name", "site2", names(pair.list))
  
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

site1 <- unique(final_pair_data[,1, drop=FALSE])
site2 <- unique(final_pair_data[,2, drop=FALSE])
colnames(site1)[1] <- "site"
colnames(site2)[1] <- "site"
site_list <- rbind(site1, site2)
site_list <- unique(site_list) ## 487 sites

## save data
write.csv(final_pair_data, file="../Data/MetOffice_data/final_pair_data_mean_temp.csv", row.names=FALSE)

##### CALCULATE SEASONAL RAINFALL SYNCHRONY

final_pair_data2 <- NULL
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
    # Rain.matrix<-matrix(c(subset_10yr_data$mean_rainfall), nrow=length(year.list.temp))
    # ncol(Rain.matrix)
    # rownames(Rain.matrix)<-year.list.temp
    # colnames(Rain.matrix)<-site.list
    
    Rain.matrix <- reshape2::acast(subset_10yr_data, year ~ site, value.var="mean_rainfall")
    
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
    
    ###  correct site numbers ###
    pair.list <- merge(pair.list, site_match, by.x="site1", by.y="TS_name")
    pair.list <- pair.list[,c("site_name", "site2")] ## get rid of site1 (which are not real site names)
    names(pair.list) <- gsub("site_name", "site1", names(pair.list))
    pair.list <- merge(pair.list, site_match, by.x="site2", by.y="TS_name")
    pair.list <- pair.list[,c("site1", "site_name")] ## get rid of site2 (which are not real site names)
    names(pair.list) <- gsub("site_name", "site2", names(pair.list)) ## now site1 and site2 are REAL site names
    
    ## merge pair.list with site data
    pair_attr2 <- unique(pair_attr[,c(2,3,7)]) ## site1, site2 and start year
    pair_attr2[] <- lapply(pair_attr2, factor)
    pair_attr2 <- pair_attr2[pair_attr2$start.year==i,] ## subset by start year i
    site_data_pair1 <- merge(pair.list, pair_attr2, by=c("site1", "site2"), all=FALSE)
    pair_attr_reverse <- pair_attr2
    names(pair_attr_reverse)[1:3] <- c("site2", "site1", "start.year")
    site_data_pair2 <- merge(pair.list, pair_attr_reverse, by=c("site1", "site2"), all=FALSE)
    site_data_pair <- rbind(site_data_pair1, site_data_pair2)
    site_data_pair <- unique(site_data_pair)
    site_data_pair <- site_data_pair[,c(1,2)]
    
    ## merge pair.list with site_data pair so pair.list only contains sites which are in the correct 10 year moving window
    pair.list<-merge(pair.list,site_data_pair,by.x=c("site1", "site2"),by.y=c("site1", "site2")) ## pair list becomes same length as site_data_pair
    
    # change site names back to numbers
    pair.list <- merge(pair.list, site_match, by.x="site1", by.y="site_name")
    pair.list <- pair.list[,c("TS_name", "site2")]
    names(pair.list) <- gsub("TS_name", "site1", names(pair.list))
    pair.list <- merge(pair.list, site_match, by.x="site2", by.y="site_name")
    pair.list <- pair.list[,c("site1", "TS_name")]
    names(pair.list) <- gsub("TS_name", "site2", names(pair.list))
    
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
    final_pair_data2 <- rbind(final_pair_data2, all_pair_attr)
    
  } ## end in year.list
} ## end in season

head(final_pair_data2)

## save file
write.csv(final_pair_data2, file="../Data/MetOffice_data/final_pair_data_mean_rainfall.csv", row.names=FALSE)

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

rm(list=ls()) # clear R

## read in final pair data for temp and rainfall
final_pair_data_rain <- read.csv("../Data/MetOffice_data/final_pair_data_mean_rainfall.csv", header=TRUE)
final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp.csv", header=TRUE)  

## change values that are >0 to just above 0
## otherwise logit transformation doesn't work 
final_pair_data_rain$lag0[final_pair_data_rain$lag0<0 ] <- 0.000001

summ_data_temp <- final_pair_data_temp %>% 
  group_by(mid.year,season) %>% 
  summarise(freq = dplyr::n())

summ_data_rain <- final_pair_data_rain %>% 
  group_by(mid.year,season) %>% 
  summarise(freq = n()) 

######################
###### Rainfall ######
######################

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

logitTransform <- function(p) { log(p/(1-p)) }
winter_rainfall$lag0_logit <- logitTransform(winter_rainfall$lag0)

winter_rainfall_model <- lmer(log(lag0/(1-lag0)) ~ mid.year + (1|pair.id)-1, data=winter_rainfall)
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
## add 10 to each value so scale remains the same, but values are all positive
results_table_winter_rain$new_sync <- results_table_winter_rain$synchrony + 10 
## scale using new_sync
results_table_winter_rain$rescaled_sync <- results_table_winter_rain$new_sync*(100/results_table_winter_rain$new_sync[1]) ## rescale to 100
results_table_winter_rain$rescaled_sd <- results_table_winter_rain$SD*(100/results_table_winter_rain$new_sync[1])
results_table_winter_rain$rescaled_ci <- results_table_winter_rain$rescaled_sd*1.96
## remove new_sync column
results_table_winter_rain <- subset(results_table_winter_rain, select = -c(new_sync))
## save final results table ##
write.csv(results_table_winter_rain, file = "../Results/Climate_results/winter_rainfall_synchrony.csv", row.names=FALSE)

###### SPRING RAINFALL ######
spring_rainfall_model <- lmer(log(lag0/(1-lag0)) ~ mid.year + (1|pair.id)-1, data=spring_rainfall)
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
## add 10 to each value so scale remains the same, but values are all positive
results_table_spring_rain$new_sync <- results_table_spring_rain$synchrony + 10 
## scale using new_sync
results_table_spring_rain$rescaled_sync <- results_table_spring_rain$new_sync*(100/results_table_spring_rain$new_sync[1]) ## rescale to 100
results_table_spring_rain$rescaled_sd <- results_table_spring_rain$SD*(100/results_table_spring_rain$new_sync[1])
results_table_spring_rain$rescaled_ci <- results_table_spring_rain$rescaled_sd*1.96
## remove new_sync column
results_table_spring_rain <- subset(results_table_spring_rain, select = -c(new_sync))
## save final results table ##
write.csv(results_table_spring_rain, file = "../Results/Climate_results/spring_rainfall_synchrony.csv", row.names=FALSE)

###### SUMMER RAINFALL ######
summer_rainfall_model <- lmer(log(lag0/(1-lag0)) ~ mid.year + (1|pair.id)-1, data=summer_rainfall)
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
## add 10 to each value so scale remains the same, but values are all positive
results_table_summer_rain$new_sync <- results_table_summer_rain$synchrony + 10 
## scale using new_sync
results_table_summer_rain$rescaled_sync <- results_table_summer_rain$new_sync*(100/results_table_summer_rain$new_sync[1]) ## rescale to 100
results_table_summer_rain$rescaled_sd <- results_table_summer_rain$SD*(100/results_table_summer_rain$new_sync[1])
results_table_summer_rain$rescaled_ci <- results_table_summer_rain$rescaled_sd*1.96
## remove new_sync column
results_table_summer_rain <- subset(results_table_summer_rain, select = -c(new_sync))
## save final results table ##
write.csv(results_table_summer_rain, file = "../Results/Climate_results/summer_rainfall_synchrony.csv", row.names=FALSE)

###### AUTUMN RAINFALL ######
autumn_rainfall_model <- lmer(log(lag0/(1-lag0)) ~ mid.year + (1|pair.id)-1, data=autumn_rainfall)
summary(autumn_rainfall_model)
plot(autumn_rainfall_model)
qqnorm(residuals(autumn_rainfall_model))
## save results
results_table_autumn_rain <- data.frame(summary(autumn_rainfall_model)$coefficients[,1:3])
## change names and add in parameter column ##
names(results_table_autumn_rain) <- c("synchrony", "SD", "t")
results_table_autumn_rain$parameter <- paste(row.names(results_table_autumn_rain))
rownames(results_table_autumn_rain) <- 1:nrow(results_table_autumn_rain)
## change parameter names to year
results_table_autumn_rain$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
## add 10 to each value so scale remains the same, but values are all positive
results_table_autumn_rain$new_sync <- results_table_autumn_rain$synchrony + 10 
## scale using new_sync
results_table_autumn_rain$rescaled_sync <- results_table_autumn_rain$new_sync*(100/results_table_autumn_rain$new_sync[1]) ## rescale to 100
results_table_autumn_rain$rescaled_sd <- results_table_autumn_rain$SD*(100/results_table_autumn_rain$new_sync[1])
results_table_autumn_rain$rescaled_ci <- results_table_autumn_rain$rescaled_sd*1.96
## remove new_sync column
results_table_autumn_rain <- subset(results_table_autumn_rain, select = -c(new_sync))
## save final results table ##
write.csv(results_table_autumn_rain, file = "../Results/Climate_results/autumn_rainfall_synchrony.csv", row.names=FALSE)

##### plots ######
## load data
winter_rain <- read.csv("../Results/Climate_results/winter_rainfall_synchrony.csv", header=TRUE)
spring_rain <- read.csv("../Results/Climate_results/spring_rainfall_synchrony.csv", header=TRUE)
summer_rain <- read.csv("../Results/Climate_results/summer_rainfall_synchrony.csv", header=TRUE)
autumn_rain <- read.csv("../Results/Climate_results/autumn_rainfall_synchrony.csv", header=TRUE)

winter_rain_plot <- ggplot(winter_rain, aes(x = parameter, y = rescaled_sync)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Rainfall synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
winter_rain_plot

spring_rain_plot <- ggplot(spring_rain, aes(x = parameter, y = rescaled_sync)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Rainfall synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
spring_rain_plot

summer_rain_plot <- ggplot(summer_rain, aes(x = parameter, y = rescaled_sync)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Rainfall synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
summer_rain_plot

autumn_rain_plot <- ggplot(autumn_rain, aes(x = parameter, y = rescaled_sync)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Rainfall synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
autumn_rain_plot

## save all 4 rainfall plots
library(ggpubr)
rain_plots <- ggarrange(winter_rain_plot, spring_rain_plot, summer_rain_plot, autumn_rain_plot,
                        hjust = 0, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)"))
rain_plots
ggsave("../Graphs/Climate/seasonal_rainfall_synchrony.png", plot = rain_plots, width=12, height=10)

#########################
###### Temperature ######
#########################

## create PairID column
final_pair_data_temp$pair.id <- paste("ID", final_pair_data_temp$site1, final_pair_data_temp$site2, sep = "_")

final_pair_data_temp$mid.year <- as.factor(final_pair_data_temp$mid.year)
final_pair_data_temp$pair.id <- as.character(final_pair_data_temp$pair.id)
final_pair_data_temp$season <- as.factor(final_pair_data_temp$season)

## split into 4 dataframes (one for each season)
winter_temp <- final_pair_data_temp[final_pair_data_temp$season=="a",]
winter_temp <- droplevels(winter_temp)
spring_temp <- final_pair_data_temp[final_pair_data_temp$season=="b",]
spring_temp <- droplevels(spring_temp)
summer_temp <- final_pair_data_temp[final_pair_data_temp$season=="c",]
summer_temp <- droplevels(summer_temp)
autumn_temp <- final_pair_data_temp[final_pair_data_temp$season=="d",]
autumn_temp <- droplevels(autumn_temp)

## run mixed effects model to extract coefficients for each dataframe
library(lme4)

###### WINTER TEMPERATURE ######
winter_temp_model <- lmer(log(lag0/(1-lag0)) ~ mid.year + (1|pair.id)-1, data=winter_temp)
summary(winter_temp_model)
## save results
results_table_winter_temp <- data.frame(summary(winter_temp_model)$coefficients[,1:3])
## change names and add in parameter column ##
names(results_table_winter_temp) <- c("synchrony", "SD", "t")
results_table_winter_temp$parameter <- paste(row.names(results_table_winter_temp))
rownames(results_table_winter_temp) <- 1:nrow(results_table_winter_temp)
## change parameter names to year
results_table_winter_temp$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
## add 10 to each value so scale remains the same, but values are all positive
results_table_winter_temp$new_sync <- results_table_winter_temp$synchrony + 10 
## scale using new_sync
results_table_winter_temp$rescaled_sync <- results_table_winter_temp$new_sync*(100/results_table_winter_temp$new_sync[1]) ## rescale to 100
results_table_winter_temp$rescaled_sd <- results_table_winter_temp$SD*(100/results_table_winter_temp$new_sync[1])
results_table_winter_temp$rescaled_ci <- results_table_winter_temp$rescaled_sd*1.96
## remove new_sync column
results_table_winter_temp <- subset(results_table_winter_temp, select = -c(new_sync))
## save final results table ##
write.csv(results_table_winter_temp, file = "../Results/Climate_results/winter_temp_synchrony.csv", row.names=FALSE)

###### SPRING TEMPERATURE ######
spring_temp_model <- lmer(log(lag0/(1-lag0)) ~ mid.year + (1|pair.id)-1, data=spring_temp)
summary(spring_temp_model)
## save results
results_table_spring_temp <- data.frame(summary(spring_temp_model)$coefficients[,1:3])
## change names and add in parameter column ##
names(results_table_spring_temp) <- c("synchrony", "SD", "t")
results_table_spring_temp$parameter <- paste(row.names(results_table_spring_temp))
rownames(results_table_spring_temp) <- 1:nrow(results_table_spring_temp)
## change parameter names to year
results_table_spring_temp$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
## add 10 to each value so scale remains the same, but values are all positive
results_table_spring_temp$new_sync <- results_table_spring_temp$synchrony + 10 
## scale using new_sync
results_table_spring_temp$rescaled_sync <- results_table_spring_temp$new_sync*(100/results_table_spring_temp$new_sync[1]) ## rescale to 100
results_table_spring_temp$rescaled_sd <- results_table_spring_temp$SD*(100/results_table_spring_temp$new_sync[1])
results_table_spring_temp$rescaled_ci <- results_table_spring_temp$rescaled_sd*1.96
## remove new_sync column
results_table_spring_temp <- subset(results_table_spring_temp, select = -c(new_sync))
## save final results table ##
write.csv(results_table_spring_temp, file = "../Results/Climate_results/spring_temp_synchrony.csv", row.names=FALSE)

###### SUMMER TEMPERATURE ######
# logitTransform <- function(p) { log(p/(1-p)) }
# summer_temp$lag0_logit <- logitTransform(summer_temp$lag0)
# 
# hist(summer_temp$lag0)
# hist(summer_temp$lag0_logit)
# 
# summer_temp_model1 <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=summer_temp)
summer_temp_model <- lmer(log(lag0/(1-lag0)) ~ mid.year + (1|pair.id)-1, data=summer_temp)

# qqnorm(residuals(summer_temp_model1))
# qqline(residuals(summer_temp_model1))
# plot(summer_temp_model1)
# 
# qqnorm(residuals(summer_temp_model))
# qqline(residuals(summer_temp_model))
# plot(summer_temp_model)


summary(summer_temp_model)
## save results
results_table_summer_temp <- data.frame(summary(summer_temp_model)$coefficients[,1:3])
## change names and add in parameter column ##
names(results_table_summer_temp) <- c("synchrony", "SD", "t")
results_table_summer_temp$parameter <- paste(row.names(results_table_summer_temp))
rownames(results_table_summer_temp) <- 1:nrow(results_table_summer_temp)
## change parameter names to year
results_table_summer_temp$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
## add 10 to each value so scale remains the same, but values are all positive
results_table_summer_temp$new_sync <- results_table_summer_temp$synchrony + 10 
## scale using new_sync
results_table_summer_temp$rescaled_sync <- results_table_summer_temp$new_sync*(100/results_table_summer_temp$new_sync[1]) ## rescale to 100
results_table_summer_temp$rescaled_sd <- results_table_summer_temp$SD*(100/results_table_summer_temp$new_sync[1])
results_table_summer_temp$rescaled_ci <- results_table_summer_temp$rescaled_sd*1.96
## remove new_sync column
results_table_summer_temp <- subset(results_table_summer_temp, select = -c(new_sync))
## save final results table ##
write.csv(results_table_summer_temp, file = "../Results/Climate_results/summer_temp_synchrony.csv", row.names=FALSE)

###### AUTUMN TEMPERATURE ######
autumn_temp_model <- lmer(log(lag0/(1-lag0)) ~ mid.year + (1|pair.id)-1, data=autumn_temp)
summary(autumn_temp_model)
## save results
results_table_autumn_temp <- data.frame(summary(autumn_temp_model)$coefficients[,1:3])
## change names and add in parameter column ##
names(results_table_autumn_temp) <- c("synchrony", "SD", "t")
results_table_autumn_temp$parameter <- paste(row.names(results_table_autumn_temp))
rownames(results_table_autumn_temp) <- 1:nrow(results_table_autumn_temp)
## change parameter names to year
results_table_autumn_temp$parameter <- rep(1985:2012)

### rescale estimate, SD and CI ### 
## add 10 to each value so scale remains the same, but values are all positive
results_table_autumn_temp$new_sync <- results_table_autumn_temp$synchrony + 10 
## scale using new_sync
results_table_autumn_temp$rescaled_sync <- results_table_autumn_temp$new_sync*(100/results_table_autumn_temp$new_sync[1]) ## rescale to 100
results_table_autumn_temp$rescaled_sd <- results_table_autumn_temp$SD*(100/results_table_autumn_temp$new_sync[1])
results_table_autumn_temp$rescaled_ci <- results_table_autumn_temp$rescaled_sd*1.96
## remove new_sync column
results_table_autumn_temp <- subset(results_table_autumn_temp, select = -c(new_sync))
## save final results table ##
write.csv(results_table_autumn_temp, file = "../Results/Climate_results/autumn_temp_synchrony.csv", row.names=FALSE)

##### plots ######
## load data

winter_temp <- read.csv("../Results/Climate_results/winter_temp_synchrony.csv", header=TRUE)
spring_temp <- read.csv("../Results/Climate_results/spring_temp_synchrony.csv", header=TRUE)
summer_temp <- read.csv("../Results/Climate_results/summer_temp_synchrony.csv", header=TRUE)
autumn_temp <- read.csv("../Results/Climate_results/autumn_temp_synchrony.csv", header=TRUE)

winter_temp_plot <- ggplot(winter_temp, aes(x = parameter, y = rescaled_sync)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Temperature synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
winter_temp_plot

spring_temp_plot <- ggplot(spring_temp, aes(x = parameter, y = rescaled_sync)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Temperature synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
spring_temp_plot

summer_temp_plot <- ggplot(summer_temp, aes(x = parameter, y = rescaled_sync)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Temperature synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
summer_temp_plot

autumn_temp_plot <- ggplot(autumn_temp, aes(x = parameter, y = rescaled_sync)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
  geom_point(size=2) + 
  labs(x = "Mid-year of moving window", y = "Temperature synchrony") +
  #scale_y_continuous(breaks=seq(40,160,10)) +
  scale_x_continuous(breaks=seq(1985,2012,3)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=3) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
autumn_temp_plot

## save all 4 temperature plots
library(ggpubr)
temp_plots <- ggarrange(winter_temp_plot, spring_temp_plot, summer_temp_plot, autumn_temp_plot,
          hjust = 0, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)"))
temp_plots
ggsave("../Graphs/Climate/seasonal_temperature_synchrony.png", plot = temp_plots, width=12, height=10)



##############################################################################################################
##############################################################################################################
##############################################################################################################

rm(list=ls()) # clear R
library(lme4)
library(lmerTest)

## read in data
final_pair_data_rain <- read.csv("../Data/MetOffice_data/final_pair_data_mean_rainfall.csv", header=TRUE)
final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp.csv", header=TRUE)  

summ_data_temp <- final_pair_data_temp %>% 
  group_by(mid.year,season) %>% 
  summarise(freq = n())

summ_data_rain <- final_pair_data_rain %>% 
  group_by(mid.year,season) %>% 
  summarise(freq = n()) 

## create PairID column
final_pair_data_rain$pair.id <- paste("ID", final_pair_data_rain$site1, final_pair_data_rain$site2, sep = "_")
final_pair_data_temp$pair.id <- paste("ID", final_pair_data_temp$site1, final_pair_data_temp$site2, sep = "_")

final_pair_data_rain$mid.year <- as.factor(final_pair_data_rain$mid.year)
final_pair_data_rain$pair.id <- as.character(final_pair_data_rain$pair.id)
final_pair_data_rain$season <- as.factor(final_pair_data_rain$season)

final_pair_data_temp$site1 <- as.factor(final_pair_data_temp$site1)
final_pair_data_temp$site2 <- as.factor(final_pair_data_temp$site2)
final_pair_data_temp$mid.year <- as.factor(final_pair_data_temp$mid.year)
final_pair_data_temp$pair.id <- as.character(final_pair_data_temp$pair.id)
final_pair_data_temp$season <- as.factor(final_pair_data_temp$season)

######## TEMPERATURE SIGNIFICANCE TESTING

###### create 3 new pair_attr files whicih compares early, late and overall
final_pair_1985 <- final_pair_data_temp[final_pair_data_temp$mid.year==1984.5,]
final_pair_2000 <- final_pair_data_temp[final_pair_data_temp$mid.year==1999.5,]
final_pair_2012 <- final_pair_data_temp[final_pair_data_temp$mid.year==2011.5,]

final_pair_early_temp <- rbind(final_pair_1985, final_pair_2000) # comparison of early years 1980 & 1995
final_pair_late_temp <- rbind(final_pair_2000, final_pair_2012) # comparison of late years 1995 & 2007

## early model (85-00) for each season
season_early <- unique(final_pair_data_temp$season)

# 
# logitTransform <- function(p) { log(p/(1-p)) }
# final_pair_data_temp$lag0_logit <- logitTransform(final_pair_data_temp$lag0)
# 
# hist(final_pair_data_temp$lag0)
# hist(final_pair_data_temp$lag0_logit)

results_table_early<-NULL
for (i in season_early){
  print(i)
  
  ## create unique pair_attr for each species
  final_pair_season <- final_pair_early_temp[final_pair_early_temp$season==i,]
    
    # early_model_2 <- lmer(lag0 ~ mid.year + (1|pair.id), data = final_pair_season)
    early_model <- lmer(log(lag0/(1-lag0)) ~ mid.year + (1|pair.id), data = final_pair_season)
    
    # qqnorm(residuals(early_model))
    # qqline(residuals(early_model))
    # plot(early_model)
    # 
    # qqnorm(residuals(early_model_2))
    # qqline(residuals(early_model_2))
    # plot(early_model_2)
    
    
    summary(early_model)
    anova(early_model)
    
    ### save and plot the results ###
    results_table_temp <- data.frame(summary(early_model)$coefficients[,1:5],i)
    results_table_early <-rbind(results_table_early,results_table_temp)
    
  }
## all non-significant
## save table
write.csv(results_table_early, file="../Results/Climate_results/temp_85_00_UKBMS.csv", row.names=FALSE)

## late model (00-12) for each season
season_late <- unique(final_pair_data_temp$season)

results_table_late<-NULL
for (i in season_late){
  print(i)
  
  ## create unique pair_attr for each species
  final_pair_season <- final_pair_late_temp[final_pair_late_temp$season==i,]
  
  late_model <- (lmer(log(lag0/(1-lag0)) ~ mid.year + (1|pair.id), data = final_pair_season))
  summary(late_model)
  anova(late_model)
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(late_model)$coefficients[,1:5],i)
  results_table_late <-rbind(results_table_late,results_table_temp)
  
}
## all non-significant
## save table
write.csv(results_table_late, file="../Results/Climate_results/temp_00_12_UKBMS.csv", row.names=FALSE)


######## RAINFALL SIGNIFICANCE TESTING

###### create 3 new pair_attr files whicih compares early, late and overall
final_pair_1985 <- final_pair_data_rain[final_pair_data_rain$mid.year==1984.5,]
final_pair_2000 <- final_pair_data_rain[final_pair_data_rain$mid.year==1999.5,]
final_pair_2012 <- final_pair_data_rain[final_pair_data_rain$mid.year==2011.5,]

final_pair_early_rain <- rbind(final_pair_1985, final_pair_2000) # comparison of early years 1980 & 1995
final_pair_late_rain <- rbind(final_pair_2000, final_pair_2012) # comparison of late years 1995 & 2007

## early model (85-00) for each season
season_early <- unique(final_pair_data_rain$season)

results_table_early2<-NULL
for (i in season_early){
  print(i)
  
  ## create unique pair_attr for each species
  final_pair_season <- final_pair_early_rain[final_pair_early_rain$season==i,]
  
  early_model <- (lmer(lag0 ~ mid.year + (1|pair.id), data = final_pair_season))
  summary(early_model)
  anova(early_model)
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(early_model)$coefficients[,1:5],i)
  results_table_early2 <-rbind(results_table_early2,results_table_temp)
  
}
## all non-significant
## save table
write.csv(results_table_early2, file="../Results/Climate_results/rain_85_00_UKBMS.csv", row.names=FALSE)

## late model (00-12) for each season
season_late <- unique(final_pair_data_rain$season)

results_table_late2<-NULL
for (i in season_late){
  print(i)
  
  ## create unique pair_attr for each species
  final_pair_season <- final_pair_late_rain[final_pair_late_rain$season==i,]
  
  late_model <- (lmer(lag0 ~ mid.year + (1|pair.id), data = final_pair_season))
  summary(late_model)
  anova(late_model)
  
  ### save and plot the results ###
  results_table_temp <- data.frame(summary(late_model)$coefficients[,1:5],i)
  results_table_late2 <-rbind(results_table_late2,results_table_temp)
  
}
## all non-significant
## save table
write.csv(results_table_late2, file="../Results/Climate_results/rain_00_12_UKBMS.csv", row.names=FALSE)


temp1 <- final_pair_data_temp[(final_pair_data_temp$season=="c" & final_pair_data_temp$mid.year==1984.5),]
hist(temp1$lag0)
temp2 <- final_pair_data_temp[(final_pair_data_temp$season=="c" & final_pair_data_temp$mid.year==1999.5),]
hist(temp2$lag0)
temp3 <- final_pair_data_temp[(final_pair_data_temp$season=="c" & final_pair_data_temp$mid.year==2011.5),]
hist(temp3$lag0)

temp_sum <- final_pair_data_temp[(final_pair_data_temp$season=="c"),]
hist(temp_sum$lag0)
temp_aut <- final_pair_data_temp[(final_pair_data_temp$season=="d"),]
hist(temp_aut$lag0)
temp_win <- final_pair_data_temp[(final_pair_data_temp$season=="a"),]
hist(temp_win$lag0)
temp_spr <- final_pair_data_temp[(final_pair_data_temp$season=="b"),]
hist(temp_spr$lag0)

## number of sites
site1 <- unique(final_pair_data_temp[,1, drop=FALSE])
site2 <- unique(final_pair_data_temp[,2, drop=FALSE])
colnames(site1)[1] <- "site"
colnames(site2)[1] <- "site"
site_list <- rbind(site1, site2)
site_list <- unique(site_list) ## 487 sites

site1 <- unique(final_pair_data_rain[,1, drop=FALSE])
site2 <- unique(final_pair_data_rain[,2, drop=FALSE])
colnames(site1)[1] <- "site"
colnames(site2)[1] <- "site"
site_list2 <- rbind(site1, site2)
site_list2 <- unique(site_list) ## 487 sites
