###########################################################
## Title: Moving window climate synchrony calculation BBS
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2019
##########################################################

rm(list=ls()) # clear R
library(dplyr)
library(maptools)
library(ggplot2)
options(scipen=999)

## read in data
mean_temp <- read.csv("../Data/MetOffice_data/Mean_temp_1980_2016_final.csv", header=TRUE)
rainfall <- read.csv("../Data/MetOffice_data/Mean_rainfall_1980_2016_final.csv", header=TRUE)
pair_attr_BBS <- read.csv("../Data/Bird_sync_data/pair_attr_BBS.csv", header=TRUE)

##### CALCULATE SEASONAL MEAN TEMPERATURE SYNCHRONY

## find unique BBS pair_attr sites to merge with mean_temp sites
site1 <- unique(subset(pair_attr_BBS[c(2,9,10)]))
site2 <- unique(subset(pair_attr_BBS[c(3,11,12)]))
colnames(site1) <- c("site", "east", "north")
colnames(site2) <- c("site", "east", "north")
site_list <- rbind(site1, site2)
site_list <- unique(site_list) ## 2499 sites

## change coordinates from 1km to 5km
nrow(site_list)   
for (i in 1:length(site_list$east)){
  ifelse(site_list$east[i]%%10000<5000,y<-2500,y<-7500) ## %% is the remainder after division
  site_list$east.5k[i]<-(site_list$east[i]%/%10000)*10000+y ## %/% number of times x fits into y
} ## 544500 changes to 542500
for (i in 1:length(site_list$north)){
  ifelse(site_list$north[i]%%10000<5000,y<-2500,y<-7500)
  site_list$north.5k[i]<-(site_list$north[i]%/%10000)*10000+y
} ## 206900 changes to 207500

site_list2<-site_list[,c(1,4,5)] ## 2499 sites
## remove duplicates where two sites have the same 5km easting and northing
#site_list <- site_list[!duplicated(t(apply(site_list[2:3], 1, sort))),] ## 1861 sites

## merge pair_attr site data with mean_temp (left with only BBS sites)
mean_temp <- merge(mean_temp, site_list2, by.x=c("easting", "northing"), by.y=c("east.5k", "north.5k"))
length(unique(mean_temp$site)) # 2490 sites which match with climate and BBS data 

## merge pair_attr site data with rainfall (left with only BBS sites)
rainfall <- merge(rainfall, site_list2, by.x=c("easting", "northing"), by.y=c("east.5k", "north.5k"))
length(unique(rainfall$site)) # 2490 sites 

## save site info for the 2490 5k sites
sites <- data.frame(unique(mean_temp$site))
names(sites) <- "site"
site_list <- merge(site_list, sites, by="site", all=FALSE)
## save this file
write.csv(site_list, file="../Data/MetOffice_data/site_list_5km_BBS.csv", row.names=FALSE)


##########################################################################################################
################################### CALCULATE SYNCHRONY ##################################################
##########################################################################################################

################## TEMPERATURE ###################

final_pair_data <- NULL
season <- unique(mean_temp$ag)
### split based on season ###
for (g in season){ # loop through each season
  season_data <- mean_temp[mean_temp$ag==g,]
  total_comp <- NULL
  print(paste("season",g))    
  
  year.list<-1994:2007
  
  for (i in year.list){ # loop through year.list
    start.year<-i
    mid.year<-i+4.5
    print(paste("mid.year=",mid.year)) 
    end.year<-i+9
    subset_10yr_data<-season_data[season_data$year>=start.year&season_data$year<=end.year,]
    subset_10yr_data<-subset_10yr_data[order(subset_10yr_data$year),]
    ################################
    # create a matrix to be filled #
    site.list <- unique(subset_10yr_data$site)
    year.list.temp <- min(subset_10yr_data$year):max(subset_10yr_data$year)
    
    Temp.matrix <- reshape2::acast(subset_10yr_data, year ~ site, value.var="mean_temp")
    
    # Temp.matrix<-matrix(c(subset_10yr_data$mean_temp), nrow=length(year.list.temp))
    # ncol(Temp.matrix)
    # rownames(Temp.matrix)<-year.list.temp
    # colnames(Temp.matrix)<-site.list
    
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
    
    ###  correct site numbers ###
    pair.list <- merge(pair.list, site_match, by.x="site1", by.y="TS_name")
    pair.list <- pair.list[,c("site_name", "site2")] ## get rid of site1 (which are not real site names)
    names(pair.list) <- gsub("site_name", "site1", names(pair.list))
    pair.list <- merge(pair.list, site_match, by.x="site2", by.y="TS_name")
    pair.list <- pair.list[,c("site1", "site_name")] ## get rid of site2 (which are not real site names)
    names(pair.list) <- gsub("site_name", "site2", names(pair.list)) ## now site1 and site2 are REAL site names
    
    ## merge pair.list with site data
    pair_attr_BBS2 <- unique(pair_attr_BBS[,c(2,3,7)]) ## site1, site2 and start year
    pair_attr_BBS2[] <- lapply(pair_attr_BBS2, factor)
    pair_attr_BBS2 <- pair_attr_BBS2[pair_attr_BBS2$start.year==i,] ## subset by start year i
    site_data_pair1 <- merge(pair.list, pair_attr_BBS2, by=c("site1", "site2"), all=FALSE)
    pair_attr_BBS_reverse <- pair_attr_BBS2
    names(pair_attr_BBS_reverse)[1:3] <- c("site2", "site1", "start.year")
    site_data_pair2 <- merge(pair.list, pair_attr_BBS_reverse, by=c("site1", "site2"), all=FALSE)
    site_data_pair <- rbind(site_data_pair1, site_data_pair2)
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
    final_pair_data <- rbind(final_pair_data, all_pair_attr)
    
  } ## end in year.list
} ## end in season

head(final_pair_data)
## save data
write.csv(final_pair_data, file="../Data/MetOffice_data/final_pair_data_mean_temp_BBS.csv", row.names=FALSE)

################## RAINFALL ###################

#### calculate syncrhony
final_pair_data <- NULL
season <- unique(rainfall$ag)
### split based on season ###
for (g in season){ # loop through each season
  season_data <- rainfall[rainfall$ag==g,]
  total_comp <- NULL
  print(paste("season",g))    
  
  year.list<-1994:2007
  
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
    
    Rain.matrix <- reshape2::acast(subset_10yr_data, year ~ site, value.var="mean_rainfall")
    
    # Rain.matrix<-matrix(c(subset_10yr_data$mean_rainfall), nrow=length(year.list.temp))
    # ncol(Rain.matrix)
    # rownames(Rain.matrix)<-year.list.temp
    # colnames(Rain.matrix)<-site.list
    
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
    pair_attr_BBS2 <- unique(pair_attr_BBS[,c(2,3,7)]) ## site1, site2 and start year
    pair_attr_BBS2[] <- lapply(pair_attr_BBS2, factor)
    pair_attr_BBS2 <- pair_attr_BBS2[pair_attr_BBS2$start.year==i,] ## subset by start year i
    site_data_pair1 <- merge(pair.list, pair_attr_BBS2, by=c("site1", "site2"), all=FALSE)
    pair_attr_BBS_reverse <- pair_attr_BBS2
    names(pair_attr_BBS_reverse)[1:3] <- c("site2", "site1", "start.year")
    site_data_pair2 <- merge(pair.list, pair_attr_BBS_reverse, by=c("site1", "site2"), all=FALSE)
    site_data_pair <- rbind(site_data_pair1, site_data_pair2)
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
    final_pair_data <- rbind(final_pair_data, all_pair_attr)
    
  } ## end in year.list
} ## end in season

head(final_pair_data)

## save file
write.csv(final_pair_data, file="../Data/MetOffice_data/final_pair_data_mean_rainfall_BBS.csv", row.names=FALSE)

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

# rm(list=ls()) # clear R
# 
# ## read in final pair data for temp and rainfall
# final_pair_data_rain <- read.csv("../Data/MetOffice_data/final_pair_data_mean_rainfall_BBS.csv", header=TRUE)
# final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp_BBS.csv", header=TRUE)  
# 
# summ_data_temp <- final_pair_data_temp %>% 
#   group_by(mid.year,season) %>% 
#   summarise(freq = n()) ## exact same number pairs of sites for each season and moving window
# 
# summ_data_rain <- final_pair_data_rain %>% 
#   group_by(mid.year,season) %>% 
#   summarise(freq = n()) ## exact same number pairs of sites for each season and moving window
# 
# 
# ######################
# ###### Rainfall ######
# ######################
# 
# ## create PairID column
# final_pair_data_rain$pair.id <- paste("ID", final_pair_data_rain$site1, final_pair_data_rain$site2, sep = "_")
# 
# final_pair_data_rain$mid.year <- as.factor(final_pair_data_rain$mid.year)
# final_pair_data_rain$pair.id <- as.character(final_pair_data_rain$pair.id)
# final_pair_data_rain$season <- as.factor(final_pair_data_rain$season)
# 
# ## split into 4 dataframes (one for each season)
# winter_rainfall <- final_pair_data_rain[final_pair_data_rain$season=="a",]
# winter_rainfall <- droplevels(winter_rainfall)
# spring_rainfall <- final_pair_data_rain[final_pair_data_rain$season=="b",]
# spring_rainfall <- droplevels(spring_rainfall)
# summer_rainfall <- final_pair_data_rain[final_pair_data_rain$season=="c",]
# summer_rainfall <- droplevels(summer_rainfall)
# autumn_rainfall <- final_pair_data_rain[final_pair_data_rain$season=="d",]
# autumn_rainfall <- droplevels(autumn_rainfall)
# 
# ## run mixed effects model to extract coefficients for each dataframe
# library(lme4)
# 
# ###### WINTER RAINFALL ######
# winter_rainfall_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=winter_rainfall)
# summary(winter_rainfall_model)
# ## save results
# results_table_winter_rain <- data.frame(summary(winter_rainfall_model)$coefficients[,1:3])
# ## change names and add in parameter column ##
# names(results_table_winter_rain) <- c("synchrony", "SD", "t")
# results_table_winter_rain$parameter <- paste(row.names(results_table_winter_rain))
# rownames(results_table_winter_rain) <- 1:nrow(results_table_winter_rain)
# ## change parameter names to year
# results_table_winter_rain$parameter <- rep(1999:2012)
# 
# ### rescale estimate, SD and CI ### 
# ## add 10 to each value so scale remains the same, but values are all positive
# results_table_winter_rain$new_sync <- results_table_winter_rain$synchrony + 10 
# ## scale using new_sync
# results_table_winter_rain$rescaled_sync <- results_table_winter_rain$new_sync*(100/results_table_winter_rain$new_sync[1]) ## rescale to 100
# results_table_winter_rain$rescaled_sd <- results_table_winter_rain$SD*(100/results_table_winter_rain$new_sync[1])
# results_table_winter_rain$rescaled_ci <- results_table_winter_rain$rescaled_sd*1.96
# ## remove new_sync column
# results_table_winter_rain <- subset(results_table_winter_rain, select = -c(new_sync))
# ## save final results table ##
# write.csv(results_table_winter_rain, file = "../Results/Climate_results/winter_rainfall_synchrony_BBS.csv", row.names=FALSE)
# 
# ###### SPRING RAINFALL ######
# spring_rainfall_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=spring_rainfall)
# summary(spring_rainfall_model)
# ## save results
# results_table_spring_rain <- data.frame(summary(spring_rainfall_model)$coefficients[,1:3])
# ## change names and add in parameter column ##
# names(results_table_spring_rain) <- c("synchrony", "SD", "t")
# results_table_spring_rain$parameter <- paste(row.names(results_table_spring_rain))
# rownames(results_table_spring_rain) <- 1:nrow(results_table_spring_rain)
# ## change parameter names to year
# results_table_spring_rain$parameter <- rep(1999:2012)
# 
# ### rescale estimate, SD and CI ### 
# ## add 10 to each value so scale remains the same, but values are all positive
# results_table_spring_rain$new_sync <- results_table_spring_rain$synchrony + 10 
# ## scale using new_sync
# results_table_spring_rain$rescaled_sync <- results_table_spring_rain$new_sync*(100/results_table_spring_rain$new_sync[1]) ## rescale to 100
# results_table_spring_rain$rescaled_sd <- results_table_spring_rain$SD*(100/results_table_spring_rain$new_sync[1])
# results_table_spring_rain$rescaled_ci <- results_table_spring_rain$rescaled_sd*1.96
# ## remove new_sync column
# results_table_spring_rain <- subset(results_table_spring_rain, select = -c(new_sync))
# ## save final results table ##
# write.csv(results_table_spring_rain, file = "../Results/Climate_results/spring_rainfall_synchrony_BBS.csv", row.names=FALSE)
# 
# ###### SUMMER RAINFALL ######
# summer_rainfall_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=summer_rainfall)
# summary(summer_rainfall_model)
# ## save results
# results_table_summer_rain <- data.frame(summary(summer_rainfall_model)$coefficients[,1:3])
# ## change names and add in parameter column ##
# names(results_table_summer_rain) <- c("synchrony", "SD", "t")
# results_table_summer_rain$parameter <- paste(row.names(results_table_summer_rain))
# rownames(results_table_summer_rain) <- 1:nrow(results_table_summer_rain)
# ## change parameter names to year
# results_table_summer_rain$parameter <- rep(1999:2012)
# 
# ### rescale estimate, SD and CI ### 
# ## add 10 to each value so scale remains the same, but values are all positive
# results_table_summer_rain$new_sync <- results_table_summer_rain$synchrony + 10 
# ## scale using new_sync
# results_table_summer_rain$rescaled_sync <- results_table_summer_rain$new_sync*(100/results_table_summer_rain$new_sync[1]) ## rescale to 100
# results_table_summer_rain$rescaled_sd <- results_table_summer_rain$SD*(100/results_table_summer_rain$new_sync[1])
# results_table_summer_rain$rescaled_ci <- results_table_summer_rain$rescaled_sd*1.96
# ## remove new_sync column
# results_table_summer_rain <- subset(results_table_summer_rain, select = -c(new_sync))
# ## save final results table ##
# write.csv(results_table_summer_rain, file = "../Results/Climate_results/summer_rainfall_synchrony_BBS.csv", row.names=FALSE)
# 
# ###### AUTUMN RAINFALL ######
# autumn_rainfall_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=autumn_rainfall)
# summary(autumn_rainfall_model)
# plot(autumn_rainfall_model)
# qqnorm(residuals(autumn_rainfall_model))
# ## save results
# results_table_autumn_rain <- data.frame(summary(autumn_rainfall_model)$coefficients[,1:3])
# ## change names and add in parameter column ##
# names(results_table_autumn_rain) <- c("synchrony", "SD", "t")
# results_table_autumn_rain$parameter <- paste(row.names(results_table_autumn_rain))
# rownames(results_table_autumn_rain) <- 1:nrow(results_table_autumn_rain)
# ## change parameter names to year
# results_table_autumn_rain$parameter <- rep(1999:2012)
# 
# ### rescale estimate, SD and CI ### 
# ## add 10 to each value so scale remains the same, but values are all positive
# results_table_autumn_rain$new_sync <- results_table_autumn_rain$synchrony + 10 
# ## scale using new_sync
# results_table_autumn_rain$rescaled_sync <- results_table_autumn_rain$new_sync*(100/results_table_autumn_rain$new_sync[1]) ## rescale to 100
# results_table_autumn_rain$rescaled_sd <- results_table_autumn_rain$SD*(100/results_table_autumn_rain$new_sync[1])
# results_table_autumn_rain$rescaled_ci <- results_table_autumn_rain$rescaled_sd*1.96
# ## remove new_sync column
# results_table_autumn_rain <- subset(results_table_autumn_rain, select = -c(new_sync))
# ## save final results table ##
# write.csv(results_table_autumn_rain, file = "../Results/Climate_results/autumn_rainfall_synchrony_BBS.csv", row.names=FALSE)
# 
# ##### plots ######
# ## load data
# winter_rain <- read.csv("../Results/Climate_results/winter_rainfall_synchrony_BBS.csv", header=TRUE)
# spring_rain <- read.csv("../Results/Climate_results/spring_rainfall_synchrony_BBS.csv", header=TRUE)
# summer_rain <- read.csv("../Results/Climate_results/summer_rainfall_synchrony_BBS.csv", header=TRUE)
# autumn_rain <- read.csv("../Results/Climate_results/autumn_rainfall_synchrony_BBS.csv", header=TRUE)
# 
# winter_rain_plot <- ggplot(winter_rain, aes(x = parameter, y = rescaled_sync)) +
#   stat_smooth(colour="black", method=loess, se=FALSE) +
#   geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
#   geom_point(size=2) + 
#   labs(x = "Mid-year of moving window", y = "Rainfall synchrony") +
#   #scale_y_continuous(breaks=seq(40,160,10)) +
#   scale_x_continuous(breaks=seq(1999,2012,3)) +
#   geom_hline(yintercept = 100, linetype = "dashed") +
#   theme_bw() +
#   theme(text = element_text(size = 16)) +
#   labs(size=3) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# winter_rain_plot
# 
# spring_rain_plot <- ggplot(spring_rain, aes(x = parameter, y = rescaled_sync)) +
#   stat_smooth(colour="black", method=loess, se=FALSE) +
#   geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
#   geom_point(size=2) + 
#   labs(x = "Mid-year of moving window", y = "Rainfall synchrony") +
#   #scale_y_continuous(breaks=seq(40,160,10)) +
#   scale_x_continuous(breaks=seq(1999,2012,3)) +
#   geom_hline(yintercept = 100, linetype = "dashed") +
#   theme_bw() +
#   theme(text = element_text(size = 16)) +
#   labs(size=3) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# spring_rain_plot
# 
# summer_rain_plot <- ggplot(summer_rain, aes(x = parameter, y = rescaled_sync)) +
#   stat_smooth(colour="black", method=loess, se=FALSE) +
#   geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
#   geom_point(size=2) + 
#   labs(x = "Mid-year of moving window", y = "Rainfall synchrony") +
#   #scale_y_continuous(breaks=seq(40,160,10)) +
#   scale_x_continuous(breaks=seq(1999,2012,3)) +
#   geom_hline(yintercept = 100, linetype = "dashed") +
#   theme_bw() +
#   theme(text = element_text(size = 16)) +
#   labs(size=3) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# summer_rain_plot
# 
# autumn_rain_plot <- ggplot(autumn_rain, aes(x = parameter, y = rescaled_sync)) +
#   stat_smooth(colour="black", method=loess, se=FALSE) +
#   geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
#   geom_point(size=2) + 
#   labs(x = "Mid-year of moving window", y = "Rainfall synchrony") +
#   #scale_y_continuous(breaks=seq(40,160,10)) +
#   scale_x_continuous(breaks=seq(1999,2012,3)) +
#   geom_hline(yintercept = 100, linetype = "dashed") +
#   theme_bw() +
#   theme(text = element_text(size = 16)) +
#   labs(size=3) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# autumn_rain_plot
# 
# ## save all 4 rainfall plots
# library(ggpubr)
# rain_plots <- ggarrange(winter_rain_plot, spring_rain_plot, summer_rain_plot, autumn_rain_plot,
#                         hjust = 0, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)"))
# rain_plots
# ggsave("../Graphs/Climate/seasonal_rainfall_synchrony_BBS.png", plot = rain_plots, width=12, height=10)
# 
# #########################
# ###### Temperature ######
# #########################
# 
# ## create PairID column
# final_pair_data_temp$pair.id <- paste("ID", final_pair_data_temp$site1, final_pair_data_temp$site2, sep = "_")
# 
# final_pair_data_temp$mid.year <- as.factor(final_pair_data_temp$mid.year)
# final_pair_data_temp$pair.id <- as.character(final_pair_data_temp$pair.id)
# final_pair_data_temp$season <- as.factor(final_pair_data_temp$season)
# 
# ## split into 4 dataframes (one for each season)
# winter_temp <- final_pair_data_temp[final_pair_data_temp$season=="a",]
# winter_temp <- droplevels(winter_temp)
# spring_temp <- final_pair_data_temp[final_pair_data_temp$season=="b",]
# spring_temp <- droplevels(spring_temp)
# summer_temp <- final_pair_data_temp[final_pair_data_temp$season=="c",]
# summer_temp <- droplevels(summer_temp)
# autumn_temp <- final_pair_data_temp[final_pair_data_temp$season=="d",]
# autumn_temp <- droplevels(autumn_temp)
# 
# ## run mixed effects model to extract coefficients for each dataframe
# library(lme4)
# 
# ###### WINTER TEMPERATURE ######
# winter_temp_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=winter_temp)
# summary(winter_temp_model)
# ## save results
# results_table_winter_temp <- data.frame(summary(winter_temp_model)$coefficients[,1:3])
# ## change names and add in parameter column ##
# names(results_table_winter_temp) <- c("synchrony", "SD", "t")
# results_table_winter_temp$parameter <- paste(row.names(results_table_winter_temp))
# rownames(results_table_winter_temp) <- 1:nrow(results_table_winter_temp)
# ## change parameter names to year
# results_table_winter_temp$parameter <- rep(1999:2012)
# 
# ### rescale estimate, SD and CI ### 
# ## add 10 to each value so scale remains the same, but values are all positive
# results_table_winter_temp$new_sync <- results_table_winter_temp$synchrony + 10 
# ## scale using new_sync
# results_table_winter_temp$rescaled_sync <- results_table_winter_temp$new_sync*(100/results_table_winter_temp$new_sync[1]) ## rescale to 100
# results_table_winter_temp$rescaled_sd <- results_table_winter_temp$SD*(100/results_table_winter_temp$new_sync[1])
# results_table_winter_temp$rescaled_ci <- results_table_winter_temp$rescaled_sd*1.96
# ## remove new_sync column
# results_table_winter_temp <- subset(results_table_winter_temp, select = -c(new_sync))
# ## save final results table ##
# write.csv(results_table_winter_temp, file = "../Results/Climate_results/winter_temp_synchrony_BBS.csv", row.names=FALSE)
# 
# ###### SPRING TEMPERATURE ######
# spring_temp_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=spring_temp)
# summary(spring_temp_model)
# ## save results
# results_table_spring_temp <- data.frame(summary(spring_temp_model)$coefficients[,1:3])
# ## change names and add in parameter column ##
# names(results_table_spring_temp) <- c("synchrony", "SD", "t")
# results_table_spring_temp$parameter <- paste(row.names(results_table_spring_temp))
# rownames(results_table_spring_temp) <- 1:nrow(results_table_spring_temp)
# ## change parameter names to year
# results_table_spring_temp$parameter <- rep(1999:2012)
# 
# ### rescale estimate, SD and CI ### 
# ## add 10 to each value so scale remains the same, but values are all positive
# results_table_spring_temp$new_sync <- results_table_spring_temp$synchrony + 10 
# ## scale using new_sync
# results_table_spring_temp$rescaled_sync <- results_table_spring_temp$new_sync*(100/results_table_spring_temp$new_sync[1]) ## rescale to 100
# results_table_spring_temp$rescaled_sd <- results_table_spring_temp$SD*(100/results_table_spring_temp$new_sync[1])
# results_table_spring_temp$rescaled_ci <- results_table_spring_temp$rescaled_sd*1.96
# ## remove new_sync column
# results_table_spring_temp <- subset(results_table_spring_temp, select = -c(new_sync))
# ## save final results table ##
# write.csv(results_table_spring_temp, file = "../Results/Climate_results/spring_temp_synchrony_BBS.csv", row.names=FALSE)
# 
# ###### SUMMER TEMPERATURE ######
# summer_temp_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=summer_temp)
# summary(summer_temp_model)
# ## save results
# results_table_summer_temp <- data.frame(summary(summer_temp_model)$coefficients[,1:3])
# ## change names and add in parameter column ##
# names(results_table_summer_temp) <- c("synchrony", "SD", "t")
# results_table_summer_temp$parameter <- paste(row.names(results_table_summer_temp))
# rownames(results_table_summer_temp) <- 1:nrow(results_table_summer_temp)
# ## change parameter names to year
# results_table_summer_temp$parameter <- rep(1999:2012)
# 
# ### rescale estimate, SD and CI ### 
# ## add 10 to each value so scale remains the same, but values are all positive
# results_table_summer_temp$new_sync <- results_table_summer_temp$synchrony + 10 
# ## scale using new_sync
# results_table_summer_temp$rescaled_sync <- results_table_summer_temp$new_sync*(100/results_table_summer_temp$new_sync[1]) ## rescale to 100
# results_table_summer_temp$rescaled_sd <- results_table_summer_temp$SD*(100/results_table_summer_temp$new_sync[1])
# results_table_summer_temp$rescaled_ci <- results_table_summer_temp$rescaled_sd*1.96
# ## remove new_sync column
# results_table_summer_temp <- subset(results_table_summer_temp, select = -c(new_sync))
# ## save final results table ##
# write.csv(results_table_summer_temp, file = "../Results/Climate_results/summer_temp_synchrony_BBS.csv", row.names=FALSE)
# 
# ###### AUTUMN TEMPERATURE ######
# autumn_temp_model <- lmer(lag0 ~ mid.year + (1|pair.id)-1, data=autumn_temp)
# summary(autumn_temp_model)
# ## save results
# results_table_autumn_temp <- data.frame(summary(autumn_temp_model)$coefficients[,1:3])
# ## change names and add in parameter column ##
# names(results_table_autumn_temp) <- c("synchrony", "SD", "t")
# results_table_autumn_temp$parameter <- paste(row.names(results_table_autumn_temp))
# rownames(results_table_autumn_temp) <- 1:nrow(results_table_autumn_temp)
# ## change parameter names to year
# results_table_autumn_temp$parameter <- rep(1999:2012)
# 
# ### rescale estimate, SD and CI ### 
# ## add 10 to each value so scale remains the same, but values are all positive
# results_table_autumn_temp$new_sync <- results_table_autumn_temp$synchrony + 10 
# ## scale using new_sync
# results_table_autumn_temp$rescaled_sync <- results_table_autumn_temp$new_sync*(100/results_table_autumn_temp$new_sync[1]) ## rescale to 100
# results_table_autumn_temp$rescaled_sd <- results_table_autumn_temp$SD*(100/results_table_autumn_temp$new_sync[1])
# results_table_autumn_temp$rescaled_ci <- results_table_autumn_temp$rescaled_sd*1.96
# ## remove new_sync column
# results_table_autumn_temp <- subset(results_table_autumn_temp, select = -c(new_sync))
# ## save final results table ##
# write.csv(results_table_autumn_temp, file = "../Results/Climate_results/autumn_temp_synchrony_BBS.csv", row.names=FALSE)
# 
# ##### plots ######
# ## load data
# 
# winter_temp <- read.csv("../Results/Climate_results/winter_temp_synchrony_BBS.csv", header=TRUE)
# spring_temp <- read.csv("../Results/Climate_results/spring_temp_synchrony_BBS.csv", header=TRUE)
# summer_temp <- read.csv("../Results/Climate_results/summer_temp_synchrony_BBS.csv", header=TRUE)
# autumn_temp <- read.csv("../Results/Climate_results/autumn_temp_synchrony_BBS.csv", header=TRUE)
# 
# winter_temp_plot <- ggplot(winter_temp, aes(x = parameter, y = rescaled_sync)) +
#   stat_smooth(colour="black", method=loess, se=FALSE) +
#   geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
#   geom_point(size=2) + 
#   labs(x = "Mid-year of moving window", y = "Temperature synchrony") +
#   #scale_y_continuous(breaks=seq(40,160,10)) +
#   scale_x_continuous(breaks=seq(1999,2012,3)) +
#   geom_hline(yintercept = 100, linetype = "dashed") +
#   theme_bw() +
#   theme(text = element_text(size = 16)) +
#   labs(size=3) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# winter_temp_plot
# 
# spring_temp_plot <- ggplot(spring_temp, aes(x = parameter, y = rescaled_sync)) +
#   stat_smooth(colour="black", method=loess, se=FALSE) +
#   geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
#   geom_point(size=2) + 
#   labs(x = "Mid-year of moving window", y = "Temperature synchrony") +
#   #scale_y_continuous(breaks=seq(40,160,10)) +
#   scale_x_continuous(breaks=seq(1999,2012,3)) +
#   geom_hline(yintercept = 100, linetype = "dashed") +
#   theme_bw() +
#   theme(text = element_text(size = 16)) +
#   labs(size=3) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# spring_temp_plot
# 
# summer_temp_plot <- ggplot(summer_temp, aes(x = parameter, y = rescaled_sync)) +
#   stat_smooth(colour="black", method=loess, se=FALSE) +
#   geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
#   geom_point(size=2) + 
#   labs(x = "Mid-year of moving window", y = "Temperature synchrony") +
#   #scale_y_continuous(breaks=seq(40,160,10)) +
#   scale_x_continuous(breaks=seq(1999,2012,3)) +
#   geom_hline(yintercept = 100, linetype = "dashed") +
#   theme_bw() +
#   theme(text = element_text(size = 16)) +
#   labs(size=3) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# summer_temp_plot
# 
# autumn_temp_plot <- ggplot(autumn_temp, aes(x = parameter, y = rescaled_sync)) +
#   stat_smooth(colour="black", method=loess, se=FALSE) +
#   geom_errorbar(aes(ymin = rescaled_sync - rescaled_sd, ymax = rescaled_sync + rescaled_sd), width=0.2, size = 0.5) +
#   geom_point(size=2) + 
#   labs(x = "Mid-year of moving window", y = "Temperature synchrony") +
#   #scale_y_continuous(breaks=seq(40,160,10)) +
#   scale_x_continuous(breaks=seq(1999,2012,3)) +
#   geom_hline(yintercept = 100, linetype = "dashed") +
#   theme_bw() +
#   theme(text = element_text(size = 16)) +
#   labs(size=3) +
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# autumn_temp_plot
# 
# ## save all 4 temperature plots
# library(ggpubr)
# temp_plots <- ggarrange(winter_temp_plot, spring_temp_plot, summer_temp_plot, autumn_temp_plot,
#                         hjust = 0, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)"))
# temp_plots
# ggsave("../Graphs/Climate/seasonal_temperature_synchrony_BBS.png", plot = temp_plots, width=12, height=10)
# 
# ##############################################################################################################
# ######################################## SIGNIFICANCE TESTING ################################################
# ##############################################################################################################
# 
# rm(list=ls()) # clear R
# library(lme4)
# library(lmerTest)
# 
# ## read in data
# final_pair_data_rain <- read.csv("../Data/MetOffice_data/final_pair_data_mean_rainfall_BBS.csv", header=TRUE)
# final_pair_data_temp <- read.csv("../Data/MetOffice_data/final_pair_data_mean_temp_BBS.csv", header=TRUE)  
# 
# summ_data_temp <- final_pair_data_temp %>% 
#   group_by(mid.year,season) %>% 
#   summarise(freq = n())
# 
# summ_data_rain <- final_pair_data_rain %>% 
#   group_by(mid.year,season) %>% 
#   summarise(freq = n()) 
# 
# ## create PairID column
# final_pair_data_rain$pair.id <- paste("ID", final_pair_data_rain$site1, final_pair_data_rain$site2, sep = "_")
# final_pair_data_temp$pair.id <- paste("ID", final_pair_data_temp$site1, final_pair_data_temp$site2, sep = "_")
# 
# final_pair_data_rain$mid.year <- as.factor(final_pair_data_rain$mid.year)
# final_pair_data_rain$pair.id <- as.character(final_pair_data_rain$pair.id)
# final_pair_data_rain$season <- as.factor(final_pair_data_rain$season)
# 
# final_pair_data_temp$mid.year <- as.factor(final_pair_data_temp$mid.year)
# final_pair_data_temp$pair.id <- as.character(final_pair_data_temp$pair.id)
# final_pair_data_temp$season <- as.factor(final_pair_data_temp$season)
# 
# ######## TEMPERATURE SIGNIFICANCE TESTING
# 
# ###### create 3 new pair_attr files whicih compares early, late and overall
# final_pair_1999 <- final_pair_data_temp[final_pair_data_temp$mid.year==1998.5,]
# final_pair_2012 <- final_pair_data_temp[final_pair_data_temp$mid.year==2011.5,]
# 
# final_pair_BBS <- rbind(final_pair_1999, final_pair_2012) # comparison of early years 1985 and 1996
# 
# ## 85-96 for each season
# season <- unique(final_pair_data_temp$season)
# 
# results_table_bbs<-NULL
# for (i in season){
#   print(i)
#   
#   ## create unique pair_attr for each species
#   final_pair_season <- final_pair_BBS[final_pair_BBS$season==i,]
#   
#   bbs_model <- (lmer(lag0 ~ mid.year + (1|pair.id), REML=FALSE, data = final_pair_season))
#   summary(bbs_model)
#   anova(bbs_model)
#   
#   ### save and plot the results ###
#   results_table_temp <- data.frame(summary(bbs_model)$coefficients[,1:5],i)
#   results_table_bbs <-rbind(results_table_bbs,results_table_temp)
#   
# }
# ## seasonal temperature significantly increases in synchrony betweeen 1999 and 2012
# ## save table
# write.csv(results_table_bbs, file="../Results/Climate_results/temp_99_12_BBS.csv", row.names=FALSE)
# 
# 
# ######## RAINFALL SIGNIFICANCE TESTING
# 
# ###### create 2 new pair_attr files whicih compares 99-12
# final_pair_1999 <- final_pair_data_rain[final_pair_data_rain$mid.year==1998.5,]
# final_pair_2012 <- final_pair_data_rain[final_pair_data_rain$mid.year==2011.5,]
# 
# final_pair_BBS <- rbind(final_pair_1999, final_pair_2012) # comparison of early years 1980 & 1996
# 
# ## model for each season
# season <- unique(final_pair_data_rain$season)
# 
# results_table_bbs2<-NULL
# for (i in season){
#   print(i)
#   
#   ## create unique pair_attr for each species
#   final_pair_season <- final_pair_BBS[final_pair_BBS$season==i,]
#   
#   bbs_model2 <- (lmer(lag0 ~ mid.year + (1|pair.id), REML=FALSE, data = final_pair_season))
#   summary(bbs_model2)
#   anova(bbs_model2)
#   
#   ### save and plot the results ###
#   results_table_temp <- data.frame(summary(bbs_model2)$coefficients[,1:5],i)
#   results_table_bbs2 <-rbind(results_table_bbs2,results_table_temp)
#   
# }
# ## all non-significant
# ## save table
# write.csv(results_table_bbs2, file="../Results/Climate_results/rain_99_12_BBS.csv", row.names=FALSE)
# 
# 
# ## number of sites
# site1 <- unique(final_pair_data_temp[,1, drop=FALSE])
# site2 <- unique(final_pair_data_temp[,2, drop=FALSE])
# colnames(site1)[1] <- "site"
# colnames(site2)[1] <- "site"
# site_list <- rbind(site1, site2)
# site_list <- unique(site_list) ## 2490 sites
# 
# site1 <- unique(final_pair_data_rain[,1, drop=FALSE])
# site2 <- unique(final_pair_data_rain[,2, drop=FALSE])
# colnames(site1)[1] <- "site"
# colnames(site2)[1] <- "site"
# site_list <- rbind(site1, site2)
# site_list <- unique(site_list) ## 2490 sites
