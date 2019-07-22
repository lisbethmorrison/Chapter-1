###########################################################
## Title: Calculate pop sync for meadow brown sites 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: June 2019
##########################################################

rm(list=ls()) # clear R

## load data
final_data <- read.csv("../Data/Butterfly_sync_data/final_data_all_spp.csv", header = TRUE) # add growth rate data
site_data <- read.csv("../Data/UKBMS_data/meadow_brown_sites.csv", header = TRUE)

## remove all species except meadow brown
final_data<-final_data[(final_data$name==75),]
## merge with sites of interest
site_data <- site_data[,-(1:4)]
site_data <- as.data.frame(site_data)
final_data <- merge(final_data, site_data, by.x="site", by.y="site_data", all.y=TRUE)
final_data<-final_data[final_data$year>=1980,] ## only look at year 1980 onwards

### run synchrony 
## end up with one value of syncrhony for each pair of sites
final_pair_data <- NULL

################################
# create a matrix to be filled #

site.list <- unique(final_data$site)
year.list<-(min(final_data$year):max(final_data$year))
Grow.matrix<-matrix(c(final_data$gr),nrow=length(year.list))
ncol(Grow.matrix)
rownames(Grow.matrix)<-year.list
colnames(Grow.matrix)<-site.list

num.ts <- length(site.list)   # number of time series
TS <- matrix(Grow.matrix, ncol=num.ts)    # stores Grow data in a matrix     (simply removes names!)

# create site match table used to correct sites names after calculating synchrony
site_match <- data.frame(TS_name = 1:ncol(TS), site_name = colnames(Grow.matrix))

###############################
# cross-correlation functions #

pair.list <- t(combn(1:num.ts,2))
nrow(pair.list) # 136 pairs of sites

colnames(pair.list) <- c("site1","site2")

## calculate cross-correlation functions at different time lags
## (NB - we are only be interested in the current year, i.e. lag = 0)
max.lag <- 0 # maximum time lag over which cross-correlations are assessed
CCF <- data.frame(matrix(NA, ncol = 2*max.lag+2, nrow = nrow(pair.list)))	## matrix to store cross-correlation coefficients
names(CCF) <- c(paste("lag", -max.lag:max.lag, sep=""), "numYears") ## adds column names

for(i in 1:nrow(pair.list)){
  print(pair.list[i,])
  ## calculate number of years to base correlation on...
  CCF$numYears[i] <- length(na.omit((TS[,pair.list[i,1]]+TS[,pair.list[i,2]]))) ## add up number of times when sum can be done (i.e. neither is NA) then only include site with > 6 common years of data
  
  ## attempt calculation of CCF...
  try(CCF[i,-ncol(CCF)] <- ccf(TS[,pair.list[i,1]], TS[,pair.list[i,2]], lag.max=max.lag, na.action=na.exclude, plot=F, type="correlation")$acf, silent=T) ## use try to prevent crashes if no data
} ## this will take a long time!!!

#hist(CCF$lag0)
#hist(CCF$numYears)

site.attr<-unique(final_data$site)

pair.attr <- pair.list ## matrix to hold attributes of pairs...
colnames(pair.attr) <- c("site1","site2")
pair.attr <- cbind(pair.attr, CCF)     ### add in correlation scores and number of comparisons at each site.
head(pair.attr)

###  correct site numbers ### 
pair.attr <- merge(pair.attr, site_match, by.x="site1", by.y="TS_name")  
pair.attr <- pair.attr[,c("site_name", "site2", "lag0", "numYears")]
names(pair.attr) <- gsub("site_name", "site1", names(pair.attr))
pair.attr <- merge(pair.attr, site_match, by.x="site2", by.y="TS_name")  
pair.attr <- pair.attr[,c("site1", "site_name", "lag0", "numYears")]
names(pair.attr) <- gsub("site_name", "site2", names(pair.attr))

## save file
write.csv(pair.attr, file="../Data/Butterfly_sync_data/meadow_brown_synchrony.csv", row.names=FALSE)



