###########################################################
## Title: Calculate BTO site data info  
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: November 2017
##########################################################

rm(list=ls()) # clear R

load("Bird_scripts/gr_let2num.Rdata") # load function to convert gridrefs to easting and northing

####################################
#### add woodland bird cbc data ####
####################################

woodland_cbc <- read.csv("../Data/Bird_sync_data/cbc_woodland_birds.csv", header=TRUE)

head(woodland_cbc)
summary(woodland_cbc)

##################################################
#### change grid refs to easting and northing ####
##################################################

#### use function "gr_let2num" from R script "CONVERT_gridrefs_to_east_norths_Colin_H" ####
easting_northing_values <- gr_let2num(woodland_cbc$Gridref)
woodland_cbc <- cbind(woodland_cbc, easting_northing_values)

#### create unique site list #####
site_data <- woodland_cbc[,c("Pcode", "Gridref", "EASTING", "NORTHING")] # create temp dataframe with only site number, easting and northing
site_data <- unique(site_data) # take out duplicated sites, left with 411 sites
names(site_data)[1:4] <- c("Site_code", "Gridref", "Easting", "Northing")

#### save site list info ####
write.csv(site_data, file = "../Data/BTO_data/site_list_east_north.csv", row.names=FALSE)

## save woodland_cbc data which has 411 sites
write.csv(woodland_cbc, file = "../Data/Bird_sync_data/cbc_woodland_birds.csv", row.names=FALSE)

######################################################
#### create list of all pairwise site comparisons ####
######################################################
site.list <- unique(as.character(site_data$Site_code))
num.ts <- length(site.list)

pair.list <- NULL

pair.list <- t(combn(site.list,2))
pair_attr <- data.frame(pair.list)
nrow(pair_attr) ## 84255 rows 
names(pair_attr) <- c("site_a", "site_b")

## add xy data for site_a and site_b
pair_attr <- merge(pair_attr, site_data, by.x="site_a", by.y="Site_code")
names(pair_attr)[4:5] <- c("site_a_EAST", "site_a_NORTH") # rename them
pair_attr <- merge(pair_attr, site_data, by.x="site_b", by.y="Site_code")
names(pair_attr)[7:8] <- c("site_b_EAST", "site_b_NORTH") # rename them
## remove grid ref columns
pair_attr <- pair_attr[-c(3,6)]
# sort out column order #
pair_attr <- pair_attr[,c("site_a", "site_b", "site_a_EAST", "site_a_NORTH", "site_b_EAST", "site_b_NORTH")]

## Add habitat code info ##
## create site code vector which has site number and habitat code ##
site_codes <- woodland_cbc[,-c(1,3:5,7:12)]
site_codes <- unique(site_codes) ## 411 rows

pair_attr <- merge(pair_attr, site_codes, by.x="site_a", by.y="Pcode") # merge by site_a
names(pair_attr)[7] <- "site_a_hab_code" # rename the column
pair_attr <- merge(pair_attr, site_codes, by.x="site_b", by.y="Pcode") # merge by site_b
names(pair_attr)[8] <- "site_b_hab_code" # rename the column

## sort column order
pair_attr <- pair_attr[,c("site_a", "site_b", "site_a_EAST", "site_a_NORTH", "site_b_EAST", "site_b_NORTH", "site_a_hab_code", "site_b_hab_code")]

#### create binary habitat similarity column ####
pair_attr$hab_sim <- ifelse(pair_attr$site_a_hab_code==pair_attr$site_b_hab_code, 1, 0) # matching habitat type = 1, not matching = 0

# take average northing between the pairwise comparisons #
pair_attr$mean_northing <-rowMeans(subset(pair_attr, select = c(site_a_NORTH, site_b_NORTH)), na.rm = TRUE)

### Estimate distance between sites ###
pair_attr$distance <- sqrt((pair_attr$site_a_NORTH-pair_attr$site_b_NORTH)^2 + (pair_attr$site_a_EAST-pair_attr$site_b_EAST)^2) ## pythagoras equation
pair_attr$distance <- pair_attr$distance/1000
## this gives distance in km

summary(pair_attr)

### save pair_attr ###
write.csv(pair_attr, file="../Data/BTO_data/pair_attr_mean_north_dist_hab_sim_CBC.csv", row.names=FALSE)
